//! Rust port of [KSW2](https://github.com/lh3/ksw2), a library for pairwise
//! biological sequence alignment based on dynamic programming.
//!
//! The port keeps a 1:1 correspondence with the C sources that live in the
//! `ksw2/` directory of this repository. The crate is split into two layers:
//!
//! * the idiomatic public API below (this documentation);
//! * the literal line-by-line ports of the C sources, which live in the
//!   private modules listed here and keep the C-shaped signatures (length
//!   arguments, `-1` sentinels, raw flag bits) so they can be diffed against
//!   the C files:
//!
//! | C source           | Rust module (private) | public function  |
//! |--------------------|-----------------------|------------------|
//! | `ksw2_gg.c`        | `gg`                  | [`gg`]           |
//! | `ksw2_gg2.c`       | `gg2`                 | [`gg2`]          |
//! | `ksw2_gg2_sse.c`   | `sse::gg2_sse`        | [`gg2_sse`]      |
//! | `ksw2_extz.c`      | `extz`                | [`extz`]         |
//! | `ksw2_extz2_sse.c` | `sse::extz2_sse`      | [`extz2_sse`]    |
//! | `ksw2_extd.c`      | `extd`                | [`extd`]         |
//! | `ksw2_extd2_sse.c` | `sse::extd2_sse`      | [`extd2_sse`]    |
//! | `ksw2_extf2_sse.c` | `sse::extf2_sse`      | [`extf2_sse`]    |
//! | `ksw2_exts2_sse.c` | `sse::exts2_sse`      | [`exts2_sse`]    |
//!
//! The `*_sse` entry points are literal ports of the SSE4.1 code paths of the
//! C sources (the ones selected by `-march=native` on any x86 CPU released
//! after 2008); they produce bit-identical scores, coordinates and CIGARs.
//! They require the `sse4.1` CPU feature at run time and are only compiled on
//! `x86_64`, exactly like the `#ifdef __SSE2__` guards of the C sources. On
//! other targets use the scalar entry points.
//!
//! Sequences are passed as slices of encoded symbols (`&[u8]`), so the length
//! arguments of the C API have no Rust equivalent — a slice carries its own
//! length. Scoring parameters travel together in a [`Scoring`]; optional
//! quantities (`band_width`, `z_drop`, `x_drop`) are `Option<i32>` where
//! `None` replaces the C convention of a negative sentinel. All alignment
//! functions write their result into a caller-owned [`ExtzResult`], which
//! also owns the CIGAR buffer so that the allocation can be recycled across
//! calls (the equivalent of the `m_cigar`/`n_cigar`/`cigar` triple of C).
//!
//! ```
//! use rsksw2::*;
//!
//! let score_matrix = gen_simple_mat(5, 2, -4); // match +2, mismatch -4, wildcard row/col 0
//! let scoring = Scoring::new(5, &score_matrix, 4, 2);
//! let encode = |s: &[u8]| -> Vec<u8> {
//!     s.iter().map(|&c| b"ACGT"[..].iter().position(|&x| x == c).unwrap() as u8).collect()
//! };
//! let query = encode(b"ACGTACGT");
//! let target = encode(b"ACGTACGT");
//! let mut result = ExtzResult::new();
//! extz2_sse(&query, &target, &scoring, None, None, 0, Flags::empty(), &mut result);
//! assert_eq!(result.score, 16);
//! assert_eq!(result.cigar, vec![8 << 4]);
//! ```

#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::identity_op)]

mod common;
mod extd;
mod extz;
mod gg;
mod gg2;
#[cfg(target_arch = "x86_64")]
mod sse;

/// A score that is "minus infinity" in the C implementation (used as the
/// initial value of [`ExtzResult::score`] and friends).
pub const NEG_INF: i32 = -0x4000_0000;

// Raw flag bits, shared between the public `Flags` type and the private
// literal ports (which keep the C-shaped `flags: i32` parameter).
pub(crate) const EZ_SCORE_ONLY: i32 = 0x01; // don't record alignment path/cigar
pub(crate) const EZ_RIGHT: i32 = 0x02; // right-align gaps
pub(crate) const EZ_GENERIC_SC: i32 = 0x04; // without this flag: match/mismatch only; last symbol is a wildcard
pub(crate) const EZ_APPROX_MAX: i32 = 0x08; // approximate max; this is faster with sse
pub(crate) const EZ_APPROX_DROP: i32 = 0x10; // approximate Z-drop; faster with sse
pub(crate) const EZ_EXTZ_ONLY: i32 = 0x40; // only perform extension
pub(crate) const EZ_REV_CIGAR: i32 = 0x80; // reverse CIGAR in the output
pub(crate) const EZ_SPLICE_FOR: i32 = 0x100;
pub(crate) const EZ_SPLICE_REV: i32 = 0x200;
pub(crate) const EZ_SPLICE_FLANK: i32 = 0x400;
pub(crate) const EZ_EQX: i32 = 0x800;

/// The CIGAR operators emitted by the ksw algorithms, with their BAM codes.
///
/// A CIGAR is a `Vec<u32>` of words encoded as `length << 4 | code`; see
/// [`ExtzResult::cigar`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum CigarOp {
    /// `M`: aligned pair, match or mismatch (unless `Flags::EQX` is set).
    Match = 0,
    /// `I`: insertion to the target.
    Ins = 1,
    /// `D`: deletion from the target.
    Del = 2,
    /// `N`: skipped region of the target (a long intron in splice alignment).
    NSkip = 3,
    /// `=`: match (only with `Flags::EQX`).
    Eq = 7,
    /// `X`: mismatch (only with `Flags::EQX`).
    X = 8,
}

impl CigarOp {
    /// The BAM operator code (the low 4 bits of an encoded CIGAR word).
    #[inline]
    pub fn code(self) -> u32 {
        self as u32
    }

    /// Decode the operator stored in the low 4 bits of an encoded CIGAR word.
    /// Returns `None` for codes the ksw code never emits (4, 5, 6 and >= 9).
    #[inline]
    pub fn from_code(code: u32) -> Option<Self> {
        match code {
            0 => Some(Self::Match),
            1 => Some(Self::Ins),
            2 => Some(Self::Del),
            3 => Some(Self::NSkip),
            7 => Some(Self::Eq),
            8 => Some(Self::X),
            _ => None,
        }
    }

    /// The operator's SAM letter: `M`, `I`, `D`, `N`, `=` or `X`.
    #[inline]
    pub fn letter(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Ins => 'I',
            Self::Del => 'D',
            Self::NSkip => 'N',
            Self::Eq => '=',
            Self::X => 'X',
        }
    }
}

/// Alignment option flags — a type-safe replacement for the raw bit mask the
/// C API takes as its `flag` parameter.
///
/// Flags combine with `|`:
///
/// ```
/// use rsksw2::Flags;
/// let flags = Flags::EXTZ_ONLY | Flags::REV_CIGAR;
/// assert!(flags.contains(Flags::EXTZ_ONLY));
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub struct Flags(u32);

impl Flags {
    /// No flags set.
    pub const fn empty() -> Self {
        Self(0)
    }

    /// Don't record the alignment path (no CIGAR); score only.
    pub const SCORE_ONLY: Self = Self(EZ_SCORE_ONLY as u32);
    /// Right-align gaps.
    pub const RIGHT: Self = Self(EZ_RIGHT as u32);
    /// Look scores up in the full [`Scoring::score_matrix`] instead of using
    /// the fast match/mismatch path (SSE extension functions only; the last
    /// alphabet symbol is then no longer a wildcard).
    pub const GENERIC_SC: Self = Self(EZ_GENERIC_SC as u32);
    /// Track the maximum approximately via diagonal differences (faster with
    /// SSE; Z-drop then requires [`Flags::APPROX_DROP`]).
    pub const APPROX_MAX: Self = Self(EZ_APPROX_MAX as u32);
    /// Z-drop in approximate mode; use together with [`Flags::APPROX_MAX`].
    pub const APPROX_DROP: Self = Self(EZ_APPROX_DROP as u32);
    /// Only perform extension: backtrack to the best cell instead of to (0, 0).
    pub const EXTZ_ONLY: Self = Self(EZ_EXTZ_ONLY as u32);
    /// Reverse the CIGAR in the output (for workflows that reverse-complement
    /// the sequences before extending).
    pub const REV_CIGAR: Self = Self(EZ_REV_CIGAR as u32);
    /// Enable forward-strand splice signals (GT donor / AG acceptor).
    pub const SPLICE_FOR: Self = Self(EZ_SPLICE_FOR as u32);
    /// Enable reverse-strand splice signals (CT donor / AC acceptor).
    pub const SPLICE_REV: Self = Self(EZ_SPLICE_REV as u32);
    /// Penalise semi-canonical splice sites by half of the non-canonical
    /// penalty ([`exts2_sse`] only).
    pub const SPLICE_FLANK: Self = Self(EZ_SPLICE_FLANK as u32);
    /// Split `M` into `=`/`X` in the CIGAR ([`extd2_sse`] only).
    pub const EQX: Self = Self(EZ_EQX as u32);

    /// The underlying bit mask.
    #[inline]
    pub const fn bits(self) -> u32 {
        self.0
    }

    /// Build a flag set from a raw bit mask, as the C API would receive it.
    /// Unknown bits are preserved (the algorithms ignore them).
    #[inline]
    pub const fn from_bits(bits: u32) -> Self {
        Self(bits)
    }

    /// `true` if no flag is set.
    #[inline]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// `true` if all flags of `other` are set in `self`.
    #[inline]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl core::ops::BitOr for Flags {
    type Output = Self;
    #[inline]
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl core::ops::BitOrAssign for Flags {
    #[inline]
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl core::fmt::Debug for Flags {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        const NAMED: &[(u32, &str)] = &[
            (EZ_SCORE_ONLY as u32, "SCORE_ONLY"),
            (EZ_RIGHT as u32, "RIGHT"),
            (EZ_GENERIC_SC as u32, "GENERIC_SC"),
            (EZ_APPROX_MAX as u32, "APPROX_MAX"),
            (EZ_APPROX_DROP as u32, "APPROX_DROP"),
            (EZ_EXTZ_ONLY as u32, "EXTZ_ONLY"),
            (EZ_REV_CIGAR as u32, "REV_CIGAR"),
            (EZ_SPLICE_FOR as u32, "SPLICE_FOR"),
            (EZ_SPLICE_REV as u32, "SPLICE_REV"),
            (EZ_SPLICE_FLANK as u32, "SPLICE_FLANK"),
            (EZ_EQX as u32, "EQX"),
        ];
        f.write_str("Flags(")?;
        let mut first = true;
        let mut rest = self.0;
        for &(bit, name) in NAMED {
            if rest & bit != 0 {
                if !first {
                    f.write_str("|")?;
                }
                f.write_str(name)?;
                first = false;
                rest &= !bit;
            }
        }
        if rest != 0 {
            if !first {
                f.write_str("|")?;
            }
            write!(f, "{rest:#x}")?;
        }
        if first {
            f.write_str("empty")?;
        }
        f.write_str(")")
    }
}

/// The scoring system shared by all alignment entry points except
/// [`extf2_sse`], which takes a plain match/mismatch/gap-extension triple.
///
/// `score_matrix` is a row-major `alphabet_size × alphabet_size` matrix:
/// `score_matrix[i * alphabet_size + j]` is the score of target symbol `i`
/// against query symbol `j`. By convention the last row/column holds the
/// wildcard (score 0) — see [`gen_simple_mat`]. Gap costs are positive
/// numbers: a gap of length `l` costs `gap_open + l * gap_extension`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scoring {
    /// Number of distinct symbols `m`; each sequence byte must be `< m`.
    pub alphabet_size: usize,
    /// Row-major `m × m` scoring matrix.
    pub score_matrix: Vec<i8>,
    /// Gap open cost.
    pub gap_open: i8,
    /// Gap extension cost.
    pub gap_extension: i8,
}

impl Scoring {
    /// Build a scoring system from an existing matrix. The matrix is copied.
    ///
    /// # Panics
    /// If `alphabet_size` is not in `1..=127` or `score_matrix` does not hold
    /// exactly `alphabet_size * alphabet_size` entries.
    pub fn new(alphabet_size: usize, score_matrix: &[i8], gap_open: i8, gap_extension: i8) -> Self {
        assert!(alphabet_size >= 1, "alphabet_size must be at least 1");
        assert!(
            alphabet_size <= i8::MAX as usize,
            "alphabet_size must fit in i8 (the SSE ports pass it as int8_t)"
        );
        assert!(
            score_matrix.len() == alphabet_size * alphabet_size,
            "score_matrix must hold alphabet_size^2 entries, got {} for alphabet_size {alphabet_size}",
            score_matrix.len()
        );
        Self {
            alphabet_size,
            score_matrix: score_matrix.to_vec(),
            gap_open,
            gap_extension,
        }
    }
}

/// Alignment result — port of `ksw_extz_t`.
///
/// `max` is stored as a 31-bit unsigned integer, mirroring the
/// `uint32_t max:31` bit-field of the C struct (the C code relies on the
/// truncation when it assigns a negative score).
///
/// All extension entry points reset the struct on entry (the equivalent of
/// `ksw_reset_extz()`), so one instance can be reused across calls while the
/// CIGAR buffer's capacity is recycled.
#[derive(Clone, Debug, Default)]
pub struct ExtzResult {
    /// Maximum cell score seen along the whole alignment.
    pub max: u32,
    /// Whether the extension stopped early due to Z-drop / X-drop.
    pub zdropped: bool,
    /// Query coordinate of `max`.
    pub max_q: i32,
    /// Target coordinate of `max`.
    pub max_t: i32,
    /// Maximum score upon reaching the end of the query.
    pub mqe: i32,
    /// Target coordinate of `mqe`.
    pub mqe_t: i32,
    /// Maximum score upon reaching the end of the target.
    pub mte: i32,
    /// Query coordinate of `mte`.
    pub mte_q: i32,
    /// Score of the cell that reaches both sequence ends, or [`NEG_INF`]
    /// when no such cell was computed.
    pub score: i32,
    /// CIGAR as BAM-encoded words (`length << 4 | CIGAR operator code`, see
    /// [`CigarOp`]).
    pub cigar: Vec<u32>,
    /// Whether the extension reached the end of the query and `mqe +
    /// end_bonus > max` (only set by [`extz2_sse`]/[`extd2_sse`] with
    /// `Flags::EXTZ_ONLY`).
    pub reach_end: bool,
}

impl ExtzResult {
    /// A zeroed result (the equivalent of `memset(&ez, 0, sizeof(ez))`
    /// followed by `ksw_reset_extz()`).
    pub fn new() -> Self {
        let mut result = Self::default();
        result.reset();
        result
    }

    /// Port of `ksw_reset_extz()`; keeps the CIGAR allocation for reuse.
    pub fn reset(&mut self) {
        common::reset_extension_result(self);
    }

    /// Number of CIGAR operations (C's `n_cigar`).
    #[inline]
    pub fn n_cigar(&self) -> usize {
        self.cigar.len()
    }

    /// Render the CIGAR as a SAM string, e.g. `10M2I5M`.
    pub fn cigar_string(&self) -> String {
        let mut rendered = String::new();
        for &encoded in &self.cigar {
            let operation = CigarOp::from_code(encoded & 0xf);
            rendered.push_str(&(encoded >> 4).to_string());
            rendered.push(operation.map_or('?', CigarOp::letter));
        }
        rendered
    }
}

/* ------------------------------------------------------------------ */
/* Global alignment                                                    */
/* ------------------------------------------------------------------ */

/// Global alignment with an affine gap cost, Green's formulation — port of
/// `ksw2_gg.c` (`ksw_gg`).
///
/// Returns the alignment score. When `cigar` is `Some`, the CIGAR is written
/// into the vector (its contents are replaced but its capacity is reused);
/// when it is `None` only the score is computed.
///
/// `band_width` limits the band half-width; `None` means unlimited. A full
/// end-to-end alignment needs `w >= |query.len() - target.len()|`.
pub fn gg(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    band_width: Option<i32>,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg::ksw_gg(
        query.len() as i32,
        query,
        target.len() as i32,
        target,
        scoring.alphabet_size as i8,
        &scoring.score_matrix,
        scoring.gap_open,
        scoring.gap_extension,
        band_width.unwrap_or(-1),
        cigar,
    )
}

/// Global alignment with an affine gap cost, Suzuki's diagonal formulation —
/// port of `ksw2_gg2.c` (`ksw_gg2`).
///
/// See [`gg`] for the shared parameters.
pub fn gg2(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    band_width: Option<i32>,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg2::ksw_gg2(
        query.len() as i32,
        query,
        target.len() as i32,
        target,
        scoring.alphabet_size as i8,
        &scoring.score_matrix,
        scoring.gap_open,
        scoring.gap_extension,
        band_width.unwrap_or(-1),
        cigar,
    )
}

/// Global alignment with an affine gap cost and SSE4.1 intrinsics — port of
/// `ksw2_gg2_sse.c` (`ksw_gg2_sse`); always produces a CIGAR.
///
/// # Panics
/// On a CPU without the `sse4.1` feature.
#[cfg(target_arch = "x86_64")]
pub fn gg2_sse(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    band_width: Option<i32>,
    cigar: &mut Vec<u32>,
) -> i32 {
    require_sse41_feature();
    unsafe {
        sse::gg2_sse::ksw_gg2_sse(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            scoring.alphabet_size as i8,
            &scoring.score_matrix,
            scoring.gap_open,
            scoring.gap_extension,
            band_width.unwrap_or(-1),
            cigar,
        )
    }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, affine gap cost                                */
/* ------------------------------------------------------------------ */

/// NW-like extension with an affine gap cost, Green's formulation — port of
/// `ksw2_extz.c` (`ksw_extz`). Results are written into `result`.
///
/// `z_drop` disables the Z-drop heuristic when `None`; with `Some(z)`, the
/// extension stops early (and sets [`ExtzResult::zdropped`]) once the score
/// falls more than `z` plus gap-extension corrections below the running
/// maximum.
pub fn extz(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    band_width: Option<i32>,
    z_drop: Option<i32>,
    flags: Flags,
    result: &mut ExtzResult,
) {
    extz::ksw_extz(
        query.len() as i32,
        query,
        target.len() as i32,
        target,
        scoring.alphabet_size as i8,
        &scoring.score_matrix,
        scoring.gap_open,
        scoring.gap_extension,
        band_width.unwrap_or(-1),
        z_drop.unwrap_or(-1),
        flags.bits() as i32,
        result,
    );
}

/// NW-like extension with an affine gap cost and SSE4.1 intrinsics — port of
/// `ksw2_extz2_sse.c` (`ksw_extz2_sse`).
///
/// `end_bonus` is added to [`ExtzResult::mqe`] when deciding whether the
/// extension reached the end of the query (`Flags::EXTZ_ONLY` only).
///
/// # Panics
/// On a CPU without the `sse4.1` feature.
#[cfg(target_arch = "x86_64")]
pub fn extz2_sse(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    band_width: Option<i32>,
    z_drop: Option<i32>,
    end_bonus: i32,
    flags: Flags,
    result: &mut ExtzResult,
) {
    require_sse41_feature();
    unsafe {
        sse::extz2_sse::ksw_extz2_sse(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            scoring.alphabet_size as i8,
            &scoring.score_matrix,
            scoring.gap_open,
            scoring.gap_extension,
            band_width.unwrap_or(-1),
            z_drop.unwrap_or(-1),
            end_bonus,
            flags.bits() as i32,
            result,
        )
    }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, two-piece affine gap cost                      */
/* ------------------------------------------------------------------ */

/// NW-like extension with a two-piece affine gap cost, Green's formulation —
/// port of `ksw2_extd.c` (`ksw_extd`).
///
/// A gap of length `l` costs
/// `min(gap_open + l*gap_extension, gap_open2 + l*gap_extension2)`.
/// `z_drop` uses `gap_extension2` as the slope.
pub fn extd(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    gap_open2: i8,
    gap_extension2: i8,
    band_width: Option<i32>,
    z_drop: Option<i32>,
    flags: Flags,
    result: &mut ExtzResult,
) {
    extd::ksw_extd(
        query.len() as i32,
        query,
        target.len() as i32,
        target,
        scoring.alphabet_size as i8,
        &scoring.score_matrix,
        scoring.gap_open,
        scoring.gap_extension,
        gap_open2,
        gap_extension2,
        band_width.unwrap_or(-1),
        z_drop.unwrap_or(-1),
        flags.bits() as i32,
        result,
    );
}

/// NW-like extension with a two-piece affine gap cost and SSE4.1 intrinsics —
/// port of `ksw2_extd2_sse.c` (`ksw_extd2_sse`); the only entry point that
/// supports `Flags::EQX`.
///
/// # Panics
/// On a CPU without the `sse4.1` feature.
#[cfg(target_arch = "x86_64")]
pub fn extd2_sse(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    gap_open2: i8,
    gap_extension2: i8,
    band_width: Option<i32>,
    z_drop: Option<i32>,
    end_bonus: i32,
    flags: Flags,
    result: &mut ExtzResult,
) {
    require_sse41_feature();
    unsafe {
        sse::extd2_sse::ksw_extd2_sse(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            scoring.alphabet_size as i8,
            &scoring.score_matrix,
            scoring.gap_open,
            scoring.gap_extension,
            gap_open2,
            gap_extension2,
            band_width.unwrap_or(-1),
            z_drop.unwrap_or(-1),
            end_bonus,
            flags.bits() as i32,
            result,
        )
    }
}

/// Extension with a linear gap cost and X-drop, score only — port of
/// `ksw2_extf2_sse.c` (`ksw_extf2_sse`).
///
/// Takes `match_score`/`mismatch_score`/`gap_extension` instead of a
/// [`Scoring`]: symbols compare by byte equality (no wildcard), a gap of
/// length `l` costs `l * gap_extension`, and `x_drop` stops the extension
/// once the current diagonal drops more than the threshold below the running
/// maximum (setting [`ExtzResult::zdropped`]). No CIGAR is produced.
///
/// # Panics
/// On a CPU without the `sse4.1` feature.
#[cfg(target_arch = "x86_64")]
pub fn extf2_sse(
    query: &[u8],
    target: &[u8],
    match_score: i8,
    mismatch_score: i8,
    gap_extension: i8,
    band_width: Option<i32>,
    x_drop: Option<i32>,
    result: &mut ExtzResult,
) {
    require_sse41_feature();
    unsafe {
        sse::extf2_sse::ksw_extf2_sse(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            match_score,
            mismatch_score,
            gap_extension,
            band_width.unwrap_or(-1),
            x_drop.unwrap_or(-1),
            result,
        )
    }
}

/// Splice-aware extension with SSE4.1 intrinsics — port of
/// `ksw2_exts2_sse.c` (`ksw_exts2_sse`).
///
/// `gap_open2` is the flat cost of a long intron (the total penalty stops
/// growing once it reaches `gap_open2`), `non_canonical_penalty` prices
/// non-canonical splice sites, and `junction_bonus` is added at hinted splice
/// sites listed in `junctions` (one byte per target base; may be empty).
///
/// # Panics
/// On a CPU without the `sse4.1` feature.
#[cfg(target_arch = "x86_64")]
pub fn exts2_sse(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    gap_open2: i8,
    non_canonical_penalty: i8,
    z_drop: Option<i32>,
    junction_bonus: i8,
    flags: Flags,
    junctions: &[u8],
    result: &mut ExtzResult,
) {
    require_sse41_feature();
    unsafe {
        sse::exts2_sse::ksw_exts2_sse(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            scoring.alphabet_size as i8,
            &scoring.score_matrix,
            scoring.gap_open,
            scoring.gap_extension,
            gap_open2,
            non_canonical_penalty,
            z_drop.unwrap_or(-1),
            junction_bonus,
            flags.bits() as i32,
            junctions,
            result,
        )
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
fn require_sse41_feature() {
    assert!(
        std::arch::is_x86_feature_detected!("sse4.1"),
        "the *_sse entry points require the sse4.1 CPU feature; \
         use the scalar extz()/extd()/gg()/gg2() variants instead"
    );
}

/// Build a simple `alphabet_size × alphabet_size` scoring matrix — port of
/// `ksw_gen_simple_mat()` from `cli.c`: `match_score` on the diagonal,
/// `mismatch_score` off the diagonal, and 0 for the wildcard row/column
/// `alphabet_size - 1`. The match score is normalised to be positive and the
/// mismatch score to be negative, like the C helper.
///
/// # Panics
/// If `alphabet_size` is 0.
pub fn gen_simple_mat(alphabet_size: usize, match_score: i8, mismatch_score: i8) -> Vec<i8> {
    assert!(alphabet_size >= 1, "alphabet_size must be at least 1");
    let match_score = if match_score < 0 {
        -match_score
    } else {
        match_score
    };
    let mismatch_score = if mismatch_score > 0 {
        -mismatch_score
    } else {
        mismatch_score
    };
    let m = alphabet_size;
    let mut matrix = vec![0i8; m * m];
    for symbol in 0..m - 1 {
        for other_symbol in 0..m - 1 {
            matrix[symbol * m + other_symbol] = if symbol == other_symbol {
                match_score
            } else {
                mismatch_score
            };
        }
    }
    matrix
}
