//! Rust port of [KSW2](https://github.com/lh3/ksw2), a library for pairwise
//! biological sequence alignment based on dynamic programming.
//!
//! The port keeps a 1:1 correspondence with the C sources that live in the
//! `ksw2/` directory of this repository:
//!
//! | C source                | Rust module                | public function     |
//! |-------------------------|----------------------------|---------------------|
//! | `ksw2_gg.c`             | [`gg`]                     | [`ksw_gg`]          |
//! | `ksw2_gg2.c`            | [`gg2`]                    | [`ksw_gg2`]         |
//! | `ksw2_gg2_sse.c`        | `sse::gg2_sse`             | [`ksw_gg2_sse`]     |
//! | `ksw2_extz.c`           | [`extz`]                   | [`ksw_extz`]        |
//! | `ksw2_extz2_sse.c`      | `sse::extz2_sse`           | [`ksw_extz2_sse`]   |
//! | `ksw2_extd.c`           | [`extd`]                   | [`ksw_extd`]        |
//! | `ksw2_extd2_sse.c`      | `sse::extd2_sse`           | [`ksw_extd2_sse`]   |
//! | `ksw2_extf2_sse.c`      | `sse::extf2_sse`           | [`ksw_extf2_sse`]   |
//! | `ksw2_exts2_sse.c`      | `sse::exts2_sse`           | [`ksw_exts2_sse`]   |
//!
//! The `*_sse` entry points are literal ports of the SSE4.1 code paths of the
//! C sources (the ones selected by `-march=native` on any x86 CPU released
//! after 2008); they produce bit-identical scores, coordinates and CIGARs.
//! They require the `sse4.1` target feature at run time and are only compiled
//! on `x86_64`, exactly like the `#ifdef __SSE2__` guards of the C sources.
//! On other targets use the scalar entry points.
//!
//! All functions write their result into a caller-owned [`KswExtz`], which
//! also owns the CIGAR buffer so that the allocation can be recycled across
//! calls (the equivalent of the `m_cigar`/`n_cigar`/`cigar` triple of C).
//!
//! ```
//! use rsksw2::*;
//!
//! // simple 5x5 scoring matrix: match +2, mismatch -4, wildcard row/col 0
//! let mat: [i8; 25] = [2, -4, -4, -4, 0,
//!                     -4, 2, -4, -4, 0,
//!                     -4, -4, 2, -4, 0,
//!                     -4, -4, -4, 2, 0,
//!                      0,  0,  0,  0, 0];
//! let query: Vec<u8> = b"ACGTACGT".iter().map(|&c| b"ACGT"[..].iter().position(|&x| x == c).unwrap() as u8).collect();
//! let target: Vec<u8> = b"ACGTACGT".iter().map(|&c| b"ACGT"[..].iter().position(|&x| x == c).unwrap() as u8).collect();
//! let mut ez = KswExtz::new();
//! ksw_extz2_sse(8, &query, 8, &target, 5, &mat, 4, 2, -1, -1, 0, 0, &mut ez);
//! assert_eq!(ez.score, 16);
//! assert_eq!(ez.cigar, vec![8 << 4]);
//! ```

#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::identity_op)]

mod common;
pub mod extd;
pub mod extz;
pub mod gg;
pub mod gg2;
#[cfg(target_arch = "x86_64")]
mod sse;

/// A score that is "minus infinity" in the C implementation.
pub const KSW_NEG_INF: i32 = -0x4000_0000;

pub const KSW_EZ_SCORE_ONLY: i32 = 0x01; // don't record alignment path/cigar
pub const KSW_EZ_RIGHT: i32 = 0x02; // right-align gaps
pub const KSW_EZ_GENERIC_SC: i32 = 0x04; // without this flag: match/mismatch only; last symbol is a wildcard
pub const KSW_EZ_APPROX_MAX: i32 = 0x08; // approximate max; this is faster with sse
pub const KSW_EZ_APPROX_DROP: i32 = 0x10; // approximate Z-drop; faster with sse
pub const KSW_EZ_EXTZ_ONLY: i32 = 0x40; // only perform extension
pub const KSW_EZ_REV_CIGAR: i32 = 0x80; // reverse CIGAR in the output
pub const KSW_EZ_SPLICE_FOR: i32 = 0x100;
pub const KSW_EZ_SPLICE_REV: i32 = 0x200;
pub const KSW_EZ_SPLICE_FLANK: i32 = 0x400;
pub const KSW_EZ_EQX: i32 = 0x800;

/// The subset of CIGAR operators used by ksw code.
pub const KSW_CIGAR_MATCH: u32 = 0;
pub const KSW_CIGAR_INS: u32 = 1;
pub const KSW_CIGAR_DEL: u32 = 2;
pub const KSW_CIGAR_N_SKIP: u32 = 3;
pub const KSW_CIGAR_EQ: u32 = 7;
pub const KSW_CIGAR_X: u32 = 8;

/// Port of `ksw_extz_t`.
///
/// `max` is stored as a 31-bit unsigned integer, mirroring the
/// `uint32_t max:31` bit-field of the C struct (the C code relies on the
/// truncation when it assigns a negative score).
#[derive(Clone, Debug, Default)]
pub struct KswExtz {
    pub max: u32,
    pub zdropped: bool,
    pub max_q: i32,
    pub max_t: i32,
    pub mqe: i32,
    pub mqe_t: i32,
    pub mte: i32,
    pub mte_q: i32,
    pub score: i32,
    pub cigar: Vec<u32>,
    pub reach_end: i32,
}

impl KswExtz {
    /// Equivalent to `memset(&ez, 0, sizeof(ez))` followed by `ksw_reset_extz()`.
    pub fn new() -> Self {
        let mut ez = Self::default();
        ez.reset();
        ez
    }

    /// Port of `ksw_reset_extz()`; keeps the CIGAR allocation for reuse.
    pub fn reset(&mut self) {
        common::reset_extz(self);
    }

    /// Number of CIGAR operations (C's `n_cigar`).
    #[inline]
    pub fn n_cigar(&self) -> usize {
        self.cigar.len()
    }

    /// Render the CIGAR as a SAM string, e.g. `10M2I5M`.
    pub fn cigar_string(&self) -> String {
        const OP: [char; 9] = ['M', 'I', 'D', 'N', '?', '?', '?', '=', 'X'];
        let mut s = String::new();
        for &c in &self.cigar {
            let op = (c & 0xf) as usize;
            s.push_str(&(c >> 4).to_string());
            s.push(if op < OP.len() { OP[op] } else { '?' });
        }
        s
    }
}

/* ------------------------------------------------------------------ */
/* Global alignment                                                    */
/* ------------------------------------------------------------------ */

/// Global alignment, Green's formulation (`ksw2_gg.c`).
///
/// Returns the alignment score. When `cigar` is `Some`, the CIGAR is written
/// into it (the buffer is cleared but its capacity is reused); when it is
/// `None` only the score is computed.
pub fn ksw_gg(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    w: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg::ksw_gg(qlen, query, tlen, target, m, mat, gapo, gape, w, cigar)
}

/// Global alignment, Suzuki's diagonal formulation (`ksw2_gg2.c`).
pub fn ksw_gg2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    w: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg2::ksw_gg2(qlen, query, tlen, target, m, mat, gapo, gape, w, cigar)
}

/// Global alignment with SSE intrinsics (`ksw2_gg2_sse.c`); always produces a CIGAR.
#[cfg(target_arch = "x86_64")]
pub fn ksw_gg2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    w: i32,
    cigar: &mut Vec<u32>,
) -> i32 {
    require_sse41();
    unsafe { sse::gg2_sse::ksw_gg2_sse(qlen, query, tlen, target, m, mat, gapo, gape, w, cigar) }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, affine gap cost                                */
/* ------------------------------------------------------------------ */

/// NW-like extension, Green's formulation (`ksw2_extz.c`).
pub fn ksw_extz(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    w: i32,
    zdrop: i32,
    flag: i32,
    ez: &mut KswExtz,
) {
    extz::ksw_extz(qlen, query, tlen, target, m, mat, q, e, w, zdrop, flag, ez)
}

/// NW-like extension with SSE intrinsics (`ksw2_extz2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extz2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    w: i32,
    zdrop: i32,
    end_bonus: i32,
    flag: i32,
    ez: &mut KswExtz,
) {
    require_sse41();
    unsafe {
        sse::extz2_sse::ksw_extz2_sse(
            qlen, query, tlen, target, m, mat, q, e, w, zdrop, end_bonus, flag, ez,
        )
    }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, two-piece affine gap cost                      */
/* ------------------------------------------------------------------ */

/// NW-like extension with a two-piece affine gap cost, Green's formulation (`ksw2_extd.c`).
pub fn ksw_extd(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    gapo2: i8,
    gape2: i8,
    w: i32,
    zdrop: i32,
    flag: i32,
    ez: &mut KswExtz,
) {
    extd::ksw_extd(
        qlen, query, tlen, target, m, mat, gapo, gape, gapo2, gape2, w, zdrop, flag, ez,
    )
}

/// NW-like extension with a two-piece affine gap cost and SSE intrinsics
/// (`ksw2_extd2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extd2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
    w: i32,
    zdrop: i32,
    end_bonus: i32,
    flag: i32,
    ez: &mut KswExtz,
) {
    require_sse41();
    unsafe {
        sse::extd2_sse::ksw_extd2_sse(
            qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez,
        )
    }
}

/// Extension with a linear gap cost and X-drop, score only (`ksw2_extf2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extf2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    mch: i8,
    mis: i8,
    e: i8,
    w: i32,
    xdrop: i32,
    ez: &mut KswExtz,
) {
    require_sse41();
    unsafe { sse::extf2_sse::ksw_extf2_sse(qlen, query, tlen, target, mch, mis, e, w, xdrop, ez) }
}

/// Splice-aware extension with SSE intrinsics (`ksw2_exts2_sse.c`).
///
/// `junc` is an optional per-target-base junction hint array (may be empty).
#[cfg(target_arch = "x86_64")]
pub fn ksw_exts2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    q2: i8,
    noncan: i8,
    zdrop: i32,
    junc_bonus: i8,
    flag: i32,
    junc: &[u8],
    ez: &mut KswExtz,
) {
    require_sse41();
    unsafe {
        sse::exts2_sse::ksw_exts2_sse(
            qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, junc_bonus, flag, junc, ez,
        )
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
fn require_sse41() {
    assert!(
        std::arch::is_x86_feature_detected!("sse4.1"),
        "the *_sse entry points require the sse4.1 target feature; \
         use the scalar ksw_extz()/ksw_extd()/ksw_gg()/ksw_gg2() variants instead"
    );
}

/// Build a simple `m`x`m` scoring matrix, mirroring `ksw_gen_simple_mat()` of `cli.c`.
pub fn gen_simple_mat(m: i32, mat: &mut [i8], a: i8, b: i8) {
    let a = if a < 0 { -a } else { a };
    let b = if b > 0 { -b } else { b };
    let mu = m as usize;
    for i in 0..mu - 1 {
        for j in 0..mu - 1 {
            mat[i * mu + j] = if i == j { a } else { b };
        }
        mat[i * mu + mu - 1] = 0;
    }
    for j in 0..mu {
        mat[(mu - 1) * mu + j] = 0;
    }
}
