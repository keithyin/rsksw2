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
//! let score_matrix: [i8; 25] = [2, -4, -4, -4, 0,
//!                              -4, 2, -4, -4, 0,
//!                              -4, -4, 2, -4, 0,
//!                              -4, -4, -4, 2, 0,
//!                               0,  0,  0,  0, 0];
//! let query: Vec<u8> = b"ACGTACGT".iter().map(|&c| b"ACGT"[..].iter().position(|&x| x == c).unwrap() as u8).collect();
//! let target: Vec<u8> = b"ACGTACGT".iter().map(|&c| b"ACGT"[..].iter().position(|&x| x == c).unwrap() as u8).collect();
//! let mut result = KswExtz::new();
//! ksw_extz2_sse(8, &query, 8, &target, 5, &score_matrix, 4, 2, -1, -1, 0, 0, &mut result);
//! assert_eq!(result.score, 16);
//! assert_eq!(result.cigar, vec![8 << 4]);
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
        const OPERATOR_BYTE: [char; 9] = ['M', 'I', 'D', 'N', '?', '?', '?', '=', 'X'];
        let mut rendered = String::new();
        for &encoded in &self.cigar {
            let operation = (encoded & 0xf) as usize;
            rendered.push_str(&(encoded >> 4).to_string());
            rendered.push(if operation < OPERATOR_BYTE.len() {
                OPERATOR_BYTE[operation]
            } else {
                '?'
            });
        }
        rendered
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
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg::ksw_gg(
        query_len,
        query,
        target_len,
        target,
        alphabet_size,
        score_matrix,
        gap_open,
        gap_extension,
        band_width,
        cigar,
    )
}

/// Global alignment, Suzuki's diagonal formulation (`ksw2_gg2.c`).
pub fn ksw_gg2(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    gg2::ksw_gg2(
        query_len,
        query,
        target_len,
        target,
        alphabet_size,
        score_matrix,
        gap_open,
        gap_extension,
        band_width,
        cigar,
    )
}

/// Global alignment with SSE intrinsics (`ksw2_gg2_sse.c`); always produces a CIGAR.
#[cfg(target_arch = "x86_64")]
pub fn ksw_gg2_sse(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    cigar: &mut Vec<u32>,
) -> i32 {
    require_sse41_feature();
    unsafe {
        sse::gg2_sse::ksw_gg2_sse(
            query_len,
            query,
            target_len,
            target,
            alphabet_size,
            score_matrix,
            gap_open,
            gap_extension,
            band_width,
            cigar,
        )
    }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, affine gap cost                                */
/* ------------------------------------------------------------------ */

/// NW-like extension, Green's formulation (`ksw2_extz.c`).
pub fn ksw_extz(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    z_drop: i32,
    flags: i32,
    result: &mut KswExtz,
) {
    extz::ksw_extz(
        query_len,
        query,
        target_len,
        target,
        alphabet_size,
        score_matrix,
        gap_open,
        gap_extension,
        band_width,
        z_drop,
        flags,
        result,
    )
}

/// NW-like extension with SSE intrinsics (`ksw2_extz2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extz2_sse(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    z_drop: i32,
    end_bonus: i32,
    flags: i32,
    result: &mut KswExtz,
) {
    require_sse41_feature();
    unsafe {
        sse::extz2_sse::ksw_extz2_sse(
            query_len,
            query,
            target_len,
            target,
            alphabet_size,
            score_matrix,
            gap_open,
            gap_extension,
            band_width,
            z_drop,
            end_bonus,
            flags,
            result,
        )
    }
}

/* ------------------------------------------------------------------ */
/* Extension alignment, two-piece affine gap cost                      */
/* ------------------------------------------------------------------ */

/// NW-like extension with a two-piece affine gap cost, Green's formulation (`ksw2_extd.c`).
pub fn ksw_extd(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    gap_open2: i8,
    gap_extension2: i8,
    band_width: i32,
    z_drop: i32,
    flags: i32,
    result: &mut KswExtz,
) {
    extd::ksw_extd(
        query_len,
        query,
        target_len,
        target,
        alphabet_size,
        score_matrix,
        gap_open,
        gap_extension,
        gap_open2,
        gap_extension2,
        band_width,
        z_drop,
        flags,
        result,
    )
}

/// NW-like extension with a two-piece affine gap cost and SSE intrinsics
/// (`ksw2_extd2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extd2_sse(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    gap_open2: i8,
    gap_extension2: i8,
    band_width: i32,
    z_drop: i32,
    end_bonus: i32,
    flags: i32,
    result: &mut KswExtz,
) {
    require_sse41_feature();
    unsafe {
        sse::extd2_sse::ksw_extd2_sse(
            query_len,
            query,
            target_len,
            target,
            alphabet_size,
            score_matrix,
            gap_open,
            gap_extension,
            gap_open2,
            gap_extension2,
            band_width,
            z_drop,
            end_bonus,
            flags,
            result,
        )
    }
}

/// Extension with a linear gap cost and X-drop, score only (`ksw2_extf2_sse.c`).
#[cfg(target_arch = "x86_64")]
pub fn ksw_extf2_sse(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    match_score: i8,
    mismatch_score: i8,
    gap_extension: i8,
    band_width: i32,
    x_drop: i32,
    result: &mut KswExtz,
) {
    require_sse41_feature();
    unsafe {
        sse::extf2_sse::ksw_extf2_sse(
            query_len,
            query,
            target_len,
            target,
            match_score,
            mismatch_score,
            gap_extension,
            band_width,
            x_drop,
            result,
        )
    }
}

/// Splice-aware extension with SSE intrinsics (`ksw2_exts2_sse.c`).
///
/// `junctions` is an optional per-target-base junction hint array (may be empty).
#[cfg(target_arch = "x86_64")]
pub fn ksw_exts2_sse(
    query_len: i32,
    query: &[u8],
    target_len: i32,
    target: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    gap_open2: i8,
    non_canonical_penalty: i8,
    z_drop: i32,
    junction_bonus: i8,
    flags: i32,
    junctions: &[u8],
    result: &mut KswExtz,
) {
    require_sse41_feature();
    unsafe {
        sse::exts2_sse::ksw_exts2_sse(
            query_len,
            query,
            target_len,
            target,
            alphabet_size,
            score_matrix,
            gap_open,
            gap_extension,
            gap_open2,
            non_canonical_penalty,
            z_drop,
            junction_bonus,
            flags,
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
        "the *_sse entry points require the sse4.1 target feature; \
         use the scalar ksw_extz()/ksw_extd()/ksw_gg()/ksw_gg2() variants instead"
    );
}

/// Build a simple `alphabet_size`x`alphabet_size` scoring matrix, mirroring
/// `ksw_gen_simple_mat()` of `cli.c`.
pub fn gen_simple_mat(
    alphabet_size: i32,
    score_matrix: &mut [i8],
    match_score: i8,
    mismatch_score: i8,
) {
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
    let alphabet_size_usize = alphabet_size as usize;
    for symbol in 0..alphabet_size_usize - 1 {
        for other_symbol in 0..alphabet_size_usize - 1 {
            score_matrix[symbol * alphabet_size_usize + other_symbol] = if symbol == other_symbol {
                match_score
            } else {
                mismatch_score
            };
        }
        score_matrix[symbol * alphabet_size_usize + alphabet_size_usize - 1] = 0;
    }
    for other_symbol in 0..alphabet_size_usize {
        score_matrix[(alphabet_size_usize - 1) * alphabet_size_usize + other_symbol] = 0;
    }
}
