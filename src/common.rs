//! Helpers shared by all implementations; direct ports of the `static inline`
//! functions at the bottom of `ksw2.h`.

use crate::{
    KSW_CIGAR_DEL, KSW_CIGAR_EQ, KSW_CIGAR_INS, KSW_CIGAR_MATCH, KSW_CIGAR_N_SKIP, KSW_CIGAR_X,
    KSW_NEG_INF, KswExtz,
};

/// Port of `ksw_push_cigar()`. The C version manages `m_cigar`/`n_cigar`
/// manually; `Vec` keeps the capacity across calls, which is equivalent.
#[inline]
pub(crate) fn push_cigar(cigar: &mut Vec<u32>, op: u32, len: i32) {
    match cigar.last_mut() {
        Some(last) if (*last & 0xf) == op => *last = last.wrapping_add((len as u32) << 4),
        _ => cigar.push(((len as u32) << 4) | op),
    }
}

/// Port of `ksw_reset_extz()`. `n_cigar` becomes 0 but the buffer is kept.
#[inline]
pub(crate) fn reset_extz(ez: &mut KswExtz) {
    ez.max_q = -1;
    ez.max_t = -1;
    ez.mqe_t = -1;
    ez.mte_q = -1;
    ez.max = 0;
    ez.score = KSW_NEG_INF;
    ez.mqe = KSW_NEG_INF;
    ez.mte = KSW_NEG_INF;
    ez.cigar.clear();
    ez.zdropped = false;
    ez.reach_end = 0;
}

/// Port of `ksw_apply_zdrop()`. Returns `true` when the extension should stop.
#[inline]
pub(crate) fn apply_zdrop(
    ez: &mut KswExtz,
    is_rot: bool,
    h: i32,
    a: i32,
    b: i32,
    zdrop: i32,
    e: i8,
) -> bool {
    let (r, t) = if is_rot { (a, b) } else { (a + b, a) };
    if h > ez.max as i32 {
        ez.max = (h as u32) & 0x7fff_ffff; // `uint32_t max:31` truncation
        ez.max_t = t;
        ez.max_q = r - t;
    } else if t >= ez.max_t && r - t >= ez.max_q {
        let tl = t - ez.max_t;
        let ql = (r - t) - ez.max_q;
        let l = if tl > ql { tl - ql } else { ql - tl };
        if zdrop >= 0 && ez.max as i32 - h > zdrop + l * e as i32 {
            ez.zdropped = true;
            return true;
        }
    }
    false
}

/// Port of `ksw_backtrack()`.
///
/// `p` is the backtrack matrix, `off`/`off_end` the per-row (or per-anti-diagonal
/// when `is_rot` is set) band boundaries, `n_col` the number of bytes per row.
pub(crate) fn backtrack(
    is_rot: bool,
    is_rev: bool,
    min_intron_len: i32,
    p: &[u8],
    off: &[i32],
    off_end: Option<&[i32]>,
    n_col: i32,
    i0: i32,
    j0: i32,
    cigar: &mut Vec<u32>,
) {
    cigar.clear();
    let n_col = n_col as usize;
    let mut i = i0;
    let mut j = j0;
    let mut state: i32 = 0;
    while i >= 0 && j >= 0 {
        // at the beginning of the loop, _state_ tells us which state to check
        let mut force_state = -1;
        let tmp: u8;
        if is_rot {
            let r = i + j;
            let ru = r as usize;
            if i < off[ru] {
                force_state = 2;
            }
            if let Some(oe) = off_end {
                if i > oe[ru] {
                    force_state = 1;
                }
            }
            tmp = if force_state < 0 {
                p[ru * n_col + (i - off[ru]) as usize]
            } else {
                0
            };
        } else {
            let iu = i as usize;
            if j < off[iu] {
                force_state = 2;
            }
            if let Some(oe) = off_end {
                if j > oe[iu] {
                    force_state = 1;
                }
            }
            tmp = if force_state < 0 {
                p[iu * n_col + (j - off[iu]) as usize]
            } else {
                0
            };
        }
        if state == 0 {
            state = (tmp & 7) as i32; // if requesting the H state, find state one maximizes it.
        } else if (tmp >> (state + 2)) & 1 == 0 {
            state = 0; // if requesting other states, _state_ stays the same if it is a continuation
        }
        if state == 0 {
            state = (tmp & 7) as i32;
        }
        if force_state >= 0 {
            state = force_state;
        }
        if state == 0 {
            push_cigar(cigar, KSW_CIGAR_MATCH, 1);
            i -= 1;
            j -= 1;
        } else if state == 1 || (state == 3 && min_intron_len <= 0) {
            push_cigar(cigar, KSW_CIGAR_DEL, 1);
            i -= 1;
        } else if state == 3 {
            push_cigar(cigar, KSW_CIGAR_N_SKIP, 1);
            i -= 1;
        } else {
            push_cigar(cigar, KSW_CIGAR_INS, 1);
            j -= 1;
        }
    }
    if i >= 0 {
        // first deletion
        let op = if min_intron_len > 0 && i >= min_intron_len {
            KSW_CIGAR_N_SKIP
        } else {
            KSW_CIGAR_DEL
        };
        push_cigar(cigar, op, i + 1);
    }
    if j >= 0 {
        // first insertion
        push_cigar(cigar, KSW_CIGAR_INS, j + 1);
    }
    if !is_rev {
        cigar.reverse(); // reverse CIGAR
    }
}

/// Port of `ksw_cigar2eqx()`: expand `M` into `=`/`X`. `ci0` must already be a
/// copy of the CIGAR to expand; `ci1` receives the result.
pub(crate) fn cigar2eqx(query: &[u8], target: &[u8], ci0: &[u32], ci1: &mut Vec<u32>) {
    ci1.clear();
    let mut x = 0usize;
    let mut y = 0usize;
    for &c in ci0 {
        let op = c & 0xf;
        let len = (c >> 4) as usize;
        if op == KSW_CIGAR_MATCH {
            for i in 0..len {
                if target[x + i] == query[y + i] {
                    push_cigar(ci1, KSW_CIGAR_EQ, 1);
                } else {
                    push_cigar(ci1, KSW_CIGAR_X, 1);
                }
            }
            x += len;
            y += len;
        } else {
            push_cigar(ci1, op, len as i32);
            if op == KSW_CIGAR_DEL || op == KSW_CIGAR_N_SKIP {
                x += len;
            } else if op == KSW_CIGAR_INS {
                y += len;
            } else if op == KSW_CIGAR_EQ || op == KSW_CIGAR_X {
                x += len;
                y += len;
            }
        }
    }
}
