//! Helpers shared by all implementations; direct ports of the `static inline`
//! functions at the bottom of `ksw2.h`.

use crate::{
    KSW_CIGAR_DEL, KSW_CIGAR_EQ, KSW_CIGAR_INS, KSW_CIGAR_MATCH, KSW_CIGAR_N_SKIP, KSW_CIGAR_X,
    KSW_NEG_INF, KswExtz,
};

/// Port of `ksw_push_cigar()`. The C version manages `m_cigar`/`n_cigar`
/// manually; `Vec` keeps the capacity across calls, which is equivalent.
#[inline]
pub(crate) fn push_cigar_operation(cigar: &mut Vec<u32>, operation: u32, length: i32) {
    match cigar.last_mut() {
        Some(last) if (*last & 0xf) == operation => {
            *last = last.wrapping_add((length as u32) << 4)
        }
        _ => cigar.push(((length as u32) << 4) | operation),
    }
}

/// Port of `ksw_reset_extz()`. The CIGAR becomes empty but its buffer is kept.
#[inline]
pub(crate) fn reset_extension_result(result: &mut KswExtz) {
    result.max_q = -1;
    result.max_t = -1;
    result.mqe_t = -1;
    result.mte_q = -1;
    result.max = 0;
    result.score = KSW_NEG_INF;
    result.mqe = KSW_NEG_INF;
    result.mte = KSW_NEG_INF;
    result.cigar.clear();
    result.zdropped = false;
    result.reach_end = 0;
}

/// Port of `ksw_apply_zdrop()`. Returns `true` when the extension should stop.
///
/// The two coordinate arguments are the ones the caller has at hand: for the
/// rotated (anti-diagonal) formulation they are `(anti_diagonal, target_index)`,
/// for the row-wise formulation `(target_index, query_index)`; `is_rotated`
/// selects which convention is in use.
#[inline]
pub(crate) fn apply_z_drop(
    result: &mut KswExtz,
    is_rotated: bool,
    score: i32,
    first_index: i32,
    second_index: i32,
    z_drop: i32,
    gap_extension: i8,
) -> bool {
    let (anti_diagonal, target_index) = if is_rotated {
        (first_index, second_index)
    } else {
        (first_index + second_index, first_index)
    };
    if score > result.max as i32 {
        result.max = (score as u32) & 0x7fff_ffff; // `uint32_t max:31` truncation
        result.max_t = target_index;
        result.max_q = anti_diagonal - target_index;
    } else if target_index >= result.max_t && anti_diagonal - target_index >= result.max_q {
        let target_delta = target_index - result.max_t;
        let query_delta = (anti_diagonal - target_index) - result.max_q;
        let gap_length = if target_delta > query_delta {
            target_delta - query_delta
        } else {
            query_delta - target_delta
        };
        if z_drop >= 0 && result.max as i32 - score > z_drop + gap_length * gap_extension as i32 {
            result.zdropped = true;
            return true;
        }
    }
    false
}

/// Port of `ksw_backtrack()`.
///
/// `backtrack_matrix` is the backtrack matrix, `band_start`/`band_end` the per-row
/// (or per-anti-diagonal when `is_rotated` is set) band boundaries, `row_width`
/// the number of bytes per row.
pub(crate) fn backtrack(
    is_rotated: bool,
    is_reversed: bool,
    min_intron_length: i32,
    backtrack_matrix: &[u8],
    band_start: &[i32],
    band_end: Option<&[i32]>,
    row_width: i32,
    start_target_index: i32,
    start_query_index: i32,
    cigar: &mut Vec<u32>,
) {
    cigar.clear();
    let row_width = row_width as usize;
    let mut target_index = start_target_index;
    let mut query_index = start_query_index;
    let mut state: i32 = 0;
    while target_index >= 0 && query_index >= 0 {
        // at the beginning of the loop, _state_ tells us which state to check
        let mut forced_state = -1;
        let cell: u8;
        if is_rotated {
            let anti_diagonal = target_index + query_index;
            let anti_diagonal_usize = anti_diagonal as usize;
            if target_index < band_start[anti_diagonal_usize] {
                forced_state = 2;
            }
            if let Some(band_end_row) = band_end {
                if target_index > band_end_row[anti_diagonal_usize] {
                    forced_state = 1;
                }
            }
            cell = if forced_state < 0 {
                let row_offset = anti_diagonal_usize * row_width;
                let column_offset = (target_index - band_start[anti_diagonal_usize]) as usize;
                backtrack_matrix[row_offset + column_offset]
            } else {
                0
            };
        } else {
            let target_index_usize = target_index as usize;
            if query_index < band_start[target_index_usize] {
                forced_state = 2;
            }
            if let Some(band_end_row) = band_end {
                if query_index > band_end_row[target_index_usize] {
                    forced_state = 1;
                }
            }
            cell = if forced_state < 0 {
                let row_offset = target_index_usize * row_width;
                let column_offset = (query_index - band_start[target_index_usize]) as usize;
                backtrack_matrix[row_offset + column_offset]
            } else {
                0
            };
        }
        if state == 0 {
            state = (cell & 7) as i32; // if requesting the H state, find state one maximizes it.
        } else if (cell >> (state + 2)) & 1 == 0 {
            state = 0; // if requesting other states, _state_ stays the same if it is a continuation
        }
        if state == 0 {
            state = (cell & 7) as i32;
        }
        if forced_state >= 0 {
            state = forced_state;
        }
        if state == 0 {
            push_cigar_operation(cigar, KSW_CIGAR_MATCH, 1);
            target_index -= 1;
            query_index -= 1;
        } else if state == 1 || (state == 3 && min_intron_length <= 0) {
            push_cigar_operation(cigar, KSW_CIGAR_DEL, 1);
            target_index -= 1;
        } else if state == 3 {
            push_cigar_operation(cigar, KSW_CIGAR_N_SKIP, 1);
            target_index -= 1;
        } else {
            push_cigar_operation(cigar, KSW_CIGAR_INS, 1);
            query_index -= 1;
        }
    }
    if target_index >= 0 {
        // first deletion
        let operation = if min_intron_length > 0 && target_index >= min_intron_length {
            KSW_CIGAR_N_SKIP
        } else {
            KSW_CIGAR_DEL
        };
        push_cigar_operation(cigar, operation, target_index + 1);
    }
    if query_index >= 0 {
        // first insertion
        push_cigar_operation(cigar, KSW_CIGAR_INS, query_index + 1);
    }
    if !is_reversed {
        cigar.reverse(); // reverse CIGAR
    }
}

/// Port of `ksw_cigar2eqx()`: expand `M` into `=`/`X`. `input_cigar` must
/// already be a copy of the CIGAR to expand; `output_cigar` receives the result.
pub(crate) fn expand_cigar_matches(
    query: &[u8],
    target: &[u8],
    input_cigar: &[u32],
    output_cigar: &mut Vec<u32>,
) {
    output_cigar.clear();
    let mut target_position = 0usize;
    let mut query_position = 0usize;
    for &encoded in input_cigar {
        let operation = encoded & 0xf;
        let length = (encoded >> 4) as usize;
        if operation == KSW_CIGAR_MATCH {
            for offset in 0..length {
                if target[target_position + offset] == query[query_position + offset] {
                    push_cigar_operation(output_cigar, KSW_CIGAR_EQ, 1);
                } else {
                    push_cigar_operation(output_cigar, KSW_CIGAR_X, 1);
                }
            }
            target_position += length;
            query_position += length;
        } else {
            push_cigar_operation(output_cigar, operation, length as i32);
            if operation == KSW_CIGAR_DEL || operation == KSW_CIGAR_N_SKIP {
                target_position += length;
            } else if operation == KSW_CIGAR_INS {
                query_position += length;
            } else if operation == KSW_CIGAR_EQ || operation == KSW_CIGAR_X {
                target_position += length;
                query_position += length;
            }
        }
    }
}
