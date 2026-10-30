//! Port of `ksw2_extz.c`: global/extension alignment with an affine gap cost,
//! Green's standard formulation.

use crate::common::{apply_z_drop, backtrack, reset_extension_result};
use crate::{
    KSW_EZ_EXTZ_ONLY, KSW_EZ_REV_CIGAR, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF, KswExtz,
};

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
    reset_extension_result(result);
    if alphabet_size <= 0 || query_len <= 0 || target_len <= 0 {
        return;
    }

    let alphabet_size_usize = alphabet_size as usize;
    let gap_extension_i32 = gap_extension as i32;
    let gap_open_plus_extension = gap_open as i32 + gap_extension_i32;
    let band_width = if band_width < 0 {
        if target_len > query_len {
            target_len
        } else {
            query_len
        }
    } else {
        band_width
    };
    // maximum number of columns of the backtrack matrix
    let backtrack_columns = if query_len < 2 * band_width + 1 {
        query_len
    } else {
        2 * band_width + 1
    };
    let record_cigar = (flags & KSW_EZ_SCORE_ONLY) == 0;
    let query_len_usize = query_len as usize;
    let target_len_usize = target_len as usize;
    let backtrack_columns_usize = backtrack_columns as usize;

    // allocate memory
    let mut query_profile = vec![0i8; query_len_usize * alphabet_size_usize];
    let mut column_scores = vec![[0i32; 2]; query_len_usize + 1]; // { h, e }
    let mut backtrack_matrix: Vec<u8> = Vec::new();
    let mut band_start_per_row: Vec<i32> = Vec::new();
    if record_cigar {
        backtrack_matrix = vec![0u8; backtrack_columns_usize * target_len_usize];
        band_start_per_row = vec![0i32; target_len_usize];
    }

    // generate the query profile
    let mut profile_index = 0usize;
    for symbol in 0..alphabet_size_usize {
        let matrix_row_start = symbol * alphabet_size_usize;
        let matrix_row = &score_matrix[matrix_row_start..matrix_row_start + alphabet_size_usize];
        for query_index in 0..query_len_usize {
            query_profile[profile_index] = matrix_row[query[query_index] as usize];
            profile_index += 1;
        }
    }

    // fill the first row
    column_scores[0] = [0, -gap_open_plus_extension - gap_open_plus_extension];
    let mut query_index = 1i32;
    while query_index <= query_len && query_index <= band_width {
        column_scores[query_index as usize] = [
            -(gap_open_plus_extension + gap_extension_i32 * (query_index - 1)),
            -(gap_open_plus_extension + gap_open_plus_extension + gap_extension_i32 * query_index),
        ];
        query_index += 1;
    }
    while query_index <= query_len {
        // everything is -inf outside the band
        column_scores[query_index as usize] = [KSW_NEG_INF, KSW_NEG_INF];
        query_index += 1;
    }

    let query_profile_ptr = query_profile.as_ptr();
    let column_scores_ptr = column_scores.as_mut_ptr();
    let backtrack_matrix_ptr = backtrack_matrix.as_mut_ptr();
    let band_start_ptr = band_start_per_row.as_mut_ptr();
    let target_ptr = target.as_ptr();
    let right_align = (flags & KSW_EZ_RIGHT) != 0;
    let mut max_score_query_index = 0i32;

    // DP loop (target sequence is in the outer loop)
    for target_index in 0..target_len {
        let mut max_score = KSW_NEG_INF;
        let band_start = if target_index > band_width {
            target_index - band_width
        } else {
            0
        };
        let band_end = if target_index + band_width < query_len - 1 {
            target_index + band_width
        } else {
            query_len - 1
        };
        // a gap of length l costs -(gap_open + l*gap_extension)
        let leading_gap_score = -(gap_open_plus_extension + gap_extension_i32 * target_index);
        let mut h_from_left = if band_start > 0 {
            KSW_NEG_INF
        } else {
            leading_gap_score
        };
        let mut f_gap_score = if band_start > 0 {
            KSW_NEG_INF
        } else {
            leading_gap_score - gap_open_plus_extension
        };
        unsafe {
            let profile_row =
                query_profile_ptr.add(*target_ptr.add(target_index as usize) as usize * query_len_usize);
            let mut cell_ptr = column_scores_ptr.add(band_start as usize);
            if !record_cigar {
                for query_index in band_start..=band_end {
                    let cell = &mut *cell_ptr;
                    let h_diagonal = cell[0];
                    let mut e_gap_score = cell[1];
                    cell[0] = h_from_left;
                    let mut h = h_diagonal + *profile_row.add(query_index as usize) as i32;
                    h = h.max(e_gap_score);
                    h = h.max(f_gap_score);
                    h_from_left = h;
                    let is_greater = max_score > h;
                    max_score_query_index = if is_greater { max_score_query_index } else { query_index };
                    max_score = if is_greater { max_score } else { h };
                    h -= gap_open_plus_extension;
                    e_gap_score -= gap_extension_i32;
                    e_gap_score = e_gap_score.max(h);
                    cell[1] = e_gap_score;
                    f_gap_score -= gap_extension_i32;
                    f_gap_score = f_gap_score.max(h);
                    cell_ptr = cell_ptr.add(1);
                }
            } else if !right_align {
                let mut backtrack_ptr =
                    backtrack_matrix_ptr.add(target_index as usize * backtrack_columns_usize);
                *band_start_ptr.add(target_index as usize) = band_start;
                for query_index in band_start..=band_end {
                    let cell = &mut *cell_ptr;
                    let h_diagonal = cell[0];
                    let mut e_gap_score = cell[1];
                    cell[0] = h_from_left;
                    let mut h = h_diagonal + *profile_row.add(query_index as usize) as i32;
                    let mut backtrack_byte: u8 = if h >= e_gap_score { 0 } else { 1 };
                    h = h.max(e_gap_score);
                    backtrack_byte = if h >= f_gap_score { backtrack_byte } else { 2 };
                    h = h.max(f_gap_score);
                    h_from_left = h;
                    let is_greater = max_score > h;
                    max_score_query_index = if is_greater { max_score_query_index } else { query_index };
                    max_score = if is_greater { max_score } else { h };
                    h -= gap_open_plus_extension;
                    e_gap_score -= gap_extension_i32;
                    backtrack_byte |= if e_gap_score > h { 0x08 } else { 0 };
                    e_gap_score = e_gap_score.max(h);
                    cell[1] = e_gap_score;
                    f_gap_score -= gap_extension_i32;
                    backtrack_byte |= if f_gap_score > h { 0x10 } else { 0 };
                    f_gap_score = f_gap_score.max(h);
                    // the backtrack byte keeps h for the current cell and e/f for the next cell
                    *backtrack_ptr = backtrack_byte;
                    backtrack_ptr = backtrack_ptr.add(1);
                    cell_ptr = cell_ptr.add(1);
                }
            } else {
                let mut backtrack_ptr =
                    backtrack_matrix_ptr.add(target_index as usize * backtrack_columns_usize);
                *band_start_ptr.add(target_index as usize) = band_start;
                for query_index in band_start..=band_end {
                    let cell = &mut *cell_ptr;
                    let h_diagonal = cell[0];
                    let mut e_gap_score = cell[1];
                    cell[0] = h_from_left;
                    let mut h = h_diagonal + *profile_row.add(query_index as usize) as i32;
                    let mut backtrack_byte: u8 = if h > e_gap_score { 0 } else { 1 };
                    h = h.max(e_gap_score);
                    backtrack_byte = if h > f_gap_score { backtrack_byte } else { 2 };
                    h = h.max(f_gap_score);
                    h_from_left = h;
                    let is_greater_or_equal = max_score >= h;
                    max_score_query_index = if is_greater_or_equal {
                        max_score_query_index
                    } else {
                        query_index
                    };
                    max_score = if is_greater_or_equal { max_score } else { h };
                    h -= gap_open_plus_extension;
                    e_gap_score -= gap_extension_i32;
                    backtrack_byte |= if e_gap_score >= h { 0x08 } else { 0 };
                    e_gap_score = e_gap_score.max(h);
                    cell[1] = e_gap_score;
                    f_gap_score -= gap_extension_i32;
                    backtrack_byte |= if f_gap_score >= h { 0x10 } else { 0 };
                    f_gap_score = f_gap_score.max(h);
                    *backtrack_ptr = backtrack_byte;
                    backtrack_ptr = backtrack_ptr.add(1);
                    cell_ptr = cell_ptr.add(1);
                }
            }
            let next_column = ((band_end + 1) as usize).min(query_len_usize);
            (*column_scores_ptr.add(next_column))[0] = h_from_left;
            (*column_scores_ptr.add(next_column))[1] = KSW_NEG_INF;
        }
        // update the result
        if band_end == query_len - 1 && column_scores[query_len_usize][0] > result.mqe {
            result.mqe = column_scores[query_len_usize][0];
            result.mqe_t = target_index;
        }
        if target_index == target_len - 1 {
            result.mte = max_score;
            result.mte_q = max_score_query_index;
        }
        if apply_z_drop(
            result,
            false,
            max_score,
            target_index,
            max_score_query_index,
            z_drop,
            gap_extension,
        ) {
            break;
        }
        if target_index == target_len - 1 && band_end == query_len - 1 {
            result.score = column_scores[query_len_usize][0];
        }
    }

    if record_cigar {
        let reverse_cigar = (flags & KSW_EZ_REV_CIGAR) != 0;
        if !result.zdropped && (flags & KSW_EZ_EXTZ_ONLY) == 0 {
            backtrack(
                false,
                reverse_cigar,
                0,
                &backtrack_matrix,
                &band_start_per_row,
                None,
                backtrack_columns,
                target_len - 1,
                query_len - 1,
                &mut result.cigar,
            );
        } else if result.max_t >= 0 && result.max_q >= 0 {
            backtrack(
                false,
                reverse_cigar,
                0,
                &backtrack_matrix,
                &band_start_per_row,
                None,
                backtrack_columns,
                result.max_t,
                result.max_q,
                &mut result.cigar,
            );
        }
    }
}
