//! Port of `ksw2_gg2.c`: global alignment, Suzuki's diagonal formulation.

use crate::common::backtrack;

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
    if alphabet_size <= 0 || query_len <= 0 || target_len <= 0 {
        return 0;
    }
    let alphabet_size_usize = alphabet_size as usize;
    let gap_open_plus_extension = gap_open as i32 + gap_extension as i32;
    let twice_gap_open_plus_extension = gap_open_plus_extension + gap_open_plus_extension;
    let band_width = if band_width < 0 {
        if target_len > query_len {
            target_len
        } else {
            query_len
        }
    } else {
        band_width
    };
    let backtrack_columns = if band_width + 1 < target_len {
        band_width + 1
    } else {
        target_len
    };
    let query_len_usize = query_len as usize;
    let target_len_usize = target_len as usize;
    let backtrack_columns_usize = backtrack_columns as usize;

    let mut u = vec![0i8; target_len_usize + 1];
    let mut v = vec![0i8; target_len_usize + 1];
    let mut x = vec![0i8; target_len_usize + 1];
    let mut y = vec![0i8; target_len_usize + 1];
    let mut s = vec![0i8; target_len_usize];
    let mut reverse_query = vec![0u8; query_len_usize];
    let mut backtrack_matrix: Vec<u8> = Vec::new();
    let mut band_start_per_row: Vec<i32> = Vec::new();
    let record_cigar = cigar.is_some();
    if record_cigar {
        backtrack_matrix = vec![0u8; (query_len_usize + target_len_usize) * backtrack_columns_usize];
        band_start_per_row = vec![0i32; query_len_usize + target_len_usize];
    }

    for query_index in 0..query_len_usize {
        reverse_query[query_index] = query[query_len_usize - 1 - query_index];
    }

    let mut diagonal_score = 0i32;
    let mut diagonal_score_target_index = 0i32;
    let u_ptr = u.as_mut_ptr();
    let v_ptr = v.as_mut_ptr();
    let x_ptr = x.as_mut_ptr();
    let y_ptr = y.as_mut_ptr();
    let s_ptr = s.as_mut_ptr();
    let reverse_query_ptr = reverse_query.as_ptr();
    let score_matrix_ptr = score_matrix.as_ptr();
    let target_ptr = target.as_ptr();
    let backtrack_matrix_ptr = backtrack_matrix.as_mut_ptr();
    let band_start_ptr = band_start_per_row.as_mut_ptr();

    for anti_diagonal in 0..query_len + target_len - 1 {
        let mut band_start = 0i32;
        let mut band_end = target_len - 1;
        // find the boundaries
        if band_start < anti_diagonal - query_len + 1 {
            band_start = anti_diagonal - query_len + 1;
        }
        if band_end > anti_diagonal {
            band_end = anti_diagonal;
        }
        if band_start < (anti_diagonal - band_width + 1) >> 1 {
            band_start = (anti_diagonal - band_width + 1) >> 1; // take the ceil
        }
        if band_end > (anti_diagonal + band_width) >> 1 {
            band_end = (anti_diagonal + band_width) >> 1; // take the floor
        }
        let band_start_usize = band_start as usize;
        let band_end_usize = band_end as usize;
        let reverse_query_offset = query_len_usize as isize - 1 - anti_diagonal as isize;
        unsafe {
            // set boundary conditions
            let mut previous_x: i8;
            let mut previous_v: i8;
            if band_start != 0 {
                if anti_diagonal > band_start + band_start + band_width - 1 {
                    previous_x = 0;
                    previous_v = 0;
                } else {
                    previous_x = *x_ptr.add(band_start_usize - 1);
                    previous_v = *v_ptr.add(band_start_usize - 1); // (anti_diagonal-1, band_start-1) in the band
                }
            } else {
                previous_x = 0;
                previous_v = if anti_diagonal != 0 { gap_open } else { 0 };
            }
            if band_end != anti_diagonal {
                if anti_diagonal < band_end + band_end - band_width - 1 {
                    *y_ptr.add(band_end_usize) = 0;
                    *u_ptr.add(band_end_usize) = 0; // (anti_diagonal-1,band_end) out of the band
                }
            } else {
                *y_ptr.add(anti_diagonal as usize) = 0;
                *u_ptr.add(anti_diagonal as usize) = if anti_diagonal != 0 { gap_open } else { 0 };
            }
            // loop fission: set scores first
            for target_index in band_start_usize..=band_end_usize {
                *s_ptr.add(target_index) = *score_matrix_ptr.add(
                    *target_ptr.add(target_index) as usize * alphabet_size_usize
                        + *reverse_query_ptr.offset(target_index as isize + reverse_query_offset) as usize,
                );
            }
            // core loop
            if record_cigar {
                let backtrack_row_offset = anti_diagonal as usize * backtrack_columns_usize;
                *band_start_ptr.add(anti_diagonal as usize) = band_start;
                for target_index in band_start_usize..=band_end_usize {
                    /* candidate_x = x(r-1,t-1) + v(r-1,t-1)
                       candidate_y = y(r-1,t)   + u(r-1,t)
                       max_score   = max{ S(t,r-t) + 2q + 2r, candidate_x, candidate_y } */
                    let mut max_score: i8 = (*s_ptr.add(target_index) as i32 + twice_gap_open_plus_extension) as i8;
                    let mut candidate_x: i8 = (previous_x as i32 + previous_v as i32) as i8;
                    let candidate_y: i8 = (*y_ptr.add(target_index) as i32 + *u_ptr.add(target_index) as i32) as i8;
                    let mut backtrack_byte: u8 = if candidate_x > max_score { 1 } else { 0 };
                    max_score = max_score.max(candidate_x);
                    backtrack_byte = if candidate_y > max_score { 2 } else { backtrack_byte };
                    max_score = max_score.max(candidate_y);
                    let u_before = *u_ptr.add(target_index);
                    *u_ptr.add(target_index) = (max_score as i32 - previous_v as i32) as i8;
                    previous_v = *v_ptr.add(target_index);
                    *v_ptr.add(target_index) = (max_score as i32 - u_before as i32) as i8;
                    let score_minus_gap_open = (max_score as i32 - gap_open as i32) as i8;
                    candidate_x = (candidate_x as i32 - score_minus_gap_open as i32) as i8;
                    let candidate_y = (candidate_y as i32 - score_minus_gap_open as i32) as i8;
                    previous_x = *x_ptr.add(target_index);
                    backtrack_byte |= if candidate_x > 0 { 0x08 } else { 0 };
                    *x_ptr.add(target_index) = if candidate_x > 0 { candidate_x } else { 0 };
                    backtrack_byte |= if candidate_y > 0 { 0x10 } else { 0 };
                    *y_ptr.add(target_index) = if candidate_y > 0 { candidate_y } else { 0 };
                    *backtrack_matrix_ptr
                        .add(backtrack_row_offset + target_index - band_start_usize) = backtrack_byte;
                }
            } else {
                for target_index in band_start_usize..=band_end_usize {
                    let mut max_score: i8 = (*s_ptr.add(target_index) as i32 + twice_gap_open_plus_extension) as i8;
                    let mut candidate_x: i8 = (previous_x as i32 + previous_v as i32) as i8;
                    let candidate_y: i8 = (*y_ptr.add(target_index) as i32 + *u_ptr.add(target_index) as i32) as i8;
                    max_score = max_score.max(candidate_x);
                    max_score = max_score.max(candidate_y);
                    let u_before = *u_ptr.add(target_index);
                    *u_ptr.add(target_index) = (max_score as i32 - previous_v as i32) as i8;
                    previous_v = *v_ptr.add(target_index);
                    *v_ptr.add(target_index) = (max_score as i32 - u_before as i32) as i8;
                    let score_minus_gap_open = (max_score as i32 - gap_open as i32) as i8;
                    candidate_x = (candidate_x as i32 - score_minus_gap_open as i32) as i8;
                    let candidate_y = (candidate_y as i32 - score_minus_gap_open as i32) as i8;
                    previous_x = *x_ptr.add(target_index);
                    *x_ptr.add(target_index) = if candidate_x > 0 { candidate_x } else { 0 };
                    *y_ptr.add(target_index) = if candidate_y > 0 { candidate_y } else { 0 };
                }
            }
        }
        if anti_diagonal > 0 {
            if diagonal_score_target_index >= band_start && diagonal_score_target_index <= band_end {
                diagonal_score += v[diagonal_score_target_index as usize] as i32 - gap_open_plus_extension;
            } else {
                diagonal_score_target_index += 1;
                diagonal_score += u[diagonal_score_target_index as usize] as i32 - gap_open_plus_extension;
            }
        } else {
            diagonal_score = v[0] as i32 - gap_open_plus_extension - gap_open_plus_extension;
            diagonal_score_target_index = 0;
        }
    }

    if let Some(cigar) = cigar {
        backtrack(
            true,
            false,
            0,
            &backtrack_matrix,
            &band_start_per_row,
            None,
            backtrack_columns,
            target_len - 1,
            query_len - 1,
            cigar,
        );
    }
    diagonal_score
}
