//! Port of `ksw2_gg2_sse.c` (SSE4.1 path).

use core::arch::x86_64::*;

use crate::common::backtrack;

// Backtrack byte layout (see `ksw2.h`): the low 3 bits hold the state that
// maximises H (0 = match, 1 = deletion/E, 2 = insertion/F, 3 = intron,
// 4 = second gap set), while bit `state + 2` records whether that state
// continues. The vectors below encode those values.
const DELETION_STATE: i8 = 1;
const INSERTION_STATE: i8 = 2;
const DELETION_CONTINUATION: i8 = 0x08;
const INSERTION_CONTINUATION: i8 = 0x10;

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_gg2_sse(
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
    unsafe {
        let alphabet_size_usize = alphabet_size as usize;
        let gap_open_i32 = gap_open as i32;
        let gap_open_plus_extension = gap_open_i32 + gap_extension as i32;
        let query_len_usize = query_len as usize;
        let target_len_usize = target_len as usize;

        let zero_vector = _mm_setzero_si128();
        let gap_open_vector = _mm_set1_epi8(gap_open);
        let twice_gap_open_plus_extension_vector =
            _mm_set1_epi8((gap_open_plus_extension * 2) as i8);
        let deletion_state_vector = _mm_set1_epi8(DELETION_STATE);
        let insertion_state_vector = _mm_set1_epi8(INSERTION_STATE);
        let deletion_continuation_vector = _mm_set1_epi8(DELETION_CONTINUATION);
        let insertion_continuation_vector = _mm_set1_epi8(INSERTION_CONTINUATION);

        let band_width = if band_width < 0 {
            if target_len > query_len {
                target_len
            } else {
                query_len
            }
        } else {
            band_width
        };
        // number of columns in the backtrack matrix
        let backtrack_columns = if band_width + 1 < target_len {
            band_width + 1
        } else {
            target_len
        };
        let target_block_count = ((target_len + 15) / 16) as usize;
        let backtrack_blocks = ((backtrack_columns + 15) / 16 + 1) as usize;
        let backtrack_row_bytes = (backtrack_blocks * 16) as i32;

        let mut u = vec![zero_vector; target_block_count];
        let mut v = vec![zero_vector; target_block_count];
        let mut x = vec![zero_vector; target_block_count];
        let mut y = vec![zero_vector; target_block_count];
        let mut s = vec![zero_vector; target_block_count];
        let mut padded_target = vec![0u8; target_block_count * 16];
        let mut reverse_query = vec![0u8; query_len_usize];
        let row_count = (query_len + target_len - 1) as usize;
        let mut backtrack_matrix = vec![0u8; (row_count * backtrack_blocks + 1) * 16];
        let mut band_start_per_row = vec![0i32; row_count];

        for query_index in 0..query_len_usize {
            reverse_query[query_index] = query[query_len_usize - 1 - query_index];
        }
        padded_target[..target_len_usize].copy_from_slice(&target[..target_len_usize]);

        let u_ptr = u.as_mut_ptr() as *mut u8;
        let v_ptr = v.as_mut_ptr() as *mut u8;
        let s_ptr = s.as_mut_ptr() as *mut u8;
        let padded_target_ptr = padded_target.as_ptr();
        let reverse_query_ptr = reverse_query.as_ptr();
        let score_matrix_ptr = score_matrix.as_ptr();

        let mut diagonal_score = 0i32;
        let mut diagonal_score_target_index = 0i32;
        let mut previous_band_start = -1i32;
        let mut previous_band_end = -1i32;

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
            let unpadded_band_start = band_start;
            let unpadded_band_end = band_end;
            band_start = band_start / 16 * 16;
            band_end = (band_end + 16) / 16 * 16 - 1;
            band_start_per_row[anti_diagonal as usize] = band_start;
            // set boundary conditions
            let previous_x: i8;
            let previous_v: i8;
            if band_start > 0 {
                if band_start - 1 >= previous_band_start && band_start - 1 <= previous_band_end {
                    // (anti_diagonal-1, band_start-1) calculated in the last round
                    previous_x = *(x.as_ptr() as *const i8).add((band_start - 1) as usize);
                    previous_v = *v_ptr.add((band_start - 1) as usize) as i8;
                } else {
                    // not calculated; set to zeros
                    previous_x = 0;
                    previous_v = 0;
                }
            } else {
                previous_x = 0;
                previous_v = if anti_diagonal != 0 { gap_open } else { 0 };
            }
            if band_end >= anti_diagonal {
                *(y.as_ptr() as *mut u8).add(anti_diagonal as usize) = 0;
                *u_ptr.add(anti_diagonal as usize) = if anti_diagonal != 0 {
                    gap_open as u8
                } else {
                    0
                };
            }
            // loop fission: set scores first
            let reverse_query_base = query_len_usize as isize - 1 - anti_diagonal as isize;
            let mut target_index = unpadded_band_start;
            while target_index <= unpadded_band_end {
                let target_index_usize = target_index as usize;
                *s_ptr.add(target_index_usize) = *score_matrix_ptr.add(
                    *padded_target_ptr.add(target_index_usize) as usize * alphabet_size_usize
                        + *reverse_query_ptr.offset(reverse_query_base + target_index as isize)
                            as usize,
                ) as u8;
                target_index += 1;
            }
            // core loop
            let mut previous_x_vector = _mm_cvtsi32_si128(previous_x as i32);
            let mut previous_v_vector = _mm_cvtsi32_si128(previous_v as i32);
            let block_start = (band_start / 16) as usize;
            let block_end = (band_end / 16) as usize;
            let backtrack_row_offset =
                anti_diagonal as usize * backtrack_blocks * 16 - band_start as usize;
            let mut block = block_start;
            while block <= block_end {
                let mut max_score_vector =
                    _mm_add_epi8(_mm_load_si128(s.as_ptr().add(block)), twice_gap_open_plus_extension_vector);

                // x_previous_block <- x[r-1][t-1..t+14], the shifted-out lane is carried
                // over in previous_x_vector
                let mut x_previous_block = _mm_load_si128(x.as_ptr().add(block));
                let mut temporary = _mm_srli_si128::<15>(x_previous_block); // x[r-1][t+15]
                x_previous_block = _mm_or_si128(_mm_slli_si128::<1>(x_previous_block), previous_x_vector);
                previous_x_vector = temporary;
                let mut v_previous_block = _mm_load_si128(v.as_ptr().add(block));
                temporary = _mm_srli_si128::<15>(v_previous_block); // v[r-1][t+15]
                v_previous_block = _mm_or_si128(_mm_slli_si128::<1>(v_previous_block), previous_v_vector);
                previous_v_vector = temporary;
                // candidate_from_x <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]
                let mut candidate_from_x = _mm_add_epi8(x_previous_block, v_previous_block);

                let u_current_block = _mm_load_si128(u.as_ptr().add(block));
                // candidate_from_y <- y[r-1][t..t+15] + u[r-1][t..t+15]
                let mut candidate_from_y =
                    _mm_add_epi8(_mm_load_si128(y.as_ptr().add(block)), u_current_block);

                let mut backtrack_byte_vector = _mm_and_si128(
                    _mm_cmpgt_epi8(candidate_from_x, max_score_vector),
                    deletion_state_vector,
                );
                max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x);
                temporary = _mm_cmpgt_epi8(candidate_from_y, max_score_vector);
                backtrack_byte_vector =
                    _mm_blendv_epi8(backtrack_byte_vector, insertion_state_vector, temporary);
                // max(z, candidate_from_y); this works because both are non-negative
                max_score_vector = _mm_max_epu8(max_score_vector, candidate_from_y);
                // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
                _mm_store_si128(
                    u.as_mut_ptr().add(block),
                    _mm_sub_epi8(max_score_vector, v_previous_block),
                );
                // v[r][t..t+15] <- z - u[r-1][t..t+15]
                _mm_store_si128(
                    v.as_mut_ptr().add(block),
                    _mm_sub_epi8(max_score_vector, u_current_block),
                );

                max_score_vector = _mm_sub_epi8(max_score_vector, gap_open_vector);
                candidate_from_x = _mm_sub_epi8(candidate_from_x, max_score_vector);
                candidate_from_y = _mm_sub_epi8(candidate_from_y, max_score_vector);
                temporary = _mm_cmpgt_epi8(candidate_from_x, zero_vector);
                backtrack_byte_vector =
                    _mm_or_si128(backtrack_byte_vector, _mm_and_si128(deletion_continuation_vector, temporary));
                _mm_store_si128(
                    x.as_mut_ptr().add(block),
                    _mm_and_si128(candidate_from_x, temporary),
                );
                temporary = _mm_cmpgt_epi8(candidate_from_y, zero_vector);
                backtrack_byte_vector =
                    _mm_or_si128(backtrack_byte_vector, _mm_and_si128(insertion_continuation_vector, temporary));
                _mm_store_si128(
                    y.as_mut_ptr().add(block),
                    _mm_and_si128(candidate_from_y, temporary),
                );
                _mm_storeu_si128(
                    backtrack_matrix
                        .as_mut_ptr()
                        .add(backtrack_row_offset + block * 16) as *mut __m128i,
                    backtrack_byte_vector,
                );
                block += 1;
            }
            if anti_diagonal > 0 {
                if diagonal_score_target_index >= unpadded_band_start
                    && diagonal_score_target_index <= unpadded_band_end
                {
                    diagonal_score +=
                        *v_ptr.add(diagonal_score_target_index as usize) as i32 - gap_open_plus_extension;
                } else {
                    diagonal_score_target_index += 1;
                    diagonal_score +=
                        *u_ptr.add(diagonal_score_target_index as usize) as i32 - gap_open_plus_extension;
                }
            } else {
                diagonal_score = *v_ptr as i32 - 2 * gap_open_plus_extension;
                diagonal_score_target_index = 0;
            }
            previous_band_start = band_start;
            previous_band_end = band_end;
        }

        backtrack(
            true,
            false,
            0,
            &backtrack_matrix[..row_count * backtrack_blocks * 16],
            &band_start_per_row,
            None,
            backtrack_row_bytes,
            target_len - 1,
            query_len - 1,
            cigar,
        );
        diagonal_score
    }
}
