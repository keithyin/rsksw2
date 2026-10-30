//! Port of `ksw2_extf2_sse.c` (SSE4.1 path): score-only extension with a
//! linear gap cost and X-drop.

use core::arch::x86_64::*;

use crate::common::reset_extension_result;
use crate::KswExtz;

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_extf2_sse(
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
    unsafe {
        reset_extension_result(result);
        if query_len <= 0 || target_len <= 0 {
            return;
        }
        let query_len_usize = query_len as usize;
        let target_len_usize = target_len as usize;
        let gap_extension_i32 = gap_extension as i32;

        let twice_gap_extension_vector = _mm_set1_epi8((gap_extension_i32 * 2) as i8);
        let match_score_vector = _mm_set1_epi8(match_score);
        let mismatch_score_vector =
            _mm_set1_epi8(if mismatch_score < 0 { mismatch_score } else { -mismatch_score });
        let target_block_count = ((target_len + 15) / 16) as usize;
        let query_block_count = ((query_len + 15) / 16) as usize;
        let band_width = if band_width < 0 {
            if target_len > query_len {
                target_len
            } else {
                query_len
            }
        } else {
            band_width
        };

        let zero_vector = _mm_setzero_si128();
        let mut u = vec![zero_vector; target_block_count];
        let mut v = vec![zero_vector; target_block_count];
        let mut s = vec![zero_vector; target_block_count + 1];
        let mut padded_target = vec![0u8; target_block_count * 16 + 16];
        let mut reverse_query = vec![0u8; (query_block_count + 1) * 16];

        for query_index in 0..query_len_usize {
            reverse_query[query_index] = query[query_len_usize - 1 - query_index];
        }
        padded_target[..target_len_usize].copy_from_slice(&target[..target_len_usize]);

        let u_ptr = u.as_mut_ptr() as *mut u8;
        let v_ptr = v.as_mut_ptr() as *mut u8;

        let mut diagonal_score = 0i32;
        let mut diagonal_score_target_index = 0i32;
        let mut previous_band_start = -1i32;
        let mut previous_band_end = -1i32;
        let mut anti_diagonal = 0i32;
        while anti_diagonal < query_len + target_len - 1 {
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
            if band_start > band_end {
                break;
            }
            let unpadded_band_start = band_start;
            let unpadded_band_end = band_end;
            band_start = band_start / 16 * 16;
            band_end = (band_end + 16) / 16 * 16 - 1;
            // set boundary conditions
            let previous_v: u8 = if band_start > 0
                && band_start - 1 >= previous_band_start
                && band_start - 1 <= previous_band_end
            {
                *v_ptr.add((band_start - 1) as usize)
            } else {
                0
            };
            if band_end >= anti_diagonal {
                *u_ptr.add(anti_diagonal as usize) = 0;
            }
            // core loop
            let mut previous_v_vector = _mm_cvtsi32_si128(previous_v as i32);
            let block_start = (band_start / 16) as usize;
            let block_end = (band_end / 16) as usize;
            let reverse_query_base = query_len_usize as isize - 1 - anti_diagonal as isize;
            let mut target_index = unpadded_band_start;
            while target_index <= unpadded_band_end {
                let target_block =
                    _mm_loadu_si128(padded_target.as_ptr().add(target_index as usize) as *const __m128i);
                let query_block = _mm_loadu_si128(
                    reverse_query.as_ptr().offset(reverse_query_base + target_index as isize)
                        as *const __m128i,
                );
                let mut score_block = _mm_cmpeq_epi8(target_block, query_block);
                score_block = _mm_blendv_epi8(mismatch_score_vector, match_score_vector, score_block);
                _mm_storeu_si128(
                    (s.as_mut_ptr() as *mut u8).add(target_index as usize) as *mut __m128i,
                    score_block,
                );
                target_index += 16;
            }
            let mut block = block_start;
            while block <= block_end {
                let mut max_score_vector =
                    _mm_add_epi8(_mm_load_si128(s.as_ptr().add(block)), twice_gap_extension_vector);
                // v_previous_block <- v[r-1][t-1..t+14]
                let mut v_previous_block = _mm_load_si128(v.as_ptr().add(block));
                let temporary = _mm_srli_si128::<15>(v_previous_block); // v[r-1][t+15]
                v_previous_block =
                    _mm_or_si128(_mm_slli_si128::<1>(v_previous_block), previous_v_vector);
                previous_v_vector = temporary;
                let u_current_block = _mm_load_si128(u.as_ptr().add(block));
                max_score_vector = _mm_max_epi8(max_score_vector, v_previous_block);
                // this works because both are non-negative
                max_score_vector = _mm_max_epu8(max_score_vector, u_current_block);
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
                block += 1;
            }
            if anti_diagonal > 0 {
                if diagonal_score_target_index >= unpadded_band_start
                    && diagonal_score_target_index <= unpadded_band_end
                    && diagonal_score_target_index + 1 >= unpadded_band_start
                    && diagonal_score_target_index + 1 <= unpadded_band_end
                {
                    let score_from_v = *v_ptr.add(diagonal_score_target_index as usize) as i32 - gap_extension_i32;
                    let score_from_u =
                        *u_ptr.add(diagonal_score_target_index as usize + 1) as i32 - gap_extension_i32;
                    if score_from_v > score_from_u {
                        diagonal_score += score_from_v;
                    } else {
                        diagonal_score += score_from_u;
                        diagonal_score_target_index += 1;
                    }
                } else if diagonal_score_target_index >= unpadded_band_start
                    && diagonal_score_target_index <= unpadded_band_end
                {
                    diagonal_score +=
                        *v_ptr.add(diagonal_score_target_index as usize) as i32 - gap_extension_i32;
                } else {
                    diagonal_score_target_index += 1;
                    diagonal_score +=
                        *u_ptr.add(diagonal_score_target_index as usize) as i32 - gap_extension_i32;
                }
                if diagonal_score > result.max as i32 {
                    result.max = (diagonal_score as u32) & 0x7fff_ffff;
                    result.max_t = diagonal_score_target_index;
                    result.max_q = anti_diagonal - diagonal_score_target_index;
                } else if x_drop >= 0 && result.max as i32 - diagonal_score > x_drop {
                    break;
                }
            } else {
                diagonal_score = *v_ptr as i32 - gap_extension_i32 - gap_extension_i32;
                diagonal_score_target_index = 0;
            }
            previous_band_start = band_start;
            previous_band_end = band_end;
            anti_diagonal += 1;
        }
        if anti_diagonal == query_len + target_len - 1 {
            result.score = diagonal_score;
        } else {
            result.zdropped = true;
        }
    }
}
