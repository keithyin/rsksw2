//! Port of `ksw2_exts2_sse.c` (SSE4.1 path): splice-aware extension.

use core::arch::x86_64::*;

use crate::common::{apply_z_drop, backtrack, reset_extension_result};
use crate::{
    KSW_EZ_APPROX_DROP, KSW_EZ_APPROX_MAX, KSW_EZ_EXTZ_ONLY, KSW_EZ_GENERIC_SC,
    KSW_EZ_REV_CIGAR, KSW_EZ_RIGHT,
    KSW_EZ_SCORE_ONLY, KSW_EZ_SPLICE_FLANK, KSW_EZ_SPLICE_FOR, KSW_EZ_SPLICE_REV, KSW_NEG_INF,
    KswExtz,
};

// Backtrack byte layout (see `ksw2.h`): the low 3 bits hold the state that
// maximises H (0 = match, 1 = deletion, 2 = insertion, 3 = intron), while bit
// `state + 2` records whether that state continues.
const DELETION_STATE: i8 = 1;
const INSERTION_STATE: i8 = 2;
const INTRON_STATE: i8 = 3;
const DELETION_CONTINUATION: i8 = 0x08;
const INSERTION_CONTINUATION: i8 = 0x10;
const INTRON_CONTINUATION: i8 = 0x20;

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_exts2_sse(
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
    unsafe {
        reset_extension_result(result);
        if alphabet_size <= 1
            || query_len <= 0
            || target_len <= 0
            || gap_open2 as i32 <= gap_open as i32 + gap_extension as i32
        {
            return;
        }
        let alphabet_size_usize = alphabet_size as usize;
        let query_len_usize = query_len as usize;
        let target_len_usize = target_len as usize;
        let gap_open_i32 = gap_open as i32;
        let gap_extension_i32 = gap_extension as i32;
        let gap_open2_i32 = gap_open2 as i32;
        let non_canonical_penalty_i32 = non_canonical_penalty as i32;
        let gap_open_plus_extension = gap_open_i32 + gap_extension_i32;

        let zero_vector = _mm_setzero_si128();
        let gap_open_vector = _mm_set1_epi8(gap_open);
        let gap_open2_vector = _mm_set1_epi8(gap_open2);
        let gap_open_plus_extension_vector = _mm_set1_epi8(gap_open_plus_extension as i8);
        let match_score_vector = _mm_set1_epi8(score_matrix[0]);
        let mismatch_score_vector = _mm_set1_epi8(score_matrix[1]);
        let wildcard_score = score_matrix[alphabet_size_usize * alphabet_size_usize - 1];
        let ambiguous_score_vector = if wildcard_score == 0 {
            _mm_set1_epi8(-gap_extension)
        } else {
            _mm_set1_epi8(wildcard_score)
        };
        // wildcard symbol, i.e. the last row/column of the scoring matrix
        let wildcard_symbol_vector = _mm_set1_epi8((alphabet_size as i32 - 1) as i8);

        let target_block_count = ((target_len + 15) / 16) as usize;
        let max_backtrack_columns = if query_len < target_len {
            query_len
        } else {
            target_len
        };
        let backtrack_blocks = ((max_backtrack_columns + 15) / 16 + 1) as usize;
        let query_block_count = ((query_len + 15) / 16) as usize;
        let mut max_matrix_score = score_matrix[0] as i32;
        let mut min_matrix_score = score_matrix.get(1).copied().unwrap_or(0) as i32;
        for matrix_index in 1..alphabet_size_usize * alphabet_size_usize {
            let matrix_score = score_matrix[matrix_index] as i32;
            if matrix_score > max_matrix_score {
                max_matrix_score = matrix_score;
            }
            if matrix_score < min_matrix_score {
                min_matrix_score = matrix_score;
            }
        }
        let _ = max_matrix_score;
        if -min_matrix_score > 2 * gap_open_plus_extension {
            return; // otherwise, we won't see any mismatches
        }

        let mut long_gap_threshold = (gap_open2_i32 - gap_open_i32) / gap_extension_i32 - 1;
        if gap_open2_i32 > gap_open_i32 + gap_extension_i32 + long_gap_threshold * gap_extension_i32 {
            long_gap_threshold += 1;
        }
        let long_gap_difference =
            long_gap_threshold * gap_extension_i32 - (gap_open2_i32 - gap_open_i32);

        let record_cigar = (flags & KSW_EZ_SCORE_ONLY) == 0;
        let approximate_max = (flags & KSW_EZ_APPROX_MAX) != 0;
        let mut u = vec![zero_vector; target_block_count];
        let mut v = vec![zero_vector; target_block_count];
        let mut x = vec![zero_vector; target_block_count];
        let mut y = vec![zero_vector; target_block_count];
        let mut x2 = vec![zero_vector; target_block_count];
        let mut donor = vec![zero_vector; target_block_count];
        let mut acceptor = vec![zero_vector; target_block_count];
        let mut s = vec![zero_vector; target_block_count + 1];
        let mut padded_target = vec![0u8; target_block_count * 16 + 16];
        let mut reverse_query = vec![0u8; (query_block_count + 1) * 16];
        let mut score_table: Vec<i32> = if approximate_max {
            Vec::new()
        } else {
            vec![KSW_NEG_INF; target_block_count * 16 + 4]
        };
        let row_count = (query_len + target_len - 1) as usize;
        let mut backtrack_matrix: Vec<u8> = Vec::new();
        let mut band_start_per_row: Vec<i32> = Vec::new();
        let mut band_end_per_row: Vec<i32> = Vec::new();
        if record_cigar {
            backtrack_matrix = vec![0u8; (row_count * backtrack_blocks + 1) * 16];
            band_start_per_row = vec![0i32; row_count];
            band_end_per_row = vec![0i32; row_count];
        }

        for query_index in 0..query_len_usize {
            reverse_query[query_index] = query[query_len_usize - 1 - query_index];
        }
        padded_target[..target_len_usize].copy_from_slice(&target[..target_len_usize]);

        // this sets u, v, x, y (they share one array in the C code)
        {
            let initial_lane = _mm_set1_epi8((-gap_open_plus_extension) as i8);
            for block_index in 0..target_block_count {
                u[block_index] = initial_lane;
                v[block_index] = initial_lane;
                x[block_index] = initial_lane;
                y[block_index] = initial_lane;
            }
            let initial_second_lane = _mm_set1_epi8((-gap_open2_i32) as i8);
            for block_index in 0..target_block_count {
                x2[block_index] = initial_second_lane;
            }
        }

        let u_ptr = u.as_mut_ptr() as *mut i8;
        let v_ptr = v.as_mut_ptr() as *mut i8;
        let x_ptr = x.as_mut_ptr() as *mut i8;
        let x2_ptr = x2.as_mut_ptr() as *mut i8;
        let y_ptr = y.as_mut_ptr() as *mut i8;
        let s_ptr = s.as_mut_ptr() as *mut u8;
        let padded_target_ptr = padded_target.as_ptr();
        let reverse_query_ptr = reverse_query.as_ptr();

        // set the donor and acceptor arrays. TODO: this assumes 0/1/2/3 encoding!
        if (flags & (KSW_EZ_SPLICE_FOR | KSW_EZ_SPLICE_REV)) != 0 {
            let splice_forward = (flags & KSW_EZ_SPLICE_FOR) != 0;
            let splice_reverse = (flags & KSW_EZ_SPLICE_REV) != 0;
            let semi_canonical_cost: i32 = if (flags & KSW_EZ_SPLICE_FLANK) != 0 {
                -non_canonical_penalty_i32 / 2
            } else {
                0
            };
            {
                let initial_lane = _mm_set1_epi8((-non_canonical_penalty_i32) as i8);
                for block_index in 0..target_block_count {
                    donor[block_index] = initial_lane;
                    acceptor[block_index] = initial_lane;
                }
            }
            let donor_ptr = donor.as_mut_ptr() as *mut i8;
            let acceptor_ptr = acceptor.as_mut_ptr() as *mut i8;
            if (flags & KSW_EZ_REV_CIGAR) == 0 {
                for target_index in 0..(target_len as usize).saturating_sub(4) {
                    let mut canonical_type = 0i32; // 0=none, 1=GT/AG only, 2=GTr/yAG
                    if splice_forward && target[target_index + 1] == 2 && target[target_index + 2] == 3 {
                        canonical_type = 1; // GTr...
                    }
                    if splice_reverse && target[target_index + 1] == 1 && target[target_index + 2] == 3 {
                        canonical_type = 1; // CTr...
                    }
                    if canonical_type != 0 && (target[target_index + 3] == 0 || target[target_index + 3] == 2) {
                        canonical_type = 2;
                    }
                    if canonical_type != 0 {
                        *donor_ptr.add(target_index) =
                            if canonical_type == 2 { 0 } else { semi_canonical_cost as i8 };
                    }
                }
                if !junctions.is_empty() {
                    for target_index in 0..(target_len as usize).saturating_sub(1) {
                        if (splice_forward && (junctions[target_index + 1] & 1) != 0)
                            || (splice_reverse
                                && (junctions[target_index + 1] & 8) != 0)
                        {
                            *donor_ptr.add(target_index) =
                                donor_ptr.add(target_index).read().wrapping_add(junction_bonus);
                        }
                    }
                }
                for target_index in 2..target_len as usize {
                    let mut canonical_type = 0i32;
                    if splice_forward && target[target_index - 1] == 0 && target[target_index] == 2 {
                        canonical_type = 1; // ...yAG
                    }
                    if splice_reverse && target[target_index - 1] == 0 && target[target_index] == 1 {
                        canonical_type = 1; // ...yAC
                    }
                    if canonical_type != 0 && (target[target_index - 2] == 1 || target[target_index - 2] == 3) {
                        canonical_type = 2;
                    }
                    if canonical_type != 0 {
                        *acceptor_ptr.add(target_index) =
                            if canonical_type == 2 { 0 } else { semi_canonical_cost as i8 };
                    }
                }
                if !junctions.is_empty() {
                    for target_index in 0..target_len as usize {
                        if (splice_forward && (junctions[target_index] & 2) != 0)
                            || (splice_reverse && (junctions[target_index] & 4) != 0)
                        {
                            *acceptor_ptr.add(target_index) =
                                acceptor_ptr.add(target_index).read().wrapping_add(junction_bonus);
                        }
                    }
                }
            } else {
                for target_index in 0..(target_len as usize).saturating_sub(4) {
                    let mut canonical_type = 0i32;
                    if splice_forward && target[target_index + 1] == 2 && target[target_index + 2] == 0 {
                        canonical_type = 1; // GAy...
                    }
                    if splice_reverse && target[target_index + 1] == 1 && target[target_index + 2] == 0 {
                        canonical_type = 1; // CAy...
                    }
                    if canonical_type != 0 && (target[target_index + 3] == 1 || target[target_index + 3] == 3) {
                        canonical_type = 2;
                    }
                    if canonical_type != 0 {
                        *donor_ptr.add(target_index) =
                            if canonical_type == 2 { 0 } else { semi_canonical_cost as i8 };
                    }
                }
                if !junctions.is_empty() {
                    for target_index in 0..(target_len as usize).saturating_sub(1) {
                        if (splice_forward && (junctions[target_index + 1] & 2) != 0)
                            || (splice_reverse
                                && (junctions[target_index + 1] & 4) != 0)
                        {
                            *donor_ptr.add(target_index) =
                                donor_ptr.add(target_index).read().wrapping_add(junction_bonus);
                        }
                    }
                }
                for target_index in 2..target_len as usize {
                    let mut canonical_type = 0i32;
                    if splice_forward && target[target_index - 1] == 3 && target[target_index] == 2 {
                        canonical_type = 1; // ...rTG
                    }
                    if splice_reverse && target[target_index - 1] == 3 && target[target_index] == 1 {
                        canonical_type = 1; // ...rTC
                    }
                    if canonical_type != 0 && (target[target_index - 2] == 0 || target[target_index - 2] == 2) {
                        canonical_type = 2;
                    }
                    if canonical_type != 0 {
                        *acceptor_ptr.add(target_index) =
                            if canonical_type == 2 { 0 } else { semi_canonical_cost as i8 };
                    }
                }
                if !junctions.is_empty() {
                    for target_index in 0..target_len as usize {
                        if (splice_forward && (junctions[target_index] & 1) != 0)
                            || (splice_reverse && (junctions[target_index] & 8) != 0)
                        {
                            *acceptor_ptr.add(target_index) =
                                acceptor_ptr.add(target_index).read().wrapping_add(junction_bonus);
                        }
                    }
                }
            }
        }

        let mut previous_band_start = -1i32;
        let mut previous_band_end = -1i32;
        let mut diagonal_score = 0i32;
        let mut diagonal_score_target_index = 0i32;

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
            let unpadded_band_start = band_start;
            let unpadded_band_end = band_end;
            band_start = band_start / 16 * 16;
            band_end = (band_end + 16) / 16 * 16 - 1;
            // set boundary conditions
            let previous_x: i8;
            let previous_x2: i8;
            let previous_v: i8;
            if band_start > 0 {
                if band_start - 1 >= previous_band_start && band_start - 1 <= previous_band_end {
                    // (anti_diagonal-1, band_start-1) calculated in the last round
                    previous_x = *x_ptr.add((band_start - 1) as usize);
                    previous_x2 = *x2_ptr.add((band_start - 1) as usize);
                    previous_v = *v_ptr.add((band_start - 1) as usize);
                } else {
                    previous_x = (-gap_open_plus_extension) as i8;
                    previous_x2 = (-gap_open2_i32) as i8;
                    previous_v = (-gap_open_plus_extension) as i8;
                }
            } else {
                previous_x = (-gap_open_plus_extension) as i8;
                previous_x2 = (-gap_open2_i32) as i8;
                previous_v = if anti_diagonal == 0 {
                    (-gap_open_plus_extension) as i8
                } else if anti_diagonal < long_gap_threshold {
                    (-gap_extension_i32) as i8
                } else if anti_diagonal == long_gap_threshold {
                    long_gap_difference as i8
                } else {
                    0
                };
            }
            if band_end >= anti_diagonal {
                *y_ptr.add(anti_diagonal as usize) = (-gap_open_plus_extension) as i8;
                *u_ptr.add(anti_diagonal as usize) = if anti_diagonal == 0 {
                    (-gap_open_plus_extension) as i8
                } else if anti_diagonal < long_gap_threshold {
                    (-gap_extension_i32) as i8
                } else if anti_diagonal == long_gap_threshold {
                    long_gap_difference as i8
                } else {
                    0
                };
            }
            // loop fission: set scores first
            let reverse_query_base = query_len_usize as isize - 1 - anti_diagonal as isize;
            if (flags & KSW_EZ_GENERIC_SC) == 0 {
                let mut target_index = unpadded_band_start;
                while target_index <= unpadded_band_end {
                    let target_index_usize = target_index as usize;
                    let target_block =
                        _mm_loadu_si128(padded_target_ptr.add(target_index_usize) as *const __m128i);
                    let query_block = _mm_loadu_si128(
                        reverse_query_ptr.offset(reverse_query_base + target_index as isize)
                            as *const __m128i,
                    );
                    let wildcard_mask = _mm_or_si128(
                        _mm_cmpeq_epi8(target_block, wildcard_symbol_vector),
                        _mm_cmpeq_epi8(query_block, wildcard_symbol_vector),
                    );
                    let mut score_block = _mm_cmpeq_epi8(target_block, query_block);
                    score_block = _mm_blendv_epi8(mismatch_score_vector, match_score_vector, score_block);
                    score_block = _mm_blendv_epi8(score_block, ambiguous_score_vector, wildcard_mask);
                    _mm_storeu_si128(s_ptr.add(target_index_usize) as *mut __m128i, score_block);
                    target_index += 16;
                }
            } else {
                let mut target_index = unpadded_band_start;
                while target_index <= unpadded_band_end {
                    let target_index_usize = target_index as usize;
                    *s_ptr.add(target_index_usize) = score_matrix[*padded_target_ptr.add(target_index_usize)
                        as usize
                        * alphabet_size_usize
                        + *reverse_query_ptr.offset(reverse_query_base + target_index as isize) as usize]
                        as u8;
                    target_index += 1;
                }
            }
            // core loop
            let mut previous_x_vector = _mm_cvtsi32_si128((previous_x as u8) as i32); // C: (uint8_t)x1
            let mut previous_x2_vector = _mm_cvtsi32_si128((previous_x2 as u8) as i32); // C: (uint8_t)x21
            let mut previous_v_vector = _mm_cvtsi32_si128((previous_v as u8) as i32); // C: (uint8_t)v1
            let block_start = (band_start / 16) as usize;
            let block_end = (band_end / 16) as usize;
            assert!(block_end - block_start + 1 <= backtrack_blocks);
            if !record_cigar {
                // score only
                let mut block = block_start;
                while block <= block_end {
                    let mut max_score_vector = _mm_load_si128(s.as_ptr().add(block));
                    // x_previous_block <- x[r-1][t-1..t+14]
                    let mut x_previous_block = _mm_load_si128(x.as_ptr().add(block));
                    let mut temporary = _mm_srli_si128::<15>(x_previous_block); // x[r-1][t+15]
                    x_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x_previous_block), previous_x_vector);
                    previous_x_vector = temporary;
                    // v_previous_block <- v[r-1][t-1..t+14]
                    let mut v_previous_block = _mm_load_si128(v.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(v_previous_block); // v[r-1][t+15]
                    v_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(v_previous_block), previous_v_vector);
                    previous_v_vector = temporary;
                    // candidate_from_x <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]
                    let mut candidate_from_x = _mm_add_epi8(x_previous_block, v_previous_block);
                    let u_current_block = _mm_load_si128(u.as_ptr().add(block));
                    // candidate_from_y <- y[r-1][t..t+15] + u[r-1][t..t+15]
                    let mut candidate_from_y =
                        _mm_add_epi8(_mm_load_si128(y.as_ptr().add(block)), u_current_block);
                    let mut x2_previous_block = _mm_load_si128(x2.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(x2_previous_block);
                    x2_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x2_previous_block), previous_x2_vector);
                    previous_x2_vector = temporary;
                    let mut candidate_from_x2 = _mm_add_epi8(x2_previous_block, v_previous_block);
                    // only the acceptor site can start an intron
                    let candidate_from_x2_with_acceptor =
                        _mm_add_epi8(candidate_from_x2, _mm_load_si128(acceptor.as_ptr().add(block)));
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x);
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_y);
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x2_with_acceptor);
                    // save u[] and v[]; update the candidates
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
                    temporary = _mm_sub_epi8(max_score_vector, gap_open_vector);
                    candidate_from_x = _mm_sub_epi8(candidate_from_x, temporary);
                    candidate_from_y = _mm_sub_epi8(candidate_from_y, temporary);
                    candidate_from_x2 =
                        _mm_sub_epi8(candidate_from_x2, _mm_sub_epi8(max_score_vector, gap_open2_vector));
                    _mm_store_si128(
                        x.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_max_epi8(candidate_from_x, zero_vector), gap_open_plus_extension_vector),
                    );
                    _mm_store_si128(
                        y.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_max_epi8(candidate_from_y, zero_vector), gap_open_plus_extension_vector),
                    );
                    temporary = _mm_load_si128(donor.as_ptr().add(block));
                    // the donor site bounds the value inherited by the next row
                    _mm_store_si128(
                        x2.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_max_epi8(candidate_from_x2, temporary), gap_open2_vector),
                    );
                    block += 1;
                }
            } else if (flags & KSW_EZ_RIGHT) == 0 {
                // gap left-alignment
                let backtrack_row_offset =
                    anti_diagonal as usize * backtrack_blocks * 16 - band_start as usize;
                band_start_per_row[anti_diagonal as usize] = band_start;
                band_end_per_row[anti_diagonal as usize] = band_end;
                let mut block = block_start;
                while block <= block_end {
                    let mut max_score_vector = _mm_load_si128(s.as_ptr().add(block));
                    let mut x_previous_block = _mm_load_si128(x.as_ptr().add(block));
                    let mut temporary = _mm_srli_si128::<15>(x_previous_block);
                    x_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x_previous_block), previous_x_vector);
                    previous_x_vector = temporary;
                    let mut v_previous_block = _mm_load_si128(v.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(v_previous_block);
                    v_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(v_previous_block), previous_v_vector);
                    previous_v_vector = temporary;
                    let mut candidate_from_x = _mm_add_epi8(x_previous_block, v_previous_block);
                    let u_current_block = _mm_load_si128(u.as_ptr().add(block));
                    let mut candidate_from_y =
                        _mm_add_epi8(_mm_load_si128(y.as_ptr().add(block)), u_current_block);
                    let mut x2_previous_block = _mm_load_si128(x2.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(x2_previous_block);
                    x2_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x2_previous_block), previous_x2_vector);
                    previous_x2_vector = temporary;
                    let mut candidate_from_x2 = _mm_add_epi8(x2_previous_block, v_previous_block);
                    let candidate_from_x2_with_acceptor =
                        _mm_add_epi8(candidate_from_x2, _mm_load_si128(acceptor.as_ptr().add(block)));

                    let mut backtrack_byte_vector = _mm_and_si128(
                        _mm_cmpgt_epi8(candidate_from_x, max_score_vector),
                        _mm_set1_epi8(DELETION_STATE),
                    ); // 1 when candidate_from_x > z
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x);
                    backtrack_byte_vector = _mm_blendv_epi8(
                        backtrack_byte_vector,
                        _mm_set1_epi8(INSERTION_STATE),
                        _mm_cmpgt_epi8(candidate_from_y, max_score_vector),
                    ); // 2 when candidate_from_y > z
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_y);
                    backtrack_byte_vector = _mm_blendv_epi8(
                        backtrack_byte_vector,
                        _mm_set1_epi8(INTRON_STATE),
                        _mm_cmpgt_epi8(candidate_from_x2_with_acceptor, max_score_vector),
                    ); // 3 when candidate_from_x2_with_acceptor > z
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x2_with_acceptor);
                    // save u[] and v[]; update the candidates
                    _mm_store_si128(
                        u.as_mut_ptr().add(block),
                        _mm_sub_epi8(max_score_vector, v_previous_block),
                    );
                    _mm_store_si128(
                        v.as_mut_ptr().add(block),
                        _mm_sub_epi8(max_score_vector, u_current_block),
                    );
                    temporary = _mm_sub_epi8(max_score_vector, gap_open_vector);
                    candidate_from_x = _mm_sub_epi8(candidate_from_x, temporary);
                    candidate_from_y = _mm_sub_epi8(candidate_from_y, temporary);
                    candidate_from_x2 =
                        _mm_sub_epi8(candidate_from_x2, _mm_sub_epi8(max_score_vector, gap_open2_vector));

                    temporary = _mm_cmpgt_epi8(candidate_from_x, zero_vector);
                    _mm_store_si128(
                        x.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_and_si128(temporary, candidate_from_x), gap_open_plus_extension_vector),
                    );
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_and_si128(temporary, _mm_set1_epi8(DELETION_CONTINUATION)),
                    ); // the E state continues while candidate_from_x > 0
                    temporary = _mm_cmpgt_epi8(candidate_from_y, zero_vector);
                    _mm_store_si128(
                        y.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_and_si128(temporary, candidate_from_y), gap_open_plus_extension_vector),
                    );
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_and_si128(temporary, _mm_set1_epi8(INSERTION_CONTINUATION)),
                    ); // the F state continues while candidate_from_y > 0

                    let donor_block = _mm_load_si128(donor.as_ptr().add(block));
                    temporary = _mm_cmpgt_epi8(candidate_from_x2, donor_block);
                    let junction_score_vector = _mm_max_epi8(candidate_from_x2, donor_block);
                    _mm_store_si128(
                        x2.as_mut_ptr().add(block),
                        _mm_sub_epi8(junction_score_vector, gap_open2_vector),
                    );
                    // the E2 state continues while candidate_from_x2 > donor
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_and_si128(temporary, _mm_set1_epi8(INTRON_CONTINUATION)),
                    );
                    _mm_storeu_si128(
                        backtrack_matrix
                            .as_mut_ptr()
                            .add(backtrack_row_offset + block * 16) as *mut __m128i,
                        backtrack_byte_vector,
                    );
                    block += 1;
                }
            } else {
                // gap right-alignment
                let backtrack_row_offset =
                    anti_diagonal as usize * backtrack_blocks * 16 - band_start as usize;
                band_start_per_row[anti_diagonal as usize] = band_start;
                band_end_per_row[anti_diagonal as usize] = band_end;
                let mut block = block_start;
                while block <= block_end {
                    let mut max_score_vector = _mm_load_si128(s.as_ptr().add(block));
                    let mut x_previous_block = _mm_load_si128(x.as_ptr().add(block));
                    let mut temporary = _mm_srli_si128::<15>(x_previous_block);
                    x_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x_previous_block), previous_x_vector);
                    previous_x_vector = temporary;
                    let mut v_previous_block = _mm_load_si128(v.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(v_previous_block);
                    v_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(v_previous_block), previous_v_vector);
                    previous_v_vector = temporary;
                    let mut candidate_from_x = _mm_add_epi8(x_previous_block, v_previous_block);
                    let u_current_block = _mm_load_si128(u.as_ptr().add(block));
                    let mut candidate_from_y =
                        _mm_add_epi8(_mm_load_si128(y.as_ptr().add(block)), u_current_block);
                    let mut x2_previous_block = _mm_load_si128(x2.as_ptr().add(block));
                    temporary = _mm_srli_si128::<15>(x2_previous_block);
                    x2_previous_block =
                        _mm_or_si128(_mm_slli_si128::<1>(x2_previous_block), previous_x2_vector);
                    previous_x2_vector = temporary;
                    let mut candidate_from_x2 = _mm_add_epi8(x2_previous_block, v_previous_block);
                    let candidate_from_x2_with_acceptor =
                        _mm_add_epi8(candidate_from_x2, _mm_load_si128(acceptor.as_ptr().add(block)));

                    // 1 when candidate_from_x <= z
                    let mut backtrack_byte_vector = _mm_andnot_si128(
                        _mm_cmpgt_epi8(max_score_vector, candidate_from_x),
                        _mm_set1_epi8(DELETION_STATE),
                    );
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x);
                    backtrack_byte_vector = _mm_blendv_epi8(
                        _mm_set1_epi8(INSERTION_STATE),
                        backtrack_byte_vector,
                        _mm_cmpgt_epi8(max_score_vector, candidate_from_y),
                    ); // 2 unless candidate_from_y > z
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_y);
                    backtrack_byte_vector = _mm_blendv_epi8(
                        _mm_set1_epi8(INTRON_STATE),
                        backtrack_byte_vector,
                        _mm_cmpgt_epi8(max_score_vector, candidate_from_x2_with_acceptor),
                    ); // 3 unless candidate_from_x2_with_acceptor > z
                    max_score_vector = _mm_max_epi8(max_score_vector, candidate_from_x2_with_acceptor);
                    // save u[] and v[]; update the candidates
                    _mm_store_si128(
                        u.as_mut_ptr().add(block),
                        _mm_sub_epi8(max_score_vector, v_previous_block),
                    );
                    _mm_store_si128(
                        v.as_mut_ptr().add(block),
                        _mm_sub_epi8(max_score_vector, u_current_block),
                    );
                    temporary = _mm_sub_epi8(max_score_vector, gap_open_vector);
                    candidate_from_x = _mm_sub_epi8(candidate_from_x, temporary);
                    candidate_from_y = _mm_sub_epi8(candidate_from_y, temporary);
                    candidate_from_x2 =
                        _mm_sub_epi8(candidate_from_x2, _mm_sub_epi8(max_score_vector, gap_open2_vector));

                    temporary = _mm_cmpgt_epi8(zero_vector, candidate_from_x);
                    _mm_store_si128(
                        x.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_andnot_si128(temporary, candidate_from_x), gap_open_plus_extension_vector),
                    );
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_andnot_si128(temporary, _mm_set1_epi8(DELETION_CONTINUATION)),
                    ); // the E state continues while candidate_from_x > 0
                    temporary = _mm_cmpgt_epi8(zero_vector, candidate_from_y);
                    _mm_store_si128(
                        y.as_mut_ptr().add(block),
                        _mm_sub_epi8(_mm_andnot_si128(temporary, candidate_from_y), gap_open_plus_extension_vector),
                    );
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_andnot_si128(temporary, _mm_set1_epi8(INSERTION_CONTINUATION)),
                    ); // the F state continues while candidate_from_y > 0

                    let donor_block = _mm_load_si128(donor.as_ptr().add(block));
                    temporary = _mm_cmpgt_epi8(donor_block, candidate_from_x2);
                    let junction_score_vector = _mm_max_epi8(donor_block, candidate_from_x2);
                    _mm_store_si128(
                        x2.as_mut_ptr().add(block),
                        _mm_sub_epi8(junction_score_vector, gap_open2_vector),
                    );
                    backtrack_byte_vector = _mm_or_si128(
                        backtrack_byte_vector,
                        _mm_andnot_si128(temporary, _mm_set1_epi8(INTRON_CONTINUATION)),
                    );
                    _mm_storeu_si128(
                        backtrack_matrix
                            .as_mut_ptr()
                            .add(backtrack_row_offset + block * 16) as *mut __m128i,
                        backtrack_byte_vector,
                    );
                    block += 1;
                }
            }

            if !approximate_max {
                // find the exact max with a 32-bit score array
                let exact_max_score: i32;
                let exact_max_score_target_index: i32;
                if anti_diagonal > 0 {
                    let unpadded_band_end_usize = unpadded_band_end as usize;
                    let mut max_score = if unpadded_band_end > 0 {
                        score_table[unpadded_band_end_usize - 1] + *u_ptr.add(unpadded_band_end_usize) as i32
                    } else {
                        score_table[unpadded_band_end_usize] + *v_ptr.add(unpadded_band_end_usize) as i32
                    };
                    score_table[unpadded_band_end_usize] = max_score; // special casing the last element
                    let mut max_score_target_index = unpadded_band_end;
                    let mut max_score_vector = _mm_set1_epi32(max_score);
                    let mut max_score_target_index_vector = _mm_set1_epi32(max_score_target_index);
                    let vectorized_end = unpadded_band_start + (unpadded_band_end - unpadded_band_start) / 4 * 4;
                    let mut target_index = unpadded_band_start;
                    while target_index < vectorized_end {
                        let target_index_usize = target_index as usize;
                        let mut score_block =
                            _mm_loadu_si128(score_table.as_ptr().add(target_index_usize) as *const __m128i);
                        let v_scores = _mm_setr_epi32(
                            *v_ptr.add(target_index_usize) as i32,
                            *v_ptr.add(target_index_usize + 1) as i32,
                            *v_ptr.add(target_index_usize + 2) as i32,
                            *v_ptr.add(target_index_usize + 3) as i32,
                        );
                        score_block = _mm_add_epi32(score_block, v_scores);
                        _mm_storeu_si128(
                            score_table.as_mut_ptr().add(target_index_usize) as *mut __m128i,
                            score_block,
                        );
                        let target_index_vector = _mm_set1_epi32(target_index);
                        let temporary = _mm_cmpgt_epi32(score_block, max_score_vector);
                        max_score_vector = _mm_blendv_epi8(max_score_vector, score_block, temporary);
                        max_score_target_index_vector =
                            _mm_blendv_epi8(max_score_target_index_vector, target_index_vector, temporary);
                        target_index += 4;
                    }
                    let mut lane_scores = [0i32; 4];
                    let mut lane_indices = [0i32; 4];
                    _mm_storeu_si128(lane_scores.as_mut_ptr() as *mut __m128i, max_score_vector);
                    _mm_storeu_si128(
                        lane_indices.as_mut_ptr() as *mut __m128i,
                        max_score_target_index_vector,
                    );
                    for lane in 0..4 {
                        if max_score < lane_scores[lane] {
                            max_score = lane_scores[lane];
                            max_score_target_index = lane_indices[lane] + lane as i32;
                        }
                    }
                    while target_index < unpadded_band_end {
                        let target_index_usize = target_index as usize;
                        score_table[target_index_usize] += *v_ptr.add(target_index_usize) as i32;
                        if score_table[target_index_usize] > max_score {
                            max_score = score_table[target_index_usize];
                            max_score_target_index = target_index;
                        }
                        target_index += 1;
                    }
                    exact_max_score = max_score;
                    exact_max_score_target_index = max_score_target_index;
                } else {
                    score_table[0] = *v_ptr as i32 - gap_open_plus_extension; // special casing anti_diagonal == 0
                    exact_max_score = score_table[0];
                    exact_max_score_target_index = 0;
                }
                // update the result
                if unpadded_band_end == target_len - 1 && score_table[unpadded_band_end as usize] > result.mte {
                    result.mte = score_table[unpadded_band_end as usize];
                    result.mte_q = anti_diagonal - band_end;
                }
                if anti_diagonal - unpadded_band_start == query_len - 1
                    && score_table[unpadded_band_start as usize] > result.mqe
                {
                    result.mqe = score_table[unpadded_band_start as usize];
                    result.mqe_t = unpadded_band_start;
                }
                if apply_z_drop(
                    result,
                    true,
                    exact_max_score,
                    anti_diagonal,
                    exact_max_score_target_index,
                    z_drop,
                    0,
                ) {
                    break;
                }
                if anti_diagonal == query_len + target_len - 2 && unpadded_band_end == target_len - 1 {
                    result.score = score_table[target_len_usize - 1];
                }
            } else {
                // find approximate max; Z-drop might be inaccurate, too.
                if anti_diagonal > 0 {
                    if diagonal_score_target_index >= unpadded_band_start
                        && diagonal_score_target_index <= unpadded_band_end
                        && diagonal_score_target_index + 1 >= unpadded_band_start
                        && diagonal_score_target_index + 1 <= unpadded_band_end
                    {
                        let score_from_v = *v_ptr.add(diagonal_score_target_index as usize) as i32;
                        let score_from_u = *u_ptr.add(diagonal_score_target_index as usize + 1) as i32;
                        if score_from_v > score_from_u {
                            diagonal_score += score_from_v;
                        } else {
                            diagonal_score += score_from_u;
                            diagonal_score_target_index += 1;
                        }
                    } else if diagonal_score_target_index >= unpadded_band_start
                        && diagonal_score_target_index <= unpadded_band_end
                    {
                        diagonal_score += *v_ptr.add(diagonal_score_target_index as usize) as i32;
                    } else {
                        diagonal_score_target_index += 1;
                        diagonal_score += *u_ptr.add(diagonal_score_target_index as usize) as i32;
                    }
                } else {
                    diagonal_score = *v_ptr as i32 - gap_open_plus_extension;
                    diagonal_score_target_index = 0;
                }
                // note: unlike extz2, the C code calls ksw_apply_zdrop() outside
                // the "if (anti_diagonal > 0)" block here
                if (flags & KSW_EZ_APPROX_DROP) != 0
                    && apply_z_drop(
                        result,
                        true,
                        diagonal_score,
                        anti_diagonal,
                        diagonal_score_target_index,
                        z_drop,
                        0,
                    )
                {
                    break;
                }
                if anti_diagonal == query_len + target_len - 2 && unpadded_band_end == target_len - 1 {
                    result.score = diagonal_score;
                }
            }
            previous_band_start = band_start;
            previous_band_end = band_end;
        }

        if record_cigar {
            let reverse_cigar = (flags & KSW_EZ_REV_CIGAR) != 0;
            let backtrack_bytes = &backtrack_matrix[..row_count * backtrack_blocks * 16];
            let backtrack_row_bytes = backtrack_blocks as i32 * 16;
            if !result.zdropped && (flags & KSW_EZ_EXTZ_ONLY) == 0 {
                backtrack(
                    true,
                    reverse_cigar,
                    long_gap_threshold,
                    backtrack_bytes,
                    &band_start_per_row,
                    Some(&band_end_per_row),
                    backtrack_row_bytes,
                    target_len - 1,
                    query_len - 1,
                    &mut result.cigar,
                );
            } else if result.max_t >= 0 && result.max_q >= 0 {
                backtrack(
                    true,
                    reverse_cigar,
                    long_gap_threshold,
                    backtrack_bytes,
                    &band_start_per_row,
                    Some(&band_end_per_row),
                    backtrack_row_bytes,
                    result.max_t,
                    result.max_q,
                    &mut result.cigar,
                );
            }
        }
    }
}
