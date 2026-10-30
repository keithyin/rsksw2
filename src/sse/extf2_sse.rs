//! Port of `ksw2_extf2_sse.c` (SSE4.1 path): score-only extension with a
//! linear gap cost and X-drop.

use core::arch::x86_64::*;

use crate::common::reset_extz;
use crate::KswExtz;

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_extf2_sse(
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
    unsafe {
        reset_extz(ez);
        if qlen <= 0 || tlen <= 0 {
            return;
        }
        let qlen_u = qlen as usize;
        let tlen_u = tlen as usize;
        let e_i = e as i32;

        let e2_ = _mm_set1_epi8((e_i * 2) as i8);
        let sc_mch_ = _mm_set1_epi8(mch);
        let sc_mis_ = _mm_set1_epi8(if mis < 0 { mis } else { -mis });
        let tlen_ = ((tlen + 15) / 16) as usize;
        let qlen_ = ((qlen + 15) / 16) as usize;
        let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };

        let zero_ = _mm_setzero_si128();
        let mut u = vec![zero_; tlen_];
        let mut v = vec![zero_; tlen_];
        let mut s = vec![zero_; tlen_ + 1];
        let mut sf = vec![0u8; tlen_ * 16 + 16];
        let mut qr = vec![0u8; (qlen_ + 1) * 16];

        for t in 0..qlen_u {
            qr[t] = query[qlen_u - 1 - t];
        }
        sf[..tlen_u].copy_from_slice(&target[..tlen_u]);

        let up = u.as_mut_ptr() as *mut u8;
        let vp = v.as_mut_ptr() as *mut u8;

        let mut h0 = 0i32;
        let mut last_h0_t = 0i32;
        let mut last_st = -1i32;
        let mut last_en = -1i32;
        let mut r = 0i32;
        while r < qlen + tlen - 1 {
            let mut st = 0i32;
            let mut en = tlen - 1;
            // find the boundaries
            if st < r - qlen + 1 {
                st = r - qlen + 1;
            }
            if en > r {
                en = r;
            }
            if st < (r - w + 1) >> 1 {
                st = (r - w + 1) >> 1; // take the ceil
            }
            if en > (r + w) >> 1 {
                en = (r + w) >> 1; // take the floor
            }
            if st > en {
                break;
            }
            let st0 = st;
            let en0 = en;
            st = st / 16 * 16;
            en = (en + 16) / 16 * 16 - 1;
            // set boundary conditions
            let v1: u8 = if st > 0 && st - 1 >= last_st && st - 1 <= last_en {
                *vp.add((st - 1) as usize)
            } else {
                0
            };
            if en >= r {
                *up.add(r as usize) = 0;
            }
            // core loop
            let mut v1_ = _mm_cvtsi32_si128(v1 as i32);
            let st_ = (st / 16) as usize;
            let en_ = (en / 16) as usize;
            let qrr_base = qlen_u as isize - 1 - r as isize;
            let mut t = st0;
            while t <= en0 {
                let sq = _mm_loadu_si128(sf.as_ptr().add(t as usize) as *const __m128i);
                let stq = _mm_loadu_si128(
                    qr.as_ptr().offset(qrr_base + t as isize) as *const __m128i
                );
                let mut tmp = _mm_cmpeq_epi8(sq, stq);
                tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
                _mm_storeu_si128(
                    (s.as_mut_ptr() as *mut u8).add(t as usize) as *mut __m128i,
                    tmp,
                );
                t += 16;
            }
            let mut t = st_;
            while t <= en_ {
                let mut z = _mm_add_epi8(_mm_load_si128(s.as_ptr().add(t)), e2_);
                let mut vt1 = _mm_load_si128(v.as_ptr().add(t)); // vt1 <- v[r-1][t..t+15]
                let tmp = _mm_srli_si128::<15>(vt1); // tmp <- v[r-1][t+15]
                vt1 = _mm_or_si128(_mm_slli_si128::<1>(vt1), v1_); // vt1 <- v[r-1][t-1..t+14]
                v1_ = tmp;
                let ut = _mm_load_si128(u.as_ptr().add(t)); // ut <- u[t..t+15]
                z = _mm_max_epi8(z, vt1);
                z = _mm_max_epu8(z, ut); // z = max(z, b); this works because both are non-negative
                _mm_store_si128(u.as_mut_ptr().add(t), _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
                _mm_store_si128(v.as_mut_ptr().add(t), _mm_sub_epi8(z, ut)); // v[r][t..t+15] <- z - u[r-1][t..t+15]
                t += 1;
            }
            if r > 0 {
                if last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t + 1 <= en0
                {
                    let d0 = *vp.add(last_h0_t as usize) as i32 - e_i;
                    let d1 = *up.add(last_h0_t as usize + 1) as i32 - e_i;
                    if d0 > d1 {
                        h0 += d0;
                    } else {
                        h0 += d1;
                        last_h0_t += 1;
                    }
                } else if last_h0_t >= st0 && last_h0_t <= en0 {
                    h0 += *vp.add(last_h0_t as usize) as i32 - e_i;
                } else {
                    last_h0_t += 1;
                    h0 += *up.add(last_h0_t as usize) as i32 - e_i;
                }
                if h0 > ez.max as i32 {
                    ez.max = (h0 as u32) & 0x7fff_ffff;
                    ez.max_t = last_h0_t;
                    ez.max_q = r - last_h0_t;
                } else if xdrop >= 0 && ez.max as i32 - h0 > xdrop {
                    break;
                }
            } else {
                h0 = *vp as i32 - e_i - e_i;
                last_h0_t = 0;
            }
            last_st = st;
            last_en = en;
            r += 1;
        }
        if r == qlen + tlen - 1 {
            ez.score = h0;
        } else {
            ez.zdropped = true;
        }
    }
}
