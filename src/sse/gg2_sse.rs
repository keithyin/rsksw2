//! Port of `ksw2_gg2_sse.c` (SSE4.1 path).

use core::arch::x86_64::*;

use crate::common::backtrack;

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_gg2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    w: i32,
    cigar: &mut Vec<u32>,
) -> i32 {
    unsafe {
        let mu = m as usize;
        let q_i = q as i32;
        let qe = q_i + e as i32;
        let qlen_u = qlen as usize;
        let tlen_u = tlen as usize;

        let zero_ = _mm_setzero_si128();
        let q_ = _mm_set1_epi8(q);
        let qe2_ = _mm_set1_epi8((qe * 2) as i8);
        let flag1_ = _mm_set1_epi8(1);
        let flag2_ = _mm_set1_epi8(2);
        let flag8_ = _mm_set1_epi8(0x08);
        let flag16_ = _mm_set1_epi8(0x10);

        let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };
        let n_col0 = if w + 1 < tlen { w + 1 } else { tlen }; // number of columns in the backtrack matrix
        let tlen_ = ((tlen + 15) / 16) as usize;
        let n_col_ = ((n_col0 + 15) / 16 + 1) as usize;
        let n_col = (n_col_ * 16) as i32;

        let mut u = vec![zero_; tlen_];
        let mut v = vec![zero_; tlen_];
        let mut x = vec![zero_; tlen_];
        let mut y = vec![zero_; tlen_];
        let mut s = vec![zero_; tlen_];
        let mut sf = vec![0u8; tlen_ * 16];
        let mut qr = vec![0u8; qlen_u];
        let n_rows = (qlen + tlen - 1) as usize;
        let mut p = vec![0u8; (n_rows * n_col_ + 1) * 16];
        let mut off = vec![0i32; n_rows];

        for t in 0..qlen_u {
            qr[t] = query[qlen_u - 1 - t];
        }
        sf[..tlen_u].copy_from_slice(&target[..tlen_u]);

        let up = u.as_mut_ptr() as *mut u8;
        let vp = v.as_mut_ptr() as *mut u8;
        let sp = s.as_mut_ptr() as *mut u8;
        let sfp = sf.as_ptr();
        let qrp = qr.as_ptr();
        let matp = mat.as_ptr();

        let mut h0 = 0i32;
        let mut last_h0_t = 0i32;
        let mut last_st = -1i32;
        let mut last_en = -1i32;

        for r in 0..qlen + tlen - 1 {
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
            let st0 = st;
            let en0 = en;
            st = st / 16 * 16;
            en = (en + 16) / 16 * 16 - 1;
            off[r as usize] = st;
            // set boundary conditions
            let x1: i8;
            let v1: i8;
            if st > 0 {
                if st - 1 >= last_st && st - 1 <= last_en {
                    x1 = *(x.as_ptr() as *const i8).add((st - 1) as usize); // (r-1,s-1) calculated in the last round
                    v1 = *vp.add((st - 1) as usize) as i8;
                } else {
                    x1 = 0;
                    v1 = 0; // not calculated; set to zeros
                }
            } else {
                x1 = 0;
                v1 = if r != 0 { q } else { 0 };
            }
            if en >= r {
                *(y.as_ptr() as *mut u8).add(r as usize) = 0;
                *up.add(r as usize) = if r != 0 { q as u8 } else { 0 };
            }
            // loop fission: set scores first
            let qrr_base = qlen_u as isize - 1 - r as isize;
            let mut t = st0;
            while t <= en0 {
                let tu = t as usize;
                *sp.add(tu) = *matp.add(
                    *sfp.add(tu) as usize * mu + *qrp.offset(qrr_base + t as isize) as usize,
                ) as u8;
                t += 1;
            }
            // core loop
            let mut x1_ = _mm_cvtsi32_si128(x1 as i32);
            let mut v1_ = _mm_cvtsi32_si128(v1 as i32);
            let st_ = (st / 16) as usize;
            let en_ = (en / 16) as usize;
            let pr_off = r as usize * n_col_ * 16 - st as usize;
            let mut t = st_;
            while t <= en_ {
                let mut z = _mm_add_epi8(_mm_load_si128(s.as_ptr().add(t)), qe2_);

                let mut xt1 = _mm_load_si128(x.as_ptr().add(t)); // xt1 <- x[r-1][t..t+15]
                let mut tmp = _mm_srli_si128::<15>(xt1); // tmp <- x[r-1][t+15]
                xt1 = _mm_or_si128(_mm_slli_si128::<1>(xt1), x1_); // xt1 <- x[r-1][t-1..t+14]
                x1_ = tmp;
                let mut vt1 = _mm_load_si128(v.as_ptr().add(t)); // vt1 <- v[r-1][t..t+15]
                tmp = _mm_srli_si128::<15>(vt1); // tmp <- v[r-1][t+15]
                vt1 = _mm_or_si128(_mm_slli_si128::<1>(vt1), v1_); // vt1 <- v[r-1][t-1..t+14]
                v1_ = tmp;
                let mut a = _mm_add_epi8(xt1, vt1); // a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]

                let ut = _mm_load_si128(u.as_ptr().add(t)); // ut <- u[t..t+15]
                let mut b = _mm_add_epi8(_mm_load_si128(y.as_ptr().add(t)), ut); // b <- y[r-1][t..t+15] + u[r-1][t..t+15]

                let mut d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
                z = _mm_max_epi8(z, a);
                tmp = _mm_cmpgt_epi8(b, z);
                d = _mm_blendv_epi8(d, flag2_, tmp); // d = b > z? 2 : d
                z = _mm_max_epu8(z, b); // z = max(z, b); this works because both are non-negative
                _mm_store_si128(u.as_mut_ptr().add(t), _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
                _mm_store_si128(v.as_mut_ptr().add(t), _mm_sub_epi8(z, ut)); // v[r][t..t+15] <- z - u[r-1][t..t+15]

                z = _mm_sub_epi8(z, q_);
                a = _mm_sub_epi8(a, z);
                b = _mm_sub_epi8(b, z);
                tmp = _mm_cmpgt_epi8(a, zero_);
                d = _mm_or_si128(d, _mm_and_si128(flag8_, tmp));
                _mm_store_si128(x.as_mut_ptr().add(t), _mm_and_si128(a, tmp));
                tmp = _mm_cmpgt_epi8(b, zero_);
                d = _mm_or_si128(d, _mm_and_si128(flag16_, tmp));
                _mm_store_si128(y.as_mut_ptr().add(t), _mm_and_si128(b, tmp));
                _mm_storeu_si128(p.as_mut_ptr().add(pr_off + t * 16) as *mut __m128i, d);
                t += 1;
            }
            if r > 0 {
                if last_h0_t >= st0 && last_h0_t <= en0 {
                    h0 += *vp.add(last_h0_t as usize) as i32 - qe;
                } else {
                    last_h0_t += 1;
                    h0 += *up.add(last_h0_t as usize) as i32 - qe;
                }
            } else {
                h0 = *vp as i32 - 2 * qe;
                last_h0_t = 0;
            }
            last_st = st;
            last_en = en;
        }

        backtrack(
            true,
            false,
            0,
            &p[..n_rows * n_col_ * 16],
            &off,
            None,
            n_col,
            tlen - 1,
            qlen - 1,
            cigar,
        );
        h0
    }
}
