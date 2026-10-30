//! Port of `ksw2_exts2_sse.c` (SSE4.1 path): splice-aware extension.

use core::arch::x86_64::*;

use crate::common::{apply_zdrop, backtrack, reset_extz};
use crate::{
    KSW_EZ_APPROX_DROP, KSW_EZ_APPROX_MAX, KSW_EZ_EXTZ_ONLY, KSW_EZ_GENERIC_SC,
    KSW_EZ_REV_CIGAR, KSW_EZ_RIGHT,
    KSW_EZ_SCORE_ONLY, KSW_EZ_SPLICE_FLANK, KSW_EZ_SPLICE_FOR, KSW_EZ_SPLICE_REV, KSW_NEG_INF,
    KswExtz,
};

#[target_feature(enable = "sse4.1")]
pub(crate) unsafe fn ksw_exts2_sse(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    q2: i8,
    noncan: i8,
    zdrop: i32,
    junc_bonus: i8,
    flag: i32,
    junc: &[u8],
    ez: &mut KswExtz,
) {
    unsafe {
        reset_extz(ez);
        if m <= 1 || qlen <= 0 || tlen <= 0 || q2 as i32 <= q as i32 + e as i32 {
            return;
        }
        let mu = m as usize;
        let qlen_u = qlen as usize;
        let tlen_u = tlen as usize;
        let q_i = q as i32;
        let e_i = e as i32;
        let q2_i = q2 as i32;
        let noncan_i = noncan as i32;
        let qe = q_i + e_i;

        let zero_ = _mm_setzero_si128();
        let q_ = _mm_set1_epi8(q);
        let q2_ = _mm_set1_epi8(q2);
        let qe_ = _mm_set1_epi8(qe as i8);
        let sc_mch_ = _mm_set1_epi8(mat[0]);
        let sc_mis_ = _mm_set1_epi8(mat[1]);
        let sc_n = mat[mu * mu - 1];
        let sc_amb_ = if sc_n == 0 { _mm_set1_epi8(-e) } else { _mm_set1_epi8(sc_n) };
        let m1_ = _mm_set1_epi8((m as i32 - 1) as i8); // wildcard

        let tlen_ = ((tlen + 15) / 16) as usize;
        let n_col0 = if qlen < tlen { qlen } else { tlen };
        let n_col_ = ((n_col0 + 15) / 16 + 1) as usize;
        let qlen_ = ((qlen + 15) / 16) as usize;
        let mut max_sc = mat[0] as i32;
        let mut min_sc = mat.get(1).copied().unwrap_or(0) as i32;
        for t in 1..mu * mu {
            let vv = mat[t] as i32;
            if vv > max_sc {
                max_sc = vv;
            }
            if vv < min_sc {
                min_sc = vv;
            }
        }
        let _ = max_sc;
        if -min_sc > 2 * qe {
            return; // otherwise, we won't see any mismatches
        }

        let mut long_thres = (q2_i - q_i) / e_i - 1;
        if q2_i > q_i + e_i + long_thres * e_i {
            long_thres += 1;
        }
        let long_diff = long_thres * e_i - (q2_i - q_i);

        let with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0;
        let approx_max = (flag & KSW_EZ_APPROX_MAX) != 0;
        let mut u = vec![zero_; tlen_];
        let mut v = vec![zero_; tlen_];
        let mut x = vec![zero_; tlen_];
        let mut y = vec![zero_; tlen_];
        let mut x2 = vec![zero_; tlen_];
        let mut donor = vec![zero_; tlen_];
        let mut acceptor = vec![zero_; tlen_];
        let mut s = vec![zero_; tlen_ + 1];
        let mut sf = vec![0u8; tlen_ * 16 + 16];
        let mut qr = vec![0u8; (qlen_ + 1) * 16];
        let mut h: Vec<i32> = if approx_max { Vec::new() } else { vec![KSW_NEG_INF; tlen_ * 16 + 4] };
        let n_rows = (qlen + tlen - 1) as usize;
        let mut p: Vec<u8> = Vec::new();
        let mut off: Vec<i32> = Vec::new();
        let mut off_end: Vec<i32> = Vec::new();
        if with_cigar {
            p = vec![0u8; (n_rows * n_col_ + 1) * 16];
            off = vec![0i32; n_rows];
            off_end = vec![0i32; n_rows];
        }

        for t in 0..qlen_u {
            qr[t] = query[qlen_u - 1 - t];
        }
        sf[..tlen_u].copy_from_slice(&target[..tlen_u]);

        // this set u, v, x, y (because they are in the same array in the C code)
        {
            let lane = _mm_set1_epi8((-qe) as i8);
            for i in 0..tlen_ {
                u[i] = lane;
                v[i] = lane;
                x[i] = lane;
                y[i] = lane;
            }
            let lane2 = _mm_set1_epi8((-q2_i) as i8);
            for i in 0..tlen_ {
                x2[i] = lane2;
            }
        }

        let up = u.as_mut_ptr() as *mut i8;
        let vp = v.as_mut_ptr() as *mut i8;
        let xp = x.as_mut_ptr() as *mut i8;
        let x2p = x2.as_mut_ptr() as *mut i8;
        let yp = y.as_mut_ptr() as *mut i8;
        let sp = s.as_mut_ptr() as *mut u8;
        let sfp = sf.as_ptr();
        let qrp = qr.as_ptr();

        // set the donor and acceptor arrays. TODO: this assumes 0/1/2/3 encoding!
        if (flag & (KSW_EZ_SPLICE_FOR | KSW_EZ_SPLICE_REV)) != 0 {
            let semi_cost: i32 = if (flag & KSW_EZ_SPLICE_FLANK) != 0 { -noncan_i / 2 } else { 0 };
            {
                let lane = _mm_set1_epi8((-noncan_i) as i8);
                for i in 0..tlen_ {
                    donor[i] = lane;
                    acceptor[i] = lane;
                }
            }
            let dp = donor.as_mut_ptr() as *mut i8;
            let ap = acceptor.as_mut_ptr() as *mut i8;
            if (flag & KSW_EZ_REV_CIGAR) == 0 {
                for t in 0..(tlen as usize).saturating_sub(4) {
                    let mut can_type = 0i32; // 0=none, 1=GT/AG only, 2=GTr/yAG
                    if (flag & KSW_EZ_SPLICE_FOR) != 0 && target[t + 1] == 2 && target[t + 2] == 3 {
                        can_type = 1; // GTr...
                    }
                    if (flag & KSW_EZ_SPLICE_REV) != 0 && target[t + 1] == 1 && target[t + 2] == 3 {
                        can_type = 1; // CTr...
                    }
                    if can_type != 0 && (target[t + 3] == 0 || target[t + 3] == 2) {
                        can_type = 2;
                    }
                    if can_type != 0 {
                        *dp.add(t) = if can_type == 2 { 0 } else { semi_cost as i8 };
                    }
                }
                if !junc.is_empty() {
                    for t in 0..(tlen as usize).saturating_sub(1) {
                        if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t + 1] & 1) != 0)
                            || ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t + 1] & 8) != 0)
                        {
                            *dp.add(t) = dp.add(t).read().wrapping_add(junc_bonus);
                        }
                    }
                }
                for t in 2..tlen as usize {
                    let mut can_type = 0i32;
                    if (flag & KSW_EZ_SPLICE_FOR) != 0 && target[t - 1] == 0 && target[t] == 2 {
                        can_type = 1; // ...yAG
                    }
                    if (flag & KSW_EZ_SPLICE_REV) != 0 && target[t - 1] == 0 && target[t] == 1 {
                        can_type = 1; // ...yAC
                    }
                    if can_type != 0 && (target[t - 2] == 1 || target[t - 2] == 3) {
                        can_type = 2;
                    }
                    if can_type != 0 {
                        *ap.add(t) = if can_type == 2 { 0 } else { semi_cost as i8 };
                    }
                }
                if !junc.is_empty() {
                    for t in 0..tlen as usize {
                        if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t] & 2) != 0)
                            || ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t] & 4) != 0)
                        {
                            *ap.add(t) = ap.add(t).read().wrapping_add(junc_bonus);
                        }
                    }
                }
            } else {
                for t in 0..(tlen as usize).saturating_sub(4) {
                    let mut can_type = 0i32;
                    if (flag & KSW_EZ_SPLICE_FOR) != 0 && target[t + 1] == 2 && target[t + 2] == 0 {
                        can_type = 1; // GAy...
                    }
                    if (flag & KSW_EZ_SPLICE_REV) != 0 && target[t + 1] == 1 && target[t + 2] == 0 {
                        can_type = 1; // CAy...
                    }
                    if can_type != 0 && (target[t + 3] == 1 || target[t + 3] == 3) {
                        can_type = 2;
                    }
                    if can_type != 0 {
                        *dp.add(t) = if can_type == 2 { 0 } else { semi_cost as i8 };
                    }
                }
                if !junc.is_empty() {
                    for t in 0..(tlen as usize).saturating_sub(1) {
                        if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t + 1] & 2) != 0)
                            || ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t + 1] & 4) != 0)
                        {
                            *dp.add(t) = dp.add(t).read().wrapping_add(junc_bonus);
                        }
                    }
                }
                for t in 2..tlen as usize {
                    let mut can_type = 0i32;
                    if (flag & KSW_EZ_SPLICE_FOR) != 0 && target[t - 1] == 3 && target[t] == 2 {
                        can_type = 1; // ...rTG
                    }
                    if (flag & KSW_EZ_SPLICE_REV) != 0 && target[t - 1] == 3 && target[t] == 1 {
                        can_type = 1; // ...rTC
                    }
                    if can_type != 0 && (target[t - 2] == 0 || target[t - 2] == 2) {
                        can_type = 2;
                    }
                    if can_type != 0 {
                        *ap.add(t) = if can_type == 2 { 0 } else { semi_cost as i8 };
                    }
                }
                if !junc.is_empty() {
                    for t in 0..tlen as usize {
                        if ((flag & KSW_EZ_SPLICE_FOR) != 0 && (junc[t] & 1) != 0)
                            || ((flag & KSW_EZ_SPLICE_REV) != 0 && (junc[t] & 8) != 0)
                        {
                            *ap.add(t) = ap.add(t).read().wrapping_add(junc_bonus);
                        }
                    }
                }
            }
        }

        let mut last_st = -1i32;
        let mut last_en = -1i32;
        let mut h0 = 0i32;
        let mut last_h0_t = 0i32;

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
            let st0 = st;
            let en0 = en;
            st = st / 16 * 16;
            en = (en + 16) / 16 * 16 - 1;
            // set boundary conditions
            let x1: i8;
            let x21: i8;
            let v1: i8;
            if st > 0 {
                if st - 1 >= last_st && st - 1 <= last_en {
                    x1 = *xp.add((st - 1) as usize); // (r-1,s-1) calculated in the last round
                    x21 = *x2p.add((st - 1) as usize);
                    v1 = *vp.add((st - 1) as usize);
                } else {
                    x1 = (-qe) as i8;
                    x21 = (-q2_i) as i8;
                    v1 = (-qe) as i8;
                }
            } else {
                x1 = (-qe) as i8;
                x21 = (-q2_i) as i8;
                v1 = if r == 0 {
                    (-qe) as i8
                } else if r < long_thres {
                    (-e_i) as i8
                } else if r == long_thres {
                    long_diff as i8
                } else {
                    0
                };
            }
            if en >= r {
                *yp.add(r as usize) = (-qe) as i8;
                *up.add(r as usize) = if r == 0 {
                    (-qe) as i8
                } else if r < long_thres {
                    (-e_i) as i8
                } else if r == long_thres {
                    long_diff as i8
                } else {
                    0
                };
            }
            // loop fission: set scores first
            let qrr_base = qlen_u as isize - 1 - r as isize;
            if (flag & KSW_EZ_GENERIC_SC) == 0 {
                let mut t = st0;
                while t <= en0 {
                    let tu = t as usize;
                    let sq = _mm_loadu_si128(sfp.add(tu) as *const __m128i);
                    let stq =
                        _mm_loadu_si128(qrp.offset(qrr_base + t as isize) as *const __m128i);
                    let mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(stq, m1_));
                    let mut tmp = _mm_cmpeq_epi8(sq, stq);
                    tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
                    tmp = _mm_blendv_epi8(tmp, sc_amb_, mask);
                    _mm_storeu_si128(sp.add(tu) as *mut __m128i, tmp);
                    t += 16;
                }
            } else {
                let mut t = st0;
                while t <= en0 {
                    let tu = t as usize;
                    *sp.add(tu) = mat[*sfp.add(tu) as usize * mu
                        + *qrp.offset(qrr_base + t as isize) as usize]
                        as u8;
                    t += 1;
                }
            }
            // core loop
            let mut x1_ = _mm_cvtsi32_si128((x1 as u8) as i32); // C: (uint8_t)x1
            let mut x21_ = _mm_cvtsi32_si128((x21 as u8) as i32); // C: (uint8_t)x21
            let mut v1_ = _mm_cvtsi32_si128((v1 as u8) as i32); // C: (uint8_t)v1
            let st_ = (st / 16) as usize;
            let en_ = (en / 16) as usize;
            assert!(en_ - st_ + 1 <= n_col_);
            if !with_cigar {
                // score only
                let mut t = st_;
                while t <= en_ {
                    let mut z = _mm_load_si128(s.as_ptr().add(t));
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
                    let mut x2t1 = _mm_load_si128(x2.as_ptr().add(t));
                    tmp = _mm_srli_si128::<15>(x2t1);
                    x2t1 = _mm_or_si128(_mm_slli_si128::<1>(x2t1), x21_);
                    x21_ = tmp;
                    let mut a2 = _mm_add_epi8(x2t1, vt1);
                    let a2a = _mm_add_epi8(a2, _mm_load_si128(acceptor.as_ptr().add(t)));
                    z = _mm_max_epi8(z, a);
                    z = _mm_max_epi8(z, b);
                    z = _mm_max_epi8(z, a2a);
                    // save u[] and v[]; update a, b and a2
                    _mm_store_si128(u.as_mut_ptr().add(t), _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
                    _mm_store_si128(v.as_mut_ptr().add(t), _mm_sub_epi8(z, ut)); // v[r][t..t+15] <- z - u[r-1][t..t+15]
                    tmp = _mm_sub_epi8(z, q_);
                    a = _mm_sub_epi8(a, tmp);
                    b = _mm_sub_epi8(b, tmp);
                    a2 = _mm_sub_epi8(a2, _mm_sub_epi8(z, q2_));
                    _mm_store_si128(x.as_mut_ptr().add(t), _mm_sub_epi8(_mm_max_epi8(a, zero_), qe_));
                    _mm_store_si128(y.as_mut_ptr().add(t), _mm_sub_epi8(_mm_max_epi8(b, zero_), qe_));
                    tmp = _mm_load_si128(donor.as_ptr().add(t));
                    _mm_store_si128(x2.as_mut_ptr().add(t), _mm_sub_epi8(_mm_max_epi8(a2, tmp), q2_));
                    t += 1;
                }
            } else if (flag & KSW_EZ_RIGHT) == 0 {
                // gap left-alignment
                let pr_off = r as usize * n_col_ * 16 - st as usize;
                off[r as usize] = st;
                off_end[r as usize] = en;
                let mut t = st_;
                while t <= en_ {
                    let mut z = _mm_load_si128(s.as_ptr().add(t));
                    let mut xt1 = _mm_load_si128(x.as_ptr().add(t));
                    let mut tmp = _mm_srli_si128::<15>(xt1);
                    xt1 = _mm_or_si128(_mm_slli_si128::<1>(xt1), x1_);
                    x1_ = tmp;
                    let mut vt1 = _mm_load_si128(v.as_ptr().add(t));
                    tmp = _mm_srli_si128::<15>(vt1);
                    vt1 = _mm_or_si128(_mm_slli_si128::<1>(vt1), v1_);
                    v1_ = tmp;
                    let mut a = _mm_add_epi8(xt1, vt1);
                    let ut = _mm_load_si128(u.as_ptr().add(t));
                    let mut b = _mm_add_epi8(_mm_load_si128(y.as_ptr().add(t)), ut);
                    let mut x2t1 = _mm_load_si128(x2.as_ptr().add(t));
                    tmp = _mm_srli_si128::<15>(x2t1);
                    x2t1 = _mm_or_si128(_mm_slli_si128::<1>(x2t1), x21_);
                    x21_ = tmp;
                    let mut a2 = _mm_add_epi8(x2t1, vt1);
                    let a2a = _mm_add_epi8(a2, _mm_load_si128(acceptor.as_ptr().add(t)));

                    let mut d = _mm_and_si128(_mm_cmpgt_epi8(a, z), _mm_set1_epi8(1)); // d = a  > z? 1 : 0
                    z = _mm_max_epi8(z, a);
                    d = _mm_blendv_epi8(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b, z)); // d = b  > z? 2 : d
                    z = _mm_max_epi8(z, b);
                    d = _mm_blendv_epi8(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2a, z)); // d = a2 > z? 3 : d
                    z = _mm_max_epi8(z, a2a);
                    // save u[] and v[]; update a, b and a2
                    _mm_store_si128(u.as_mut_ptr().add(t), _mm_sub_epi8(z, vt1));
                    _mm_store_si128(v.as_mut_ptr().add(t), _mm_sub_epi8(z, ut));
                    tmp = _mm_sub_epi8(z, q_);
                    a = _mm_sub_epi8(a, tmp);
                    b = _mm_sub_epi8(b, tmp);
                    a2 = _mm_sub_epi8(a2, _mm_sub_epi8(z, q2_));

                    tmp = _mm_cmpgt_epi8(a, zero_);
                    _mm_store_si128(
                        x.as_mut_ptr().add(t),
                        _mm_sub_epi8(_mm_and_si128(tmp, a), qe_),
                    );
                    d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                    tmp = _mm_cmpgt_epi8(b, zero_);
                    _mm_store_si128(
                        y.as_mut_ptr().add(t),
                        _mm_sub_epi8(_mm_and_si128(tmp, b), qe_),
                    );
                    d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0

                    let tmp2 = _mm_load_si128(donor.as_ptr().add(t));
                    tmp = _mm_cmpgt_epi8(a2, tmp2);
                    let mm = _mm_max_epi8(a2, tmp2);
                    _mm_store_si128(x2.as_mut_ptr().add(t), _mm_sub_epi8(mm, q2_));
                    d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x20)));
                    _mm_storeu_si128(p.as_mut_ptr().add(pr_off + t * 16) as *mut __m128i, d);
                    t += 1;
                }
            } else {
                // gap right-alignment
                let pr_off = r as usize * n_col_ * 16 - st as usize;
                off[r as usize] = st;
                off_end[r as usize] = en;
                let mut t = st_;
                while t <= en_ {
                    let mut z = _mm_load_si128(s.as_ptr().add(t));
                    let mut xt1 = _mm_load_si128(x.as_ptr().add(t));
                    let mut tmp = _mm_srli_si128::<15>(xt1);
                    xt1 = _mm_or_si128(_mm_slli_si128::<1>(xt1), x1_);
                    x1_ = tmp;
                    let mut vt1 = _mm_load_si128(v.as_ptr().add(t));
                    tmp = _mm_srli_si128::<15>(vt1);
                    vt1 = _mm_or_si128(_mm_slli_si128::<1>(vt1), v1_);
                    v1_ = tmp;
                    let mut a = _mm_add_epi8(xt1, vt1);
                    let ut = _mm_load_si128(u.as_ptr().add(t));
                    let mut b = _mm_add_epi8(_mm_load_si128(y.as_ptr().add(t)), ut);
                    let mut x2t1 = _mm_load_si128(x2.as_ptr().add(t));
                    tmp = _mm_srli_si128::<15>(x2t1);
                    x2t1 = _mm_or_si128(_mm_slli_si128::<1>(x2t1), x21_);
                    x21_ = tmp;
                    let mut a2 = _mm_add_epi8(x2t1, vt1);
                    let a2a = _mm_add_epi8(a2, _mm_load_si128(acceptor.as_ptr().add(t)));

                    let mut d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), _mm_set1_epi8(1)); // d = z > a?  0 : 1
                    z = _mm_max_epi8(z, a);
                    d = _mm_blendv_epi8(_mm_set1_epi8(2), d, _mm_cmpgt_epi8(z, b)); // d = z > b?  d : 2
                    z = _mm_max_epi8(z, b);
                    d = _mm_blendv_epi8(_mm_set1_epi8(3), d, _mm_cmpgt_epi8(z, a2a)); // d = z > a2? d : 3
                    z = _mm_max_epi8(z, a2a);
                    // save u[] and v[]; update a, b and a2
                    _mm_store_si128(u.as_mut_ptr().add(t), _mm_sub_epi8(z, vt1));
                    _mm_store_si128(v.as_mut_ptr().add(t), _mm_sub_epi8(z, ut));
                    tmp = _mm_sub_epi8(z, q_);
                    a = _mm_sub_epi8(a, tmp);
                    b = _mm_sub_epi8(b, tmp);
                    a2 = _mm_sub_epi8(a2, _mm_sub_epi8(z, q2_));

                    tmp = _mm_cmpgt_epi8(zero_, a);
                    _mm_store_si128(
                        x.as_mut_ptr().add(t),
                        _mm_sub_epi8(_mm_andnot_si128(tmp, a), qe_),
                    );
                    d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                    tmp = _mm_cmpgt_epi8(zero_, b);
                    _mm_store_si128(
                        y.as_mut_ptr().add(t),
                        _mm_sub_epi8(_mm_andnot_si128(tmp, b), qe_),
                    );
                    d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0

                    let tmp2 = _mm_load_si128(donor.as_ptr().add(t));
                    tmp = _mm_cmpgt_epi8(tmp2, a2);
                    let mm = _mm_max_epi8(tmp2, a2);
                    _mm_store_si128(x2.as_mut_ptr().add(t), _mm_sub_epi8(mm, q2_));
                    d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
                    _mm_storeu_si128(p.as_mut_ptr().add(pr_off + t * 16) as *mut __m128i, d);
                    t += 1;
                }
            }

            if !approx_max {
                // find the exact max with a 32-bit score array
                let max_h: i32;
                let max_t: i32;
                if r > 0 {
                    let en0u = en0 as usize;
                    let mut mh = if en0 > 0 {
                        h[en0u - 1] + *up.add(en0u) as i32
                    } else {
                        h[en0u] + *vp.add(en0u) as i32
                    };
                    h[en0u] = mh; // special casing the last element
                    let mut mt = en0;
                    let mut max_h_ = _mm_set1_epi32(mh);
                    let mut max_t_ = _mm_set1_epi32(mt);
                    let en1 = st0 + (en0 - st0) / 4 * 4;
                    let mut t = st0;
                    while t < en1 {
                        let tu = t as usize;
                        let mut h1_ = _mm_loadu_si128(h.as_ptr().add(tu) as *const __m128i);
                        let tv = _mm_setr_epi32(
                            *vp.add(tu) as i32,
                            *vp.add(tu + 1) as i32,
                            *vp.add(tu + 2) as i32,
                            *vp.add(tu + 3) as i32,
                        );
                        h1_ = _mm_add_epi32(h1_, tv);
                        _mm_storeu_si128(h.as_mut_ptr().add(tu) as *mut __m128i, h1_);
                        let tt = _mm_set1_epi32(t);
                        let tmp = _mm_cmpgt_epi32(h1_, max_h_);
                        max_h_ = _mm_blendv_epi8(max_h_, h1_, tmp);
                        max_t_ = _mm_blendv_epi8(max_t_, tt, tmp);
                        t += 4;
                    }
                    let mut hh = [0i32; 4];
                    let mut tt = [0i32; 4];
                    _mm_storeu_si128(hh.as_mut_ptr() as *mut __m128i, max_h_);
                    _mm_storeu_si128(tt.as_mut_ptr() as *mut __m128i, max_t_);
                    for i in 0..4 {
                        if mh < hh[i] {
                            mh = hh[i];
                            mt = tt[i] + i as i32;
                        }
                    }
                    while t < en0 {
                        let tu = t as usize;
                        h[tu] += *vp.add(tu) as i32;
                        if h[tu] > mh {
                            mh = h[tu];
                            mt = t;
                        }
                        t += 1;
                    }
                    max_h = mh;
                    max_t = mt;
                } else {
                    h[0] = *vp as i32 - qe; // special casing r==0
                    max_h = h[0];
                    max_t = 0;
                }
                // update ez
                if en0 == tlen - 1 && h[en0 as usize] > ez.mte {
                    ez.mte = h[en0 as usize];
                    ez.mte_q = r - en;
                }
                if r - st0 == qlen - 1 && h[st0 as usize] > ez.mqe {
                    ez.mqe = h[st0 as usize];
                    ez.mqe_t = st0;
                }
                if apply_zdrop(ez, true, max_h, r, max_t, zdrop, 0) {
                    break;
                }
                if r == qlen + tlen - 2 && en0 == tlen - 1 {
                    ez.score = h[tlen_u - 1];
                }
            } else {
                // find approximate max; Z-drop might be inaccurate, too.
                if r > 0 {
                    if last_h0_t >= st0
                        && last_h0_t <= en0
                        && last_h0_t + 1 >= st0
                        && last_h0_t + 1 <= en0
                    {
                        let d0 = *vp.add(last_h0_t as usize) as i32;
                        let d1 = *up.add(last_h0_t as usize + 1) as i32;
                        if d0 > d1 {
                            h0 += d0;
                        } else {
                            h0 += d1;
                            last_h0_t += 1;
                        }
                    } else if last_h0_t >= st0 && last_h0_t <= en0 {
                        h0 += *vp.add(last_h0_t as usize) as i32;
                    } else {
                        last_h0_t += 1;
                        h0 += *up.add(last_h0_t as usize) as i32;
                    }
                } else {
                    h0 = *vp as i32 - qe;
                    last_h0_t = 0;
                }
                // note: unlike extz2, the C code calls ksw_apply_zdrop() outside
                // the "if (r > 0)" block here
                if (flag & KSW_EZ_APPROX_DROP) != 0
                    && apply_zdrop(ez, true, h0, r, last_h0_t, zdrop, 0)
                {
                    break;
                }
                if r == qlen + tlen - 2 && en0 == tlen - 1 {
                    ez.score = h0;
                }
            }
            last_st = st;
            last_en = en;
        }

        if with_cigar {
            let rev_cigar = (flag & KSW_EZ_REV_CIGAR) != 0;
            let p_bytes = &p[..n_rows * n_col_ * 16];
            let nc = n_col_ as i32 * 16;
            if !ez.zdropped && (flag & KSW_EZ_EXTZ_ONLY) == 0 {
                backtrack(
                    true,
                    rev_cigar,
                    long_thres,
                    p_bytes,
                    &off,
                    Some(&off_end),
                    nc,
                    tlen - 1,
                    qlen - 1,
                    &mut ez.cigar,
                );
            } else if ez.max_t >= 0 && ez.max_q >= 0 {
                backtrack(
                    true,
                    rev_cigar,
                    long_thres,
                    p_bytes,
                    &off,
                    Some(&off_end),
                    nc,
                    ez.max_t,
                    ez.max_q,
                    &mut ez.cigar,
                );
            }
        }
    }
}
