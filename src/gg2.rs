//! Port of `ksw2_gg2.c`: global alignment, Suzuki's diagonal formulation.

use crate::common::backtrack;

pub fn ksw_gg2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    w: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    if m <= 0 || qlen <= 0 || tlen <= 0 {
        return 0;
    }
    let mu = m as usize;
    let qe = q as i32 + e as i32;
    let qe2 = qe + qe;
    let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };
    let n_col = if w + 1 < tlen { w + 1 } else { tlen };
    let qlen_u = qlen as usize;
    let tlen_u = tlen as usize;
    let n_col_u = n_col as usize;

    let mut u = vec![0i8; tlen_u + 1];
    let mut v = vec![0i8; tlen_u + 1];
    let mut x = vec![0i8; tlen_u + 1];
    let mut y = vec![0i8; tlen_u + 1];
    let mut s = vec![0i8; tlen_u];
    let mut qr = vec![0u8; qlen_u];
    let mut p: Vec<u8> = Vec::new();
    let mut off: Vec<i32> = Vec::new();
    let with_cigar = cigar.is_some();
    if with_cigar {
        p = vec![0u8; (qlen_u + tlen_u) * n_col_u];
        off = vec![0i32; qlen_u + tlen_u];
    }

    for t in 0..qlen_u {
        qr[t] = query[qlen_u - 1 - t];
    }

    let mut h0 = 0i32;
    let mut last_h0_t = 0i32;
    let up = u.as_mut_ptr();
    let vp = v.as_mut_ptr();
    let xp = x.as_mut_ptr();
    let yp = y.as_mut_ptr();
    let sp = s.as_mut_ptr();
    let qrp = qr.as_ptr();
    let matp = mat.as_ptr();
    let tgt = target.as_ptr();
    let pp = p.as_mut_ptr();
    let offp = off.as_mut_ptr();

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
        let st_u = st as usize;
        let en_u = en as usize;
        let qrr = qlen_u as isize - 1 - r as isize;
        unsafe {
            // set boundary conditions
            let mut x1: i8;
            let mut v1: i8;
            if st != 0 {
                if r > st + st + w - 1 {
                    x1 = 0;
                    v1 = 0;
                } else {
                    x1 = *xp.add(st_u - 1);
                    v1 = *vp.add(st_u - 1); // (r-1, st-1) in the band
                }
            } else {
                x1 = 0;
                v1 = if r != 0 { q } else { 0 };
            }
            if en != r {
                if r < en + en - w - 1 {
                    *yp.add(en_u) = 0;
                    *up.add(en_u) = 0; // (r-1,en) out of the band
                }
            } else {
                *yp.add(r as usize) = 0;
                *up.add(r as usize) = if r != 0 { q } else { 0 };
            }
            // loop fission: set scores first
            for t in st_u..=en_u {
                *sp.add(t) = *matp.add(*tgt.add(t) as usize * mu + *qrp.offset(t as isize + qrr) as usize);
            }
            // core loop
            if with_cigar {
                let pr_base = r as usize * n_col_u;
                *offp.add(r as usize) = st;
                for t in st_u..=en_u {
                    /* a = x(r-1,t-1) + v(r-1,t-1)
                       b = y(r-1,t)   + u(r-1,t)
                       z = max{ S(t,r-t) + 2q + 2r, a, b } */
                    let mut z: i8 = (*sp.add(t) as i32 + qe2) as i8;
                    let mut a: i8 = (x1 as i32 + v1 as i32) as i8;
                    let b: i8 = (*yp.add(t) as i32 + *up.add(t) as i32) as i8;
                    let mut d: u8 = if a > z { 1 } else { 0 };
                    z = z.max(a);
                    d = if b > z { 2 } else { d };
                    z = z.max(b);
                    let u1 = *up.add(t);
                    *up.add(t) = (z as i32 - v1 as i32) as i8;
                    v1 = *vp.add(t);
                    *vp.add(t) = (z as i32 - u1 as i32) as i8;
                    let zz = (z as i32 - q as i32) as i8;
                    a = (a as i32 - zz as i32) as i8;
                    let b = (b as i32 - zz as i32) as i8;
                    x1 = *xp.add(t);
                    d |= if a > 0 { 0x08 } else { 0 };
                    *xp.add(t) = if a > 0 { a } else { 0 };
                    d |= if b > 0 { 0x10 } else { 0 };
                    *yp.add(t) = if b > 0 { b } else { 0 };
                    *pp.add(pr_base + t - st_u) = d;
                }
            } else {
                for t in st_u..=en_u {
                    let mut z: i8 = (*sp.add(t) as i32 + qe2) as i8;
                    let mut a: i8 = (x1 as i32 + v1 as i32) as i8;
                    let b: i8 = (*yp.add(t) as i32 + *up.add(t) as i32) as i8;
                    z = z.max(a);
                    z = z.max(b);
                    let u1 = *up.add(t);
                    *up.add(t) = (z as i32 - v1 as i32) as i8;
                    v1 = *vp.add(t);
                    *vp.add(t) = (z as i32 - u1 as i32) as i8;
                    let zz = (z as i32 - q as i32) as i8;
                    a = (a as i32 - zz as i32) as i8;
                    let b = (b as i32 - zz as i32) as i8;
                    x1 = *xp.add(t);
                    *xp.add(t) = if a > 0 { a } else { 0 };
                    *yp.add(t) = if b > 0 { b } else { 0 };
                }
            }
        }
        if r > 0 {
            if last_h0_t >= st && last_h0_t <= en {
                h0 += v[last_h0_t as usize] as i32 - qe;
            } else {
                last_h0_t += 1;
                h0 += u[last_h0_t as usize] as i32 - qe;
            }
        } else {
            h0 = v[0] as i32 - qe - qe;
            last_h0_t = 0;
        }
    }

    if let Some(ci) = cigar {
        backtrack(true, false, 0, &p, &off, None, n_col, tlen - 1, qlen - 1, ci);
    }
    h0
}
