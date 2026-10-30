//! Port of `ksw2_extd.c`: global/extension alignment with a two-piece affine
//! gap cost, Green's standard formulation.

use crate::common::{apply_zdrop, backtrack, reset_extz};
use crate::{
    KSW_EZ_EXTZ_ONLY, KSW_EZ_REV_CIGAR, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF, KswExtz,
};

pub fn ksw_extd(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    gapo2: i8,
    gape2: i8,
    w: i32,
    zdrop: i32,
    flag: i32,
    ez: &mut KswExtz,
) {
    reset_extz(ez);
    if m <= 0 || qlen <= 0 || tlen <= 0 {
        return;
    }

    let mu = m as usize;
    let gape_i = gape as i32;
    let gape2_i = gape2 as i32;
    let gapoe = gapo as i32 + gape_i;
    let gapoe2 = gapo2 as i32 + gape2_i;
    let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };
    let n_col = if qlen < 2 * w + 1 { qlen } else { 2 * w + 1 };
    let with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0;
    let qlen_u = qlen as usize;
    let tlen_u = tlen as usize;
    let n_col_u = n_col as usize;

    let mut qp = vec![0i8; qlen_u * mu];
    let mut eh = vec![[0i32; 4]; qlen_u + 1]; // { h, e, e2, pad }
    let mut z: Vec<u8> = Vec::new();
    let mut off: Vec<i32> = Vec::new();
    if with_cigar {
        z = vec![0u8; n_col_u * tlen_u];
        off = vec![0i32; tlen_u];
    }

    // generate the query profile
    let mut i = 0usize;
    for k in 0..mu {
        let p = &mat[k * mu..k * mu + mu];
        for j in 0..qlen_u {
            qp[i] = p[query[j] as usize];
            i += 1;
        }
    }

    // fill the first row
    eh[0] = [0, -gapoe - gapoe, -gapoe2 - gapoe2, 0];
    let mut j = 1i32;
    while j <= qlen && j <= w {
        let t1 = -(gapo as i32 + gape_i * j);
        let t2 = -(gapo2 as i32 + gape2_i * j);
        let a = -(gapoe + gape_i * j);
        let b = -(gapoe2 + gape2_i * j);
        let tmp = if a > b { a } else { b };
        eh[j as usize] = [if t1 > t2 { t1 } else { t2 }, tmp - gapoe, tmp - gapoe2, 0];
        j += 1;
    }
    while j <= qlen {
        eh[j as usize] = [KSW_NEG_INF, KSW_NEG_INF, KSW_NEG_INF, 0];
        j += 1;
    }

    let qpp = qp.as_ptr();
    let ehp = eh.as_mut_ptr();
    let zp = z.as_mut_ptr();
    let offp = off.as_mut_ptr();
    let tgt = target.as_ptr();
    let right = (flag & KSW_EZ_RIGHT) != 0;
    let mut max_j = 0i32;

    for i in 0..tlen {
        let mut max = KSW_NEG_INF;
        let st = if i > w { i - w } else { 0 };
        let en = if i + w < qlen - 1 { i + w } else { qlen - 1 };
        let a = -(gapoe + gape_i * i);
        let b = -(gapoe2 + gape2_i * i);
        let tmp = if a > b { a } else { b };
        let mut h1 = if st > 0 { KSW_NEG_INF } else { tmp };
        let mut f = if st > 0 { KSW_NEG_INF } else { tmp - gapoe };
        let mut f2 = if st > 0 { KSW_NEG_INF } else { tmp - gapoe2 };
        unsafe {
            let q = qpp.add(*tgt.add(i as usize) as usize * qlen_u);
            if !with_cigar {
                let (nh1, nmax, nmax_j) = score_row(
                    ehp, q, st, en, gapoe, gape_i, gapoe2, gape2_i, h1, f, f2, max, max_j,
                );
                h1 = nh1;
                max = nmax;
                max_j = nmax_j;
            } else if !right {
                let mut zi = zp.add(i as usize * n_col_u);
                *offp.add(i as usize) = st;
                let mut p = ehp.add(st as usize);
                let mut qj = q.add(st as usize);
                for j in st..=en {
                    let mut cell = *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    let mut e2 = cell[2];
                    cell[0] = h1;
                    let mut h = h_prev + *qj as i32;
                    let mut d: u8 = if h >= e { 0 } else { 1 };
                    h = if e > h { e } else { h };
                    d = if h >= f { d } else { 2 };
                    h = if f > h { f } else { h };
                    d = if h >= e2 { d } else { 3 };
                    h = if e2 > h { e2 } else { h };
                    d = if h >= f2 { d } else { 4 };
                    h = if f2 > h { f2 } else { h };
                    h1 = h;
                    let gt = max > h;
                    max_j = if gt { max_j } else { j };
                    max = if gt { max } else { h };
                    h -= gapoe;
                    e -= gape_i;
                    d |= if e > h { 1 << 3 } else { 0 };
                    e = if h > e { h } else { e };
                    cell[1] = e;
                    f -= gape_i;
                    d |= if f > h { 1 << 4 } else { 0 };
                    f = if h > f { h } else { f };
                    let h2 = h1 - gapoe2;
                    e2 -= gape2_i;
                    d |= if e2 > h2 { 1 << 5 } else { 0 };
                    e2 = if h2 > e2 { h2 } else { e2 };
                    cell[2] = e2;
                    f2 -= gape2_i;
                    d |= if f2 > h2 { 1 << 6 } else { 0 };
                    f2 = if h2 > f2 { h2 } else { f2 };
                    *zi = d;
                    zi = zi.add(1);
                    *p = cell;
                    p = p.add(1);
                    qj = qj.add(1);
                }
            } else {
                let mut zi = zp.add(i as usize * n_col_u);
                *offp.add(i as usize) = st;
                let mut p = ehp.add(st as usize);
                let mut qj = q.add(st as usize);
                for j in st..=en {
                    let mut cell = *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    let mut e2 = cell[2];
                    cell[0] = h1;
                    let mut h = h_prev + *qj as i32;
                    let mut d: u8 = if h > e { 0 } else { 1 };
                    h = if e > h { e } else { h };
                    d = if h > f { d } else { 2 };
                    h = if f > h { f } else { h };
                    d = if h > e2 { d } else { 3 };
                    h = if e2 > h { e2 } else { h };
                    d = if h > f2 { d } else { 4 };
                    h = if f2 > h { f2 } else { h };
                    h1 = h;
                    // note: unlike ksw2_extz.c, all three branches of
                    // ksw2_extd.c use `max > h` here (no right-align variant)
                    let gt = max > h;
                    max_j = if gt { max_j } else { j };
                    max = if gt { max } else { h };
                    h -= gapoe;
                    e -= gape_i;
                    d |= if e >= h { 1 << 3 } else { 0 };
                    e = if h > e { h } else { e };
                    cell[1] = e;
                    f -= gape_i;
                    d |= if f >= h { 1 << 4 } else { 0 };
                    f = if h > f { h } else { f };
                    let h2 = h1 - gapoe2;
                    e2 -= gape2_i;
                    d |= if e2 >= h2 { 1 << 5 } else { 0 };
                    e2 = if h2 > e2 { h2 } else { e2 };
                    cell[2] = e2;
                    f2 -= gape2_i;
                    d |= if f2 >= h2 { 1 << 6 } else { 0 };
                    f2 = if h2 > f2 { h2 } else { f2 };
                    *zi = d;
                    zi = zi.add(1);
                    *p = cell;
                    p = p.add(1);
                    qj = qj.add(1);
                }
            }
            let jj = ((en + 1) as usize).min(qlen_u);
            (*ehp.add(jj))[0] = h1;
            (*ehp.add(jj))[1] = KSW_NEG_INF;
        }
        // update ez
        if en == qlen - 1 && eh[qlen_u][0] > ez.mqe {
            ez.mqe = eh[qlen_u][0];
            ez.mqe_t = i;
        }
        if i == tlen - 1 {
            ez.mte = max;
            ez.mte_q = max_j;
        }
        if apply_zdrop(ez, false, max, i, max_j, zdrop, gape2) {
            break;
        }
        if i == tlen - 1 && en == qlen - 1 {
            ez.score = eh[qlen_u][0];
        }
    }

    if with_cigar {
        let rev_cigar = (flag & KSW_EZ_REV_CIGAR) != 0;
        if !ez.zdropped && (flag & KSW_EZ_EXTZ_ONLY) == 0 {
            backtrack(
                false,
                rev_cigar,
                0,
                &z,
                &off,
                None,
                n_col,
                tlen - 1,
                qlen - 1,
                &mut ez.cigar,
            );
        } else if ez.max_t >= 0 && ez.max_q >= 0 {
            backtrack(
                false,
                rev_cigar,
                0,
                &z,
                &off,
                None,
                n_col,
                ez.max_t,
                ez.max_q,
                &mut ez.cigar,
            );
        }
    }
}

/// Score-only inner loop of [`ksw_extd`], kept in its own function so that
/// LLVM does not run out of registers in the outer loop.
#[inline(never)]
unsafe fn score_row(
    ehp: *mut [i32; 4],
    q: *const i8,
    st: i32,
    en: i32,
    gapoe: i32,
    gape: i32,
    gapoe2: i32,
    gape2: i32,
    h1_in: i32,
    f_in: i32,
    f2_in: i32,
    max_in: i32,
    max_j_in: i32,
) -> (i32, i32, i32) {
    let mut h1 = h1_in;
    let mut f = f_in;
    let mut f2 = f2_in;
    let mut max = max_in;
    let mut max_j = max_j_in;
    unsafe {
        let mut p = ehp.add(st as usize);
        let mut qj = q.add(st as usize);
        for j in st..=en {
            let mut cell = *p;
            let h_prev = cell[0];
            let mut e = cell[1];
            let mut e2 = cell[2];
            cell[0] = h1;
            let mut h = h_prev + *qj as i32;
            h = h.max(e);
            h = h.max(f);
            h = h.max(e2);
            h = h.max(f2);
            h1 = h;
            let gt = max > h;
            max_j = if gt { max_j } else { j };
            max = if gt { max } else { h };
            h -= gapoe;
            e -= gape;
            e = e.max(h);
            cell[1] = e;
            f -= gape;
            f = f.max(h);
            let h2 = h1 - gapoe2;
            e2 -= gape2;
            e2 = e2.max(h2);
            cell[2] = e2;
            f2 -= gape2;
            f2 = f2.max(h2);
            *p = cell;
            p = p.add(1);
            qj = qj.add(1);
        }
    }
    (h1, max, max_j)
}
