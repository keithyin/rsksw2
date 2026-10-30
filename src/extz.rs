//! Port of `ksw2_extz.c`: global/extension alignment with an affine gap cost,
//! Green's standard formulation.

use crate::common::{apply_zdrop, backtrack, reset_extz};
use crate::{
    KSW_EZ_EXTZ_ONLY, KSW_EZ_REV_CIGAR, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF, KswExtz,
};

pub fn ksw_extz(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
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
    let gapoe = gapo as i32 + gape_i;
    let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };
    let n_col = if qlen < 2 * w + 1 { qlen } else { 2 * w + 1 }; // maximum #columns of the backtrack matrix
    let with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0;
    let qlen_u = qlen as usize;
    let tlen_u = tlen as usize;
    let n_col_u = n_col as usize;

    // allocate memory
    let mut qp = vec![0i8; qlen_u * mu];
    let mut eh = vec![[0i32; 2]; qlen_u + 1]; // { h, e }
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
    eh[0] = [0, -gapoe - gapoe];
    let mut j = 1i32;
    while j <= qlen && j <= w {
        eh[j as usize] = [-(gapoe + gape_i * (j - 1)), -(gapoe + gapoe + gape_i * j)];
        j += 1;
    }
    while j <= qlen {
        // everything is -inf outside the band
        eh[j as usize] = [KSW_NEG_INF, KSW_NEG_INF];
        j += 1;
    }

    let qpp = qp.as_ptr();
    let ehp = eh.as_mut_ptr();
    let zp = z.as_mut_ptr();
    let offp = off.as_mut_ptr();
    let tgt = target.as_ptr();
    let right = (flag & KSW_EZ_RIGHT) != 0;
    let mut max_j = 0i32;

    // DP loop (target sequence is in the outer loop)
    for i in 0..tlen {
        let mut max = KSW_NEG_INF;
        let st = if i > w { i - w } else { 0 };
        let en = if i + w < qlen - 1 { i + w } else { qlen - 1 };
        // a gap of length l costs -(gapo + l*gape)
        let neg_h = -(gapoe + gape_i * i);
        let mut h1 = if st > 0 { KSW_NEG_INF } else { neg_h };
        let mut f = if st > 0 { KSW_NEG_INF } else { neg_h - gapoe };
        unsafe {
            let q = qpp.add(*tgt.add(i as usize) as usize * qlen_u);
            let mut p = ehp.add(st as usize);
            if !with_cigar {
                for j in st..=en {
                    let cell = &mut *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    cell[0] = h1;
                    let mut h = h_prev + *q.add(j as usize) as i32;
                    h = h.max(e);
                    h = h.max(f);
                    h1 = h;
                    let gt = max > h;
                    max_j = if gt { max_j } else { j };
                    max = if gt { max } else { h };
                    h -= gapoe;
                    e -= gape_i;
                    e = e.max(h);
                    cell[1] = e;
                    f -= gape_i;
                    f = f.max(h);
                    p = p.add(1);
                }
            } else if !right {
                let mut zi = zp.add(i as usize * n_col_u);
                *offp.add(i as usize) = st;
                for j in st..=en {
                    let cell = &mut *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    cell[0] = h1;
                    let mut h = h_prev + *q.add(j as usize) as i32;
                    let mut d: u8 = if h >= e { 0 } else { 1 };
                    h = h.max(e);
                    d = if h >= f { d } else { 2 };
                    h = h.max(f);
                    h1 = h;
                    let gt = max > h;
                    max_j = if gt { max_j } else { j };
                    max = if gt { max } else { h };
                    h -= gapoe;
                    e -= gape_i;
                    d |= if e > h { 0x08 } else { 0 };
                    e = e.max(h);
                    cell[1] = e;
                    f -= gape_i;
                    d |= if f > h { 0x10 } else { 0 };
                    f = f.max(h);
                    *zi = d; // z[i,j] keeps h for the current cell and e/f for the next cell
                    zi = zi.add(1);
                    p = p.add(1);
                }
            } else {
                let mut zi = zp.add(i as usize * n_col_u);
                *offp.add(i as usize) = st;
                for j in st..=en {
                    let cell = &mut *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    cell[0] = h1;
                    let mut h = h_prev + *q.add(j as usize) as i32;
                    let mut d: u8 = if h > e { 0 } else { 1 };
                    h = h.max(e);
                    d = if h > f { d } else { 2 };
                    h = h.max(f);
                    h1 = h;
                    let ge = max >= h;
                    max_j = if ge { max_j } else { j };
                    max = if ge { max } else { h };
                    h -= gapoe;
                    e -= gape_i;
                    d |= if e >= h { 0x08 } else { 0 };
                    e = e.max(h);
                    cell[1] = e;
                    f -= gape_i;
                    d |= if f >= h { 0x10 } else { 0 };
                    f = f.max(h);
                    *zi = d;
                    zi = zi.add(1);
                    p = p.add(1);
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
        if apply_zdrop(ez, false, max, i, max_j, zdrop, gape) {
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
