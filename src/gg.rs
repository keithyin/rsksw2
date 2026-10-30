//! Port of `ksw2_gg.c`: global alignment, Green's standard formulation.

use crate::common::backtrack;

pub fn ksw_gg(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i8,
    mat: &[i8],
    gapo: i8,
    gape: i8,
    w: i32,
    cigar: Option<&mut Vec<u32>>,
) -> i32 {
    if m <= 0 || qlen <= 0 || tlen <= 0 {
        return crate::KSW_NEG_INF;
    }
    let mu = m as usize;
    let gape_i = gape as i32;
    let gapoe = gapo as i32 + gape_i;
    let w = if w < 0 { if tlen > qlen { tlen } else { qlen } } else { w };
    let n_col = if qlen < 2 * w + 1 { qlen } else { 2 * w + 1 };
    let qlen_u = qlen as usize;
    let tlen_u = tlen as usize;
    let n_col_u = n_col as usize;

    let mut qp = vec![0i8; qlen_u * mu];
    let mut eh = vec![[0i32; 2]; qlen_u + 1]; // { h, e }
    let mut z: Vec<u8> = Vec::new();
    let mut off: Vec<i32> = Vec::new();
    let with_cigar = cigar.is_some();
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
        eh[j as usize] = [crate::KSW_NEG_INF, crate::KSW_NEG_INF];
        j += 1;
    }

    let qpp = qp.as_ptr();
    let ehp = eh.as_mut_ptr();
    let zp = z.as_mut_ptr();
    let offp = off.as_mut_ptr();
    let tgt = target.as_ptr();

    // DP loop (target sequence is in the outer loop)
    for i in 0..tlen {
        let st = if i > w { i - w } else { 0 };
        let en = if i + w + 1 < qlen { i + w + 1 } else { qlen };
        let neg_h = -(gapoe + gape_i * i);
        let mut h1 = if st > 0 { crate::KSW_NEG_INF } else { neg_h };
        let mut f = if st > 0 { crate::KSW_NEG_INF } else { neg_h - gapoe };
        unsafe {
            let q = qpp.add(*tgt.add(i as usize) as usize * qlen_u);
            let mut p = ehp.add(st as usize);
            if with_cigar {
                let mut zi = zp.add(i as usize * n_col_u);
                *offp.add(i as usize) = st;
                for j in st..en {
                    // At the beginning of the loop: cell = { H(i-1,j-1), E(i,j) },
                    // f = F(i,j) and h1 = H(i,j-1)
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
                for j in st..en {
                    let cell = &mut *p;
                    let h_prev = cell[0];
                    let mut e = cell[1];
                    cell[0] = h1;
                    let mut h = h_prev + *q.add(j as usize) as i32;
                    h = h.max(e);
                    h = h.max(f);
                    h1 = h;
                    h -= gapoe;
                    e -= gape_i;
                    e = e.max(h);
                    cell[1] = e;
                    f -= gape_i;
                    f = f.max(h);
                    p = p.add(1);
                }
            }
            let jj = (en as usize).min(qlen_u);
            (*ehp.add(jj))[0] = h1;
            (*ehp.add(jj))[1] = crate::KSW_NEG_INF;
        }
    }

    // backtrack
    let score = eh[qlen_u][0];
    if let Some(ci) = cigar {
        backtrack(false, false, 0, &z, &off, None, n_col, tlen - 1, qlen - 1, ci);
    }
    score
}
