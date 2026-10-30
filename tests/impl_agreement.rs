//! Cross-checks between the independent implementations that ship in this
//! crate: Green's scalar formulation vs. Suzuki's diagonal (scalar and SSE)
//! formulation. They must agree on the alignment score and the CIGAR.

use rsksw2::*;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn base(&mut self) -> u8 {
        (self.next() % 4) as u8
    }
    fn seq(&mut self, n: i32) -> Vec<u8> {
        (0..n).map(|_| self.base()).collect()
    }
}

fn simple_mat(a: i8, b: i8) -> [i8; 25] {
    let mut mat = [0i8; 25];
    gen_simple_mat(5, &mut mat, a, b);
    mat
}

#[test]
fn global_alignment_agrees() {
    let mut rng = Rng(0x2545F4914F6CDD1D);
    let mat = simple_mat(2, -4);
    for _ in 0..300 {
        let qlen = 1 + (rng.next() % 60) as i32;
        let tlen = 1 + (rng.next() % 60) as i32;
        let q = rng.seq(qlen);
        let t = rng.seq(tlen);
        let mut c1 = Vec::new();
        let mut c2 = Vec::new();
        let s1 = ksw_gg(qlen, &q, tlen, &t, 5, &mat, 4, 2, -1, Some(&mut c1));
        let s2 = ksw_gg2(qlen, &q, tlen, &t, 5, &mat, 4, 2, -1, Some(&mut c2));
        let mut c3 = Vec::new();
        let s3 = ksw_gg2_sse(qlen, &q, tlen, &t, 5, &mat, 4, 2, -1, &mut c3);
        assert_eq!(s1, s2, "gg vs gg2 q={q:?} t={t:?}");
        assert_eq!(s1, s3, "gg vs gg2_sse q={q:?} t={t:?}");
        assert_eq!(c1, c2, "cigar gg vs gg2 q={q:?} t={t:?}");
        assert_eq!(c1, c3, "cigar gg vs gg2_sse q={q:?} t={t:?}");
    }
}

#[test]
fn affine_extension_agrees() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mat = simple_mat(2, -4);
    for _ in 0..300 {
        let qlen = 1 + (rng.next() % 50) as i32;
        let tlen = 1 + (rng.next() % 50) as i32;
        let q = rng.seq(qlen);
        let t = rng.seq(tlen);
        let mut e1 = KswExtz::new();
        let mut e2 = KswExtz::new();
        ksw_extz(qlen, &q, tlen, &t, 5, &mat, 4, 2, -1, -1, 0, &mut e1);
        ksw_extz2_sse(qlen, &q, tlen, &t, 5, &mat, 4, 2, -1, -1, 0, 0, &mut e2);
        assert_eq!(e1.score, e2.score, "extz vs extz2_sse q={q:?} t={t:?}");
        assert_eq!(e1.cigar, e2.cigar, "cigar extz vs extz2_sse q={q:?} t={t:?}");
        assert_eq!(e1.max, e2.max, "max extz vs extz2_sse q={q:?} t={t:?}");
        assert_eq!(
            (e1.max_t, e1.max_q),
            (e2.max_t, e2.max_q),
            "max pos extz vs extz2_sse q={q:?} t={t:?}"
        );
    }
}

#[test]
fn dual_affine_extension_agrees() {
    let mut rng = Rng(0xD1B54A32D192ED03);
    let mat = simple_mat(2, -4);
    for _ in 0..300 {
        let qlen = 1 + (rng.next() % 50) as i32;
        let tlen = 1 + (rng.next() % 50) as i32;
        let q = rng.seq(qlen);
        let t = rng.seq(tlen);
        let mut e1 = KswExtz::new();
        let mut e2 = KswExtz::new();
        ksw_extd(qlen, &q, tlen, &t, 5, &mat, 4, 2, 13, 1, -1, -1, 0, &mut e1);
        ksw_extd2_sse(qlen, &q, tlen, &t, 5, &mat, 4, 2, 13, 1, -1, -1, 0, 0, &mut e2);
        assert_eq!(e1.score, e2.score, "extd vs extd2_sse q={q:?} t={t:?}");
        assert_eq!(e1.cigar, e2.cigar, "cigar extd vs extd2_sse q={q:?} t={t:?}");
    }
}

#[test]
fn scalar_ports_match_reference_example() {
    // the example from the ksw2 README
    let mat = simple_mat(1, -2);
    let enc = |s: &str| -> Vec<u8> {
        s.bytes()
            .map(|c| match c.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect()
    };
    let t = enc("ATAGCTAGCTAGCAT");
    let q = enc("AGCTAcCGCAT");
    let mut ez = KswExtz::new();
    ksw_extz(
        q.len() as i32,
        &q,
        t.len() as i32,
        &t,
        5,
        &mat,
        2,
        1,
        -1,
        -1,
        0,
        &mut ez,
    );
    assert_eq!(ez.cigar_string(), "2D7M2D4M"); // verified against ksw2-test (C)
}
