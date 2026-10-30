//! Randomised comparison of every Rust KSW2 function against the C library.
//!
//! `scripts/gen_fuzz_expected.sh` runs `tests/c/fuzz_ref.c` (400 randomised
//! cases exercising every flag the public API can take, including ones
//! `cli.c` cannot set) and stores its output in `tests/expected/fuzz.txt`.
//! This test reproduces the exact same cases with the Rust port and requires
//! byte-identical output.

use std::fmt::Write;
use rsksw2::*;

const EXPECTED: &str = include_str!("expected/fuzz.txt");

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
}

const W_LIST: [i32; 7] = [-1, 0, 1, 2, 5, 20, 100];
const Z_LIST: [i32; 5] = [-1, 0, 5, 15, 100];
const F_LIST: [i32; 12] = [0x000, 0x001, 0x040, 0x002, 0x080, 0x082, 0x003, 0x004, 0x006, 0x044, 0x018, 0x800];
const FS_LIST: [i32; 9] = [0x100, 0x200, 0x104, 0x204, 0x500, 0x600, 0x101, 0x201, 0x000];

struct Params {
    qlen: usize,
    tlen: usize,
    q: Vec<u8>,
    t: Vec<u8>,
    junc: Vec<u8>,
    mat: [i8; 25],
    a: i8,
    b: i8,
    gq: i8,
    ge: i8,
    gq2: i8,
    ge2: i8,
    w: i32,
    zdrop: i32,
    t_bonus: i32,
    flag: i32,
    sflag: i32,
    xdrop: i32,
    noncan: i8,
    jb: i8,
}

fn gen_case(id: u64) -> Params {
    let mut rng = Rng(id.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1));
    let qlen = 1 + (rng.next() % 60) as usize;
    let tlen = 1 + (rng.next() % 60) as usize;
    let q: Vec<u8> = (0..qlen).map(|_| (rng.next() % 5) as u8).collect();
    let t: Vec<u8> = (0..tlen).map(|_| (rng.next() % 5) as u8).collect();
    let a = 1 + (rng.next() % 4) as i8;
    let b = -((1 + (rng.next() % 5)) as i8);
    let mut mat = [0i8; 25];
    gen_simple_mat(5, &mut mat, a, b);
    let gq = 2 + (rng.next() % 9) as i8;
    let ge = 1 + (rng.next() % 4) as i8;
    let gq2 = 10 + (rng.next() % 11) as i8;
    let ge2 = 1 + (rng.next() % 4) as i8;
    // A band that does not cover the alignment makes the C *scalar*
    // implementations read out of bounds in ksw_backtrack() (they lack the
    // off_end[] guard of the SSE versions) -- undefined behaviour, not a
    // meaningful oracle. Keep regenerating until the band covers the ends.
    let w = loop {
        let w = W_LIST[(rng.next() % 7) as usize];
        if !(w >= 0 && (qlen as i32 - tlen as i32).abs() > w) {
            break w;
        }
    };
    let zdrop = Z_LIST[(rng.next() % 5) as usize];
    let t_bonus = (rng.next() % 6) as i32;
    let flag = F_LIST[(rng.next() % 12) as usize];
    let sflag = FS_LIST[(rng.next() % 9) as usize];
    let xdrop = if rng.next() % 2 != 0 { -1 } else { 8 + (rng.next() % 25) as i32 };
    let noncan = (rng.next() % 8) as i8;
    let jb = (rng.next() % 15) as i8;
    let junc: Vec<u8> = (0..tlen).map(|_| rng.next() as u8).collect();
    Params {
        qlen, tlen, q, t, junc, mat, a, b, gq, ge, gq2, ge2,
        w, zdrop, t_bonus, flag, sflag, xdrop, noncan, jb,
    }
}

fn print_cigar(out: &mut String, c: &[u32]) {
    out.push_str(" [");
    for (k, &x) in c.iter().enumerate() {
        if k > 0 {
            out.push(' ');
        }
        let _ = write!(out, "{}:{}", x >> 4, x & 0xf);
    }
    out.push(']');
}

fn print_ez(out: &mut String, name: &str, ez: &KswExtz) {
    let _ = write!(
        out,
        "{} sc={} mqe={} mqe_t={} mte={} mte_q={} max={} max_t={} max_q={} z={} re={}",
        name, ez.score, ez.mqe, ez.mqe_t, ez.mte, ez.mte_q,
        ez.max, ez.max_t, ez.max_q, ez.zdropped as i32, ez.reach_end
    );
    print_cigar(out, &ez.cigar);
    out.push('\n');
}

#[test]
fn matches_c_library_on_randomised_cases() {
    let mut rebuilt = String::new();
    for id in 0u64..2000 {
        let p = gen_case(id);
        let _ = write!(
            rebuilt,
            "== {} ", id
        );
        for &x in &p.q {
            rebuilt.push(b"ACGTN"[x as usize] as char);
        }
        rebuilt.push(' ');
        for &x in &p.t {
            rebuilt.push(b"ACGTN"[x as usize] as char);
        }
        let _ = write!(
            rebuilt,
            " q={} e={} q2={} e2={} w={} z={} T={} f={} sf={} x={} nc={} jb={}\n",
            p.gq, p.ge, p.gq2, p.ge2, p.w, p.zdrop, p.t_bonus, p.flag, p.sflag,
            p.xdrop, p.noncan, p.jb
        );

        let mut ci = Vec::new();
        let sc = ksw_gg(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.w, Some(&mut ci),
        );
        let _ = write!(rebuilt, "gg {}", sc);
        print_cigar(&mut rebuilt, &ci);
        rebuilt.push('\n');

        let mut ci = Vec::new();
        let sc = ksw_gg2(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.w, Some(&mut ci),
        );
        let _ = write!(rebuilt, "gg2 {}", sc);
        print_cigar(&mut rebuilt, &ci);
        rebuilt.push('\n');

        let mut ci = Vec::new();
        let sc = ksw_gg2_sse(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.w, &mut ci,
        );
        let _ = write!(rebuilt, "gg2_sse {}", sc);
        print_cigar(&mut rebuilt, &ci);
        rebuilt.push('\n');

        let mut ez = KswExtz::new();
        ksw_extz(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge,
            p.w, p.zdrop, p.flag, &mut ez,
        );
        print_ez(&mut rebuilt, "extz", &ez);

        let mut ez = KswExtz::new();
        ksw_extz2_sse(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge,
            p.w, p.zdrop, p.t_bonus, p.flag, &mut ez,
        );
        print_ez(&mut rebuilt, "extz2_sse", &ez);

        let mut ez = KswExtz::new();
        ksw_extd(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.gq2, p.ge2,
            p.w, p.zdrop, p.flag, &mut ez,
        );
        print_ez(&mut rebuilt, "extd", &ez);

        let mut ez = KswExtz::new();
        ksw_extd2_sse(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.gq2, p.ge2,
            p.w, p.zdrop, p.t_bonus, p.flag, &mut ez,
        );
        print_ez(&mut rebuilt, "extd2_sse", &ez);

        let mut ez = KswExtz::new();
        ksw_extf2_sse(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, p.a, p.b, p.ge, p.w, p.xdrop, &mut ez,
        );
        print_ez(&mut rebuilt, "extf2_sse", &ez);

        let mut ez = KswExtz::new();
        ksw_exts2_sse(
            p.qlen as i32, &p.q, p.tlen as i32, &p.t, 5, &p.mat, p.gq, p.ge, p.gq2,
            p.noncan, p.zdrop, p.jb, p.sflag, &p.junc, &mut ez,
        );
        print_ez(&mut rebuilt, "exts2_sse", &ez);
    }

    if rebuilt != EXPECTED {
        for (i, (a, b)) in EXPECTED.lines().zip(rebuilt.lines()).enumerate() {
            if a != b {
                panic!("mismatch on line {}:\n  C: {a}\n  R: {b}", i + 1);
            }
        }
        panic!("length differs: C {} vs Rust {}", EXPECTED.len(), rebuilt.len());
    }
}
