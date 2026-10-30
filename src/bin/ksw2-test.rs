//! Port of `cli.c`: a drop-in replacement for the `ksw2-test` binary.
//!
//! Output format and command line options are identical to the C program, so
//! the two binaries can be compared byte for byte.

use std::fs;
use std::io::Write;

use flate2::read::MultiGzDecoder;
use rsksw2::*;

/// Same encoding table as `cli.c` (0/1/2/3 for A/C/G/T, 4 for anything else).
const SEQ_NT4_TABLE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

struct Record {
    name: String,
    seq: String,
}

/// Minimal re-implementation of `kseq.h` for FASTA/FASTQ input.
fn read_records(path: &str) -> std::io::Result<Vec<Record>> {
    let raw = fs::read(path)?;
    let data = if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::new();
        std::io::Read::read_to_end(&mut MultiGzDecoder::new(&raw[..]), &mut out)?;
        out
    } else {
        raw
    };
    let mut recs = Vec::new();
    let mut i = 0usize;
    let n = data.len();
    while i < n {
        // skip empty lines between records
        while i < n && (data[i] == b'\n' || data[i] == b'\r') {
            i += 1;
        }
        if i >= n {
            break;
        }
        let fastq = data[i] == b'@';
        if data[i] != b'>' && data[i] != b'@' {
            // garbage; skip the line
            while i < n && data[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        i += 1; // skip '>' / '@'
        let name_start = i;
        while i < n && data[i] != b'\n' && data[i] != b'\r' && !data[i].is_ascii_whitespace() {
            i += 1;
        }
        let name = String::from_utf8_lossy(&data[name_start..i]).into_owned();
        while i < n && data[i] != b'\n' {
            i += 1;
        }
        if i < n {
            i += 1;
        }
        let mut seq = Vec::new();
        let mut seq_len = 0usize;
        loop {
            if i >= n {
                break;
            }
            let c = data[i];
            if c == b'>' || (c == b'@' && fastq && seq_len > 0) {
                break;
            }
            if c == b'+' && fastq && seq_len > 0 {
                // FASTQ quality block: skip the '+' line and seq_len bytes
                while i < n && data[i] != b'\n' {
                    i += 1;
                }
                if i < n {
                    i += 1;
                }
                let mut q = 0usize;
                while i < n && q < seq_len {
                    if data[i] != b'\n' && data[i] != b'\r' {
                        q += 1;
                    }
                    i += 1;
                }
                break;
            }
            if c == b'\n' || c == b'\r' {
                i += 1;
                continue;
            }
            seq.push(c);
            seq_len += 1;
            i += 1;
        }
        recs.push(Record {
            name,
            seq: String::from_utf8_lossy(&seq).into_owned(),
        });
    }
    Ok(recs)
}

fn encode(seq: &str) -> Vec<u8> {
    seq.bytes().map(|c| SEQ_NT4_TABLE[c as usize]).collect()
}

struct Opts {
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
    w: i32,
    zdrop: i32,
    flag: i32,
    m: i8,
    mat: [i8; 25],
}

fn global_aln(algo: &str, ez: &mut KswExtz, qseq_: &str, tseq_: &str, o: &Opts) {
    // equivalent of the field resets at the top of global_aln() in cli.c
    ez.max_q = -1;
    ez.max_t = -1;
    ez.mqe_t = -1;
    ez.mte_q = -1;
    ez.max = 0;
    ez.mqe = KSW_NEG_INF;
    ez.mte = KSW_NEG_INF;
    ez.cigar.clear();

    let qseq = encode(qseq_);
    let tseq = encode(tseq_);
    let qlen = qseq.len() as i32;
    let tlen = tseq.len() as i32;
    let (m, mat, q, e, q2, e2, w, zdrop, flag) =
        (o.m, &o.mat[..], o.q, o.e, o.q2, o.e2, o.w, o.zdrop, o.flag);
    let score_only = flag & KSW_EZ_SCORE_ONLY != 0;

    match algo {
        "gg" => {
            let ci = if score_only { None } else { Some(&mut ez.cigar) };
            ez.score = ksw_gg(qlen, &qseq, tlen, &tseq, m, mat, q, e, w, ci);
        }
        "gg2" => {
            let ci = if score_only { None } else { Some(&mut ez.cigar) };
            ez.score = ksw_gg2(qlen, &qseq, tlen, &tseq, m, mat, q, e, w, ci);
        }
        "gg2_sse" => {
            ez.score = ksw_gg2_sse(qlen, &qseq, tlen, &tseq, m, mat, q, e, w, &mut ez.cigar);
        }
        "extz" => {
            ksw_extz(qlen, &qseq, tlen, &tseq, m, mat, q, e, w, zdrop, flag, ez);
        }
        "extz2_sse" => {
            ksw_extz2_sse(qlen, &qseq, tlen, &tseq, m, mat, q, e, w, zdrop, 0, flag, ez);
        }
        "extd" => {
            ksw_extd(qlen, &qseq, tlen, &tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
        }
        "extd2_sse" => {
            ksw_extd2_sse(
                qlen, &qseq, tlen, &tseq, m, mat, q, e, q2, e2, w, zdrop, 0, flag, ez,
            );
        }
        "extf2_sse" => {
            ksw_extf2_sse(qlen, &qseq, tlen, &tseq, mat[0], mat[1], e, w, zdrop, ez);
        }
        "exts2_sse" => {
            let mut mat2 = [0i8; 25];
            gen_simple_mat(5, &mut mat2, 1, 2);
            ksw_exts2_sse(
                qlen,
                &qseq,
                tlen,
                &tseq,
                m,
                &mat2,
                2,
                1,
                32,
                4,
                zdrop,
                0,
                flag | KSW_EZ_SPLICE_FOR,
                &[],
                ez,
            );
        }
        "test" => {
            ksw_extd2_sse(
                qlen, &qseq, tlen, &tseq, m, mat, 4, 2, 24, 1, 751, 400, 0, 8, ez,
            );
        }
        _ => {
            eprintln!("ERROR: can't find algorithm '{algo}'");
            std::process::exit(1);
        }
    }
}

fn print_aln(out: &mut impl Write, tname: &str, qname: &str, ez: &KswExtz) {
    let _ = write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}",
        tname,
        qname,
        ez.score,
        ez.max as i32,
        ez.max_t,
        ez.max_q
    );
    if !ez.cigar.is_empty() {
        let _ = out.write_all(b"\t");
        for &c in &ez.cigar {
            let op = (c & 0xf) as usize;
            // cli.c prints "MID"[op], i.e. a NUL byte for op 3 (N) and an
            // out-of-bounds read for op >= 4; we reproduce the former and use
            // '='/'X' for the EQX operators.
            let ch = match op {
                0 => b'M',
                1 => b'I',
                2 => b'D',
                3 => 0u8,
                7 => b'=',
                8 => b'X',
                _ => b'?',
            };
            let _ = write!(out, "{}", c >> 4);
            let _ = out.write_all(&[ch]);
        }
    }
    let _ = out.write_all(b"\n");
}

fn usage(algo: &str, zdrop: i32, a: i8, b: i8, q: i8, e: i8, q2: i8, e2: i8) {
    eprintln!("Usage: ksw2-test [options] <DNA-target> <DNA-query>");
    eprintln!("Options:");
    eprintln!(
        "  -t STR        algorithm: gg, gg2, gg2_sse, extz, extz2_sse, extd, extd2_sse [{}]",
        algo
    );
    eprintln!("  -R INT        repeat INT times (for benchmarking) [1]");
    eprintln!("  -w INT        band width [inf]");
    eprintln!("  -z INT        Z-drop [{}]", zdrop);
    eprintln!("  -r            gap right alignment");
    eprintln!("  -s            score only");
    eprintln!("  -A INT        match score [{}]", a);
    eprintln!("  -B INT        mismatch penalty [{}]", b);
    eprintln!("  -O INT[,INT]  gap open penalty [{},{}]", q, q2);
    eprintln!("  -E INT[,INT]  gap extension penalty [{},{}]", e, e2);
    eprintln!("  -a            all vs all");
}

fn parse_pair(s: &str) -> (i32, Option<i32>) {
    let mut it = s.split(',');
    let first = it.next().unwrap_or("0").trim().parse::<i32>().unwrap_or(0);
    let second = it.next().map(|x| x.trim().parse::<i32>().unwrap_or(0));
    (first, second)
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut a: i8 = 2;
    let mut b: i8 = 4;
    let mut q: i8 = 4;
    let mut e: i8 = 2;
    let mut q2: i8 = 13;
    let mut e2: i8 = 1;
    let mut pair = true;
    let mut w = -1i32;
    let mut flag = 0i32;
    let mut rep = 1i32;
    let mut zdrop = -1i32;
    let mut algo = String::from("extd");

    let mut i = 1usize;
    let mut operands: Vec<String> = Vec::new();
    while i < argv.len() {
        let arg = &argv[i];
        if arg.len() >= 2 && arg.starts_with('-') && arg != "-" {
            // supports "-t extd" as well as "-textd"
            let opt = arg.as_bytes()[1] as char;
            let mut inline_val = if arg.len() > 2 { Some(arg[2..].to_string()) } else { None };
            let mut next_val = |i: &mut usize| -> String {
                if let Some(v) = inline_val.take() {
                    return v;
                }
                *i += 1;
                argv.get(*i).cloned().unwrap_or_default()
            };
            match opt {
                't' => algo = next_val(&mut i),
                'w' => w = next_val(&mut i).trim().parse().unwrap_or(0),
                'R' => rep = next_val(&mut i).trim().parse().unwrap_or(1),
                'z' => zdrop = next_val(&mut i).trim().parse().unwrap_or(0),
                'A' => a = next_val(&mut i).trim().parse().unwrap_or(0),
                'B' => b = next_val(&mut i).trim().parse().unwrap_or(0),
                'O' => {
                    let (x, y) = parse_pair(&next_val(&mut i));
                    q = x as i8;
                    q2 = y.unwrap_or(x) as i8;
                }
                'E' => {
                    let (x, y) = parse_pair(&next_val(&mut i));
                    e = x as i8;
                    e2 = y.unwrap_or(x) as i8;
                }
                'r' => flag |= KSW_EZ_RIGHT,
                's' => flag |= KSW_EZ_SCORE_ONLY,
                'g' => flag |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP,
                'K' => {} // no_kalloc: the Rust port has no memory pool
                'a' => pair = false,
                _ => {}
            }
        } else {
            operands.push(arg.clone());
        }
        i += 1;
    }

    if operands.len() < 2 {
        usage(&algo, zdrop, a, b, q, e, q2, e2);
        std::process::exit(1);
    }

    let mut mat = [0i8; 25];
    gen_simple_mat(5, &mut mat, a, -b);
    let opts = Opts {
        q,
        e,
        q2,
        e2,
        w,
        zdrop,
        flag,
        m: 5,
        mat,
    };

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    let mut ez = KswExtz::new();

    let f0 = read_records(&operands[0]);
    let f1 = read_records(&operands[1]);
    match (f0, f1) {
        (Err(_), Err(_)) => {
            global_aln(&algo, &mut ez, &operands[1], &operands[0], &opts);
            print_aln(&mut out, "first", "second", &ez);
        }
        (Ok(trecs), Ok(qrecs)) => {
            if pair {
                for (t, qr) in trecs.iter().zip(qrecs.iter()) {
                    for _ in 0..rep {
                        global_aln(&algo, &mut ez, &qr.seq, &t.seq, &opts);
                    }
                    print_aln(&mut out, &t.name, &qr.name, &ez);
                }
            } else {
                for qr in qrecs.iter() {
                    for t in trecs.iter() {
                        global_aln(&algo, &mut ez, &qr.seq, &t.seq, &opts);
                        print_aln(&mut out, &t.name, &qr.name, &ez);
                    }
                }
            }
        }
        _ => {
            // same as cli.c: when only one of the two files can be opened, do nothing
        }
    }
    let _ = out.flush();
}
