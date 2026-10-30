//! Port of `cli.c`: a drop-in replacement for the `ksw2-test` binary.
//!
//! Output format and command line options are identical to the C program, so
//! the two binaries can be compared byte for byte.

use std::fs;
use std::io::Write;

use flate2::read::MultiGzDecoder;
use rsksw2::*;

/// Same encoding table as `cli.c` (0/1/2/3 for A/C/G/T, 4 for anything else).
const NUCLEOTIDE_ENCODING_TABLE: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;
    table
};

struct Record {
    name: String,
    sequence: String,
}

/// Minimal re-implementation of `kseq.h` for FASTA/FASTQ input.
fn read_records(path: &str) -> std::io::Result<Vec<Record>> {
    let raw_bytes = fs::read(path)?;
    let contents = if raw_bytes.len() >= 2 && raw_bytes[0] == 0x1f && raw_bytes[1] == 0x8b {
        let mut decompressed = Vec::new();
        std::io::Read::read_to_end(&mut MultiGzDecoder::new(&raw_bytes[..]), &mut decompressed)?;
        decompressed
    } else {
        raw_bytes
    };
    let mut records = Vec::new();
    let mut position = 0usize;
    let data_len = contents.len();
    while position < data_len {
        // skip empty lines between records
        while position < data_len && (contents[position] == b'\n' || contents[position] == b'\r') {
            position += 1;
        }
        if position >= data_len {
            break;
        }
        let is_fastq = contents[position] == b'@';
        if contents[position] != b'>' && contents[position] != b'@' {
            // garbage; skip the line
            while position < data_len && contents[position] != b'\n' {
                position += 1;
            }
            continue;
        }
        position += 1; // skip '>' / '@'
        let name_start = position;
        while position < data_len
            && contents[position] != b'\n'
            && contents[position] != b'\r'
            && !contents[position].is_ascii_whitespace()
        {
            position += 1;
        }
        let name = String::from_utf8_lossy(&contents[name_start..position]).into_owned();
        while position < data_len && contents[position] != b'\n' {
            position += 1;
        }
        if position < data_len {
            position += 1;
        }
        let mut sequence_buffer = Vec::new();
        let mut sequence_length = 0usize;
        loop {
            if position >= data_len {
                break;
            }
            let byte = contents[position];
            if byte == b'>' || (byte == b'@' && is_fastq && sequence_length > 0) {
                break;
            }
            if byte == b'+' && is_fastq && sequence_length > 0 {
                // FASTQ quality block: skip the '+' line and sequence_length bytes
                while position < data_len && contents[position] != b'\n' {
                    position += 1;
                }
                if position < data_len {
                    position += 1;
                }
                let mut quality_bytes = 0usize;
                while position < data_len && quality_bytes < sequence_length {
                    if contents[position] != b'\n' && contents[position] != b'\r' {
                        quality_bytes += 1;
                    }
                    position += 1;
                }
                break;
            }
            if byte == b'\n' || byte == b'\r' {
                position += 1;
                continue;
            }
            sequence_buffer.push(byte);
            sequence_length += 1;
            position += 1;
        }
        records.push(Record {
            name,
            sequence: String::from_utf8_lossy(&sequence_buffer).into_owned(),
        });
    }
    Ok(records)
}

fn encode_sequence(sequence: &str) -> Vec<u8> {
    sequence
        .bytes()
        .map(|byte| NUCLEOTIDE_ENCODING_TABLE[byte as usize])
        .collect()
}

struct Options {
    gap_open: i8,
    gap_extension: i8,
    gap_open2: i8,
    gap_extension2: i8,
    band_width: i32,
    z_drop: i32,
    flags: i32,
    alphabet_size: i8,
    score_matrix: [i8; 25],
}

fn run_alignment(
    algorithm: &str,
    result: &mut KswExtz,
    query_sequence: &str,
    target_sequence: &str,
    options: &Options,
) {
    // equivalent of the field resets at the top of global_aln() in cli.c
    result.max_q = -1;
    result.max_t = -1;
    result.mqe_t = -1;
    result.mte_q = -1;
    result.max = 0;
    result.mqe = KSW_NEG_INF;
    result.mte = KSW_NEG_INF;
    result.cigar.clear();

    let encoded_query = encode_sequence(query_sequence);
    let encoded_target = encode_sequence(target_sequence);
    let query_len = encoded_query.len() as i32;
    let target_len = encoded_target.len() as i32;
    let alphabet_size = options.alphabet_size;
    let score_matrix = &options.score_matrix[..];
    let gap_open = options.gap_open;
    let gap_extension = options.gap_extension;
    let gap_open2 = options.gap_open2;
    let gap_extension2 = options.gap_extension2;
    let band_width = options.band_width;
    let z_drop = options.z_drop;
    let flags = options.flags;
    let score_only = flags & KSW_EZ_SCORE_ONLY != 0;

    match algorithm {
        "gg" => {
            let cigar = if score_only { None } else { Some(&mut result.cigar) };
            result.score = ksw_gg(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                band_width,
                cigar,
            );
        }
        "gg2" => {
            let cigar = if score_only { None } else { Some(&mut result.cigar) };
            result.score = ksw_gg2(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                band_width,
                cigar,
            );
        }
        "gg2_sse" => {
            result.score = ksw_gg2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                band_width,
                &mut result.cigar,
            );
        }
        "extz" => {
            ksw_extz(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                band_width,
                z_drop,
                flags,
                result,
            );
        }
        "extz2_sse" => {
            ksw_extz2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                band_width,
                z_drop,
                0,
                flags,
                result,
            );
        }
        "extd" => {
            ksw_extd(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                gap_open2,
                gap_extension2,
                band_width,
                z_drop,
                flags,
                result,
            );
        }
        "extd2_sse" => {
            ksw_extd2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extension,
                gap_open2,
                gap_extension2,
                band_width,
                z_drop,
                0,
                flags,
                result,
            );
        }
        "extf2_sse" => {
            ksw_extf2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                score_matrix[0],
                score_matrix[1],
                gap_extension,
                band_width,
                z_drop,
                result,
            );
        }
        "exts2_sse" => {
            let mut splice_score_matrix = [0i8; 25];
            gen_simple_mat(5, &mut splice_score_matrix, 1, 2);
            ksw_exts2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                &splice_score_matrix,
                2,
                1,
                32,
                4,
                z_drop,
                0,
                flags | KSW_EZ_SPLICE_FOR,
                &[],
                result,
            );
        }
        "test" => {
            ksw_extd2_sse(
                query_len,
                &encoded_query,
                target_len,
                &encoded_target,
                alphabet_size,
                score_matrix,
                4,
                2,
                24,
                1,
                751,
                400,
                0,
                8,
                result,
            );
        }
        _ => {
            eprintln!("ERROR: can't find algorithm '{algorithm}'");
            std::process::exit(1);
        }
    }
}

fn print_alignment(writer: &mut impl Write, target_name: &str, query_name: &str, result: &KswExtz) {
    let _ = write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}",
        target_name,
        query_name,
        result.score,
        result.max as i32,
        result.max_t,
        result.max_q
    );
    if !result.cigar.is_empty() {
        let _ = writer.write_all(b"\t");
        for &encoded in &result.cigar {
            let operation = (encoded & 0xf) as usize;
            // cli.c prints "MID"[op], i.e. a NUL byte for op 3 (N) and an
            // out-of-bounds read for op >= 4; we reproduce the former and use
            // '='/'X' for the EQX operators.
            let operator_byte = match operation {
                0 => b'M',
                1 => b'I',
                2 => b'D',
                3 => 0u8,
                7 => b'=',
                8 => b'X',
                _ => b'?',
            };
            let _ = write!(writer, "{}", encoded >> 4);
            let _ = writer.write_all(&[operator_byte]);
        }
    }
    let _ = writer.write_all(b"\n");
}

fn usage(
    algorithm: &str,
    z_drop: i32,
    match_score: i8,
    mismatch_penalty: i8,
    gap_open: i8,
    gap_extension: i8,
    gap_open2: i8,
    gap_extension2: i8,
) {
    eprintln!("Usage: ksw2-test [options] <DNA-target> <DNA-query>");
    eprintln!("Options:");
    eprintln!(
        "  -t STR        algorithm: gg, gg2, gg2_sse, extz, extz2_sse, extd, extd2_sse [{}]",
        algorithm
    );
    eprintln!("  -R INT        repeat INT times (for benchmarking) [1]");
    eprintln!("  -w INT        band width [inf]");
    eprintln!("  -z INT        Z-drop [{}]", z_drop);
    eprintln!("  -r            gap right alignment");
    eprintln!("  -s            score only");
    eprintln!("  -A INT        match score [{}]", match_score);
    eprintln!("  -B INT        mismatch penalty [{}]", mismatch_penalty);
    eprintln!("  -O INT[,INT]  gap open penalty [{},{}]", gap_open, gap_open2);
    eprintln!(
        "  -E INT[,INT]  gap extension penalty [{},{}]",
        gap_extension, gap_extension2
    );
    eprintln!("  -a            all vs all");
}

fn parse_penalty_pair(spec: &str) -> (i32, Option<i32>) {
    let mut parts = spec.split(',');
    let first = parts.next().unwrap_or("0").trim().parse::<i32>().unwrap_or(0);
    let second = parts.next().map(|part| part.trim().parse::<i32>().unwrap_or(0));
    (first, second)
}

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    let mut match_score: i8 = 2;
    let mut mismatch_penalty: i8 = 4;
    let mut gap_open: i8 = 4;
    let mut gap_extension: i8 = 2;
    let mut gap_open2: i8 = 13;
    let mut gap_extension2: i8 = 1;
    let mut pairwise_mode = true;
    let mut band_width = -1i32;
    let mut flags = 0i32;
    let mut repeat_count = 1i32;
    let mut z_drop = -1i32;
    let mut algorithm = String::from("extd");

    let mut position = 1usize;
    let mut operands: Vec<String> = Vec::new();
    while position < arguments.len() {
        let argument = &arguments[position];
        if argument.len() >= 2 && argument.starts_with('-') && argument != "-" {
            // supports "-t extd" as well as "-textd"
            let option = argument.as_bytes()[1] as char;
            let mut inline_value = if argument.len() > 2 {
                Some(argument[2..].to_string())
            } else {
                None
            };
            let mut next_value = |position: &mut usize| -> String {
                if let Some(value) = inline_value.take() {
                    return value;
                }
                *position += 1;
                arguments.get(*position).cloned().unwrap_or_default()
            };
            match option {
                't' => algorithm = next_value(&mut position),
                'w' => band_width = next_value(&mut position).trim().parse().unwrap_or(0),
                'R' => repeat_count = next_value(&mut position).trim().parse().unwrap_or(1),
                'z' => z_drop = next_value(&mut position).trim().parse().unwrap_or(0),
                'A' => match_score = next_value(&mut position).trim().parse().unwrap_or(0),
                'B' => mismatch_penalty = next_value(&mut position).trim().parse().unwrap_or(0),
                'O' => {
                    let (first, second) = parse_penalty_pair(&next_value(&mut position));
                    gap_open = first as i8;
                    gap_open2 = second.unwrap_or(first) as i8;
                }
                'E' => {
                    let (first, second) = parse_penalty_pair(&next_value(&mut position));
                    gap_extension = first as i8;
                    gap_extension2 = second.unwrap_or(first) as i8;
                }
                'r' => flags |= KSW_EZ_RIGHT,
                's' => flags |= KSW_EZ_SCORE_ONLY,
                'g' => flags |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP,
                'K' => {} // no_kalloc: the Rust port has no memory pool
                'a' => pairwise_mode = false,
                _ => {}
            }
        } else {
            operands.push(argument.clone());
        }
        position += 1;
    }

    if operands.len() < 2 {
        usage(
            &algorithm,
            z_drop,
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extension,
            gap_open2,
            gap_extension2,
        );
        std::process::exit(1);
    }

    let mut score_matrix = [0i8; 25];
    gen_simple_mat(5, &mut score_matrix, match_score, -mismatch_penalty);
    let options = Options {
        gap_open,
        gap_extension,
        gap_open2,
        gap_extension2,
        band_width,
        z_drop,
        flags,
        alphabet_size: 5,
        score_matrix,
    };

    let stdout = std::io::stdout();
    let mut writer = std::io::BufWriter::new(stdout.lock());
    let mut result = KswExtz::new();

    let target_file = read_records(&operands[0]);
    let query_file = read_records(&operands[1]);
    match (target_file, query_file) {
        (Err(_), Err(_)) => {
            run_alignment(&algorithm, &mut result, &operands[1], &operands[0], &options);
            print_alignment(&mut writer, "first", "second", &result);
        }
        (Ok(target_records), Ok(query_records)) => {
            if pairwise_mode {
                for (target_record, query_record) in target_records.iter().zip(query_records.iter()) {
                    for _ in 0..repeat_count {
                        run_alignment(
                            &algorithm,
                            &mut result,
                            &query_record.sequence,
                            &target_record.sequence,
                            &options,
                        );
                    }
                    print_alignment(&mut writer, &target_record.name, &query_record.name, &result);
                }
            } else {
                for query_record in query_records.iter() {
                    for target_record in target_records.iter() {
                        run_alignment(
                            &algorithm,
                            &mut result,
                            &query_record.sequence,
                            &target_record.sequence,
                            &options,
                        );
                        print_alignment(&mut writer, &target_record.name, &query_record.name, &result);
                    }
                }
            }
        }
        _ => {
            // same as cli.c: when only one of the two files can be opened, do nothing
        }
    }
    let _ = writer.flush();
}
