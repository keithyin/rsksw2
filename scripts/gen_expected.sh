#!/usr/bin/env bash
# Regenerate tests/expected/t1_q1.txt from the C reference implementation.
#
# The file is used by tests/consistency.rs, which replays exactly the same
# cases with the Rust binary and requires a byte-identical result.
set -eu
cd "$(dirname "$0")/.."

CBIN=ksw2/ksw2-test
[ -x "$CBIN" ] || (cd ksw2 && make -j"$(nproc)" >/dev/null)

OUT=tests/expected/t1_q1.txt
mkdir -p "$(dirname "$OUT")"

algos="gg gg2 gg2_sse extz extz2_sse extd extd2_sse extf2_sse exts2_sse test"
flagsets=(
    ""
    "-r"
    "-s"
    "-r -s"
    "-g"
    "-s -g"
    "-z 100"
    "-z 100 -r"
    "-w 30"
    "-w 30 -r"
    "-O 4,24 -E 2,1"
    "-O 6,26 -E 2,1 -z 200"
    "-A 1 -B 2 -O 2 -E 1"
    "-a"
)

: > "$OUT"
for algo in $algos; do
    for fl in "${flagsets[@]}"; do
        printf '@@ %s |%s|\n' "$algo" "$fl" >> "$OUT"
        # shellcheck disable=SC2086
        $CBIN -t $algo $fl ksw2/test/t1.fa ksw2/test/q1.fa >> "$OUT" 2>&1
    done
done
echo "wrote $OUT ($(wc -l < "$OUT") lines)"
