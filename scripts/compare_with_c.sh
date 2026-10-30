#!/usr/bin/env bash
# Compare the Rust port against the original C implementation on the test data.
#
#   ./scripts/compare_with_c.sh            # small + medium data sets
#   ./scripts/compare_with_c.sh large      # also the 16.5 kb / 50 kb pair
#
# The C binary is built from ./ksw2 (make) and compared byte for byte with the
# Rust binary. Any difference is reported and makes the script exit non-zero.
set -u
cd "$(dirname "$0")/.."

CBIN=ksw2/ksw2-test
RBIN=target/release/ksw2-test

if [ ! -x "$CBIN" ]; then
    echo "building the C reference implementation..."
    (cd ksw2 && make -j"$(nproc)" >/dev/null)
fi
cargo build --release --quiet || exit 1

datasets=(
    "ksw2/test/t1.fa ksw2/test/q1.fa"
    "ksw2/test/t2.fa.gz ksw2/test/q2.fa.gz"
)
if [ "${1:-}" = "large" ]; then
    datasets+=("ksw2/test/MT-human.fa ksw2/test/MT-orang.fa")
fi

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

fail=0
n=0
for ds in "${datasets[@]}"; do
    for algo in $algos; do
        for fl in "${flagsets[@]}"; do
            # shellcheck disable=SC2086
            a=$($CBIN -t $algo $fl $ds 2>&1)
            b=$($RBIN -t $algo $fl $ds 2>&1)
            n=$((n + 1))
            if [ "$a" != "$b" ]; then
                echo "DIFF: $algo [$fl] $ds"
                echo "  C: $(echo "$a" | head -2 | tr '\n' '|')"
                echo "  R: $(echo "$b" | head -2 | tr '\n' '|')"
                fail=$((fail + 1))
            fi
        done
    done
done
echo "compared $n cases, $fail mismatches"
[ "$fail" -eq 0 ]
