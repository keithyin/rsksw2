#!/usr/bin/env bash
# Compare the run time of the Rust port with the C implementation.
#
#   ./scripts/bench.sh [target.fa] [query.fa] [repeats]
#
# Default data set: ksw2/test/MT-human.fa vs ksw2/test/MT-orang.fa (~16.5 kb),
# the pair used by the KSW2 README. Times are wall clock seconds, best of 3.
set -eu
cd "$(dirname "$0")/.."

TG=${1:-ksw2/test/MT-human.fa}
QR=${2:-ksw2/test/MT-orang.fa}
REP=${3:-1}

CBIN=ksw2/ksw2-test
RBIN=target/release/ksw2-test
[ -x "$CBIN" ] || (cd ksw2 && make -j"$(nproc)" >/dev/null)
RUSTFLAGS="${RUSTFLAGS:--C target-cpu=native}" cargo build --release --quiet

now() { date +%s%N; }

# best of 3 wall clock runs, printed in nanoseconds
best() {
    local bin=$1 algo=$2 flags=$3 best=999999999999
    for _ in 1 2 3; do
        local t0 t1 dt
        t0=$(now)
        # shellcheck disable=SC2086
        $bin -t "$algo" -R "$REP" $flags "$TG" "$QR" >/dev/null 2>&1
        t1=$(now)
        dt=$((t1 - t0))
        [ "$dt" -lt "$best" ] && best=$dt
    done
    echo "$best"
}

printf '%-12s %-6s %10s %10s %8s\n' algo mode "C (s)" "Rust (s)" "Rust/C"
for algo in ${ALGOS:-extz2_sse extd2_sse extz extd gg2_sse gg2 gg}; do
    for flags in "-s" ""; do
        c=$(best "$CBIN" "$algo" "$flags")
        r=$(best "$RBIN" "$algo" "$flags")
        awk -v a="$algo" -v f="${flags:--}" -v c="$c" -v r="$r" \
            'BEGIN { printf "%-12s %-6s %10.3f %10.3f %8.2f\n", a, f, c/1e9, r/1e9, r/c }'
    done
done
