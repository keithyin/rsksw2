#!/usr/bin/env bash
# Regenerate tests/expected/fuzz.txt from the C reference implementation.
#
# tests/c/fuzz_ref.c runs 400 randomised cases (parameters, band widths,
# z-drop, flags -- including ones cli.c cannot set --) through every KSW2
# function and prints all outputs. tests/fuzz_against_c.rs replays the exact
# same cases with the Rust port and requires identical output.
#
# NOTE: the oracle is built from patched copies of the C sources in which the
# backtrack matrices (`z` in the scalar files, `mem2` in the SSE files) are
# allocated with kcalloc() instead of kmalloc(). The C library leaves those
# bytes uninitialised; cells outside the computed band are then read during
# backtracking, so the C CIGAR depends on heap garbage and is not even
# self-consistent between runs (observed with ksw_gg2_sse). kcalloc() defines
# those bytes as zero, which is exactly what the Rust port (zeroed Vec) does,
# and it does not affect any cell the DP actually computes.
set -eu
cd "$(dirname "$0")/.."

OUT=tests/expected/fuzz.txt
mkdir -p "$(dirname "$OUT")"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp ksw2/ksw2_*.c ksw2/ksw2.h "$TMP/"
python3 - "$TMP" <<'EOF'
import re, sys, pathlib
d = pathlib.Path(sys.argv[1])
for f in d.glob('ksw2_*.c'):
    s = f.read_text()
    s2 = re.sub(r'(z|mem2) = \(uint8_t\*\)kmalloc\(km, ([^;]*)\);',
                r'\1 = (uint8_t*)kcalloc(km, 1, \2);', s)
    if s2 != s:
        f.write_text(s2)
        print('patched', f.name)
EOF

gcc -O2 -msse4.1 -I "$TMP" -o "$TMP/ksw2_fuzz_ref" tests/c/fuzz_ref.c \
    "$TMP"/ksw2_gg.c "$TMP"/ksw2_gg2.c "$TMP"/ksw2_gg2_sse.c \
    "$TMP"/ksw2_extz.c "$TMP"/ksw2_extz2_sse.c \
    "$TMP"/ksw2_extd.c "$TMP"/ksw2_extd2_sse.c \
    "$TMP"/ksw2_extf2_sse.c "$TMP"/ksw2_exts2_sse.c

"$TMP/ksw2_fuzz_ref" > "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT") lines)"
