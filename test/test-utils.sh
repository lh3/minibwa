#!/usr/bin/env bash
# Tests for utility commands: bench, fastmap, version
# These test the debugging and performance evaluation tools.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
# Temporary directory for test output files.
# Uses $TMPDIR if set, otherwise falls back to a local subdirectory.
TMP_DIR="${TMPDIR:-$SCRIPT_DIR}/_minibwa_test_tmp.$$"
mkdir -p "$TMP_DIR"
BINARY="${REPO_DIR}/minibwa"
PASS=0
FAIL=0

cleanup() {
    rm -f "$TMP_DIR"/mb_utils_test* "$TMP_DIR"/mb_bench_* "$TMP_DIR"/mb_fastmap.out 2>/dev/null || true
}
trap cleanup EXIT

if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary not found or not executable: $BINARY"
    echo "Build it first with: cd $REPO_DIR && make"
    exit 1
fi
pass() { PASS=$((PASS + 1)); echo -e "  \033[0;32mPASS\033[0m: $1"; }
fail() { FAIL=$((FAIL + 1)); echo -e "  \033[0;31mFAIL\033[0m: $1"; }
skip() { echo -e "  \033[1;33mSKIP\033[0m: $1"; }

echo "============================================"
echo " Minibwa Utility Command Tests"
echo "============================================"
echo ""

# Create test index
echo "Setting up test index..."
"$BINARY" index "$DATA_DIR/chrM-human.fa.gz" "$TMP_DIR/mb_utils_test" > /dev/null 2>&1
if [ ! -f "$TMP_DIR/mb_utils_test.mbw" ]; then
    echo "ERROR: Failed to create test index"
    exit 1
fi

# Test 1: bench command (2a - rank2a benchmark)
echo "[UT 1] bench -b 2a (rank2a benchmark)"
output=$("$BINARY" bench -b 2a "$TMP_DIR/mb_utils_test.mbw" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -q "checksum"; then
    pass "bench 2a runs without error"
else
    fail "bench 2a failed (exit code $rc)"
fi

# Test 2: bench command (sa benchmark)
echo "[UT 2] bench -b sa (SA query benchmark)"
output=$("$BINARY" bench -b sa "$TMP_DIR/mb_utils_test.mbw" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -q "checksum"; then
    pass "bench sa runs without error"
else
    fail "bench sa failed (exit code $rc)"
fi

# Test 3: bench command (msa benchmark)
echo "[UT 3] bench -b msa (batched SA benchmark)"
output=$("$BINARY" bench -b msa "$TMP_DIR/mb_utils_test.mbw" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -q "checksum"; then
    pass "bench msa runs without error"
else
    fail "bench msa failed (exit code $rc)"
fi

# Test 4: bench with different interval sizes
echo "[UT 4] bench -b msa with different interval sizes (-v 50)"
output=$("$BINARY" bench -b msa -v 50 "$TMP_DIR/mb_utils_test.mbw" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -q "checksum"; then
    pass "bench msa -v 50 runs without error"
else
    fail "bench msa -v 50 failed (exit code $rc)"
fi

# Test 5: bench with single SA mode
echo "[UT 5] bench -b msa -1 (single SA mode)"
output=$("$BINARY" bench -b msa -1 "$TMP_DIR/mb_utils_test.mbw" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -q "checksum"; then
    pass "bench msa -1 runs without error"
else
    fail "bench msa -1 failed (exit code $rc)"
fi

# Test 6: bench help
echo "[UT 6] bench --help"
output=$("$BINARY" bench --help 2>&1)
if echo "$output" | grep -q "Usage:"; then
    pass "bench --help shows usage"
else
    fail "bench --help does not show usage"
fi

# Test 7: fastmap command
echo "[UT 7] fastmap - test seeding strategies"
fastmap_out="$TMP_DIR/mb_fastmap.out"
"$BINARY" fastmap "$TMP_DIR/mb_utils_test" "$DATA_DIR/chrM-read_1.fa.gz" > "$fastmap_out" 2>/dev/null
rc=$?
if [ "$rc" -eq 0 ]; then
    if [ -s "$fastmap_out" ]; then
        pass "fastmap runs without error"
        lines=$(wc -l < "$fastmap_out")
        pass "fastmap produces output ($lines lines)"
    else
        fail "fastmap produced empty output"
    fi
else
    fail "fastmap failed (exit code $rc)"
fi

# Test 8: fastmap with different options
echo "[UT 8] fastmap with -n (no base alignment)"
output=$("$BINARY" fastmap -n "$TMP_DIR/mb_utils_test" "$DATA_DIR/chrM-read_1.fa.gz" 2>&1)
rc=$?
if [ "$rc" -eq 0 ] && echo "$output" | grep -qE "(seed|anchor|map)"; then
    pass "fastmap -n produces output"
else
    fail "fastmap -n failed (exit code $rc)"
fi

# Test 9: version output format
echo "[UT 9] version output format"
output=$("$BINARY" version 2>&1)
if echo "$output" | grep -qE "^[0-9]+\.[0-9]+"; then
    pass "version outputs semantic version string"
else
    fail "version output unexpected: $output"
fi

# Test 10: version via map subcommand
echo "[UT 10] version flag in map subcommand"
output=$("$BINARY" map --version 2>&1)
if echo "$output" | grep -qE "[0-9]+\.[0-9]+"; then
    pass "map --version shows version"
else
    fail "map --version output unexpected"
fi

# Test 11: bench with verbose output (-p)
echo "[UT 11] bench with -p (print per-data-point results)"
"$BINARY" bench -b 2a -p "$TMP_DIR/mb_utils_test.mbw" > "$TMP_DIR/mb_bench_p.out" 2>&1
rc=$?
if [ "$rc" -eq 0 ] && [ -s "$TMP_DIR/mb_bench_p.out" ] && grep -qE "^[0-9]+$" "$TMP_DIR/mb_bench_p.out"; then
    pass "bench -p prints per-data-point results"
else
    fail "bench -p did not print results (exit code $rc)"
fi
rm -f "$TMP_DIR"/mb_bench_p.out

# Test 12: bench with checksum verification
echo "[UT 12] bench checksum consistency"
c1=$("$BINARY" bench -b 2a "$TMP_DIR/mb_utils_test.mbw" 2>&1 | grep -oE "checksum = [0-9a-f]+" | head -1)
c2=$("$BINARY" bench -b 2a "$TMP_DIR/mb_utils_test.mbw" 2>&1 | grep -oE "checksum = [0-9a-f]+" | head -1)
if [ "$c1" = "$c2" ] && [ -n "$c1" ]; then
    pass "bench checksum is consistent across runs"
else
    fail "bench checksum differs across runs"
fi

# Cleanup
echo "Cleaning up..."
rm -f "$TMP_DIR"/mb_utils_test* "$TMP_DIR"/mb_bench_* "$TMP_DIR"/mb_fastmap.out 2>/dev/null || true

echo ""
echo "============================================"
echo " Utility Command Test Summary"
echo "============================================"
echo -e " \033[0;32m$PASS\033[0m passed"
echo -e " \033[0;31m$FAIL\033[0m failed"
echo "============================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
