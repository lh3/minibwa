#!/usr/bin/env bash
# API example tests for minibwa
# Tests the C API examples in api-test/

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
# Temporary directory for test output files.
# Uses $TMPDIR if set, otherwise falls back to a local subdirectory.
TMP_DIR="${TMPDIR:-$SCRIPT_DIR}/_minibwa_test_tmp"
mkdir -p "$TMP_DIR"
BINARY="${REPO_DIR}/minibwa"
API_DIR="$REPO_DIR/api-test"
PASS=0
FAIL=0

cleanup() {
    rm -f "$TMP_DIR"/mb_api_test.l2b "$TMP_DIR"/mb_api_test.mbw "$TMP_DIR"/mb_api_one_out.txt "$TMP_DIR"/mb_api_batch_out.txt 2>/dev/null || true
    rm -f "$API_DIR/mbmap-one" "$API_DIR/mbmap-batch" 2>/dev/null || true
}
trap cleanup EXIT

if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary not found or not executable: $BINARY"
    echo "Build it first with: cd $REPO_DIR && make"
    exit 1
fi

pass() { PASS=$((PASS + 1)); echo -e "  \033[0;32mPASS\033[0m: $1"; }
fail() { FAIL=$((FAIL + 1)); echo -e "  \033[0;31mFAIL\033[0m: $1"; }

echo "============================================"
echo " Minibwa API Example Tests"
echo "============================================"
echo ""

# Verify API examples are built (should be built by CI before tests)
echo "[API 1] Verifying API examples"
if [ -f "$API_DIR/mbmap-one" ] && [ -f "$API_DIR/mbmap-batch" ]; then
    pass "API examples found"
else
    fail "API examples not found (build with: make -C api-test)"
    exit 1
fi

# Create test index for API
echo "[API 2] Setting up test index"
"$BINARY" index "$DATA_DIR/chrM-human.fa.gz" $TMP_DIR/mb_api_test > /dev/null 2>&1
if [ -f "$TMP_DIR/mb_api_test.l2b" ] && [ -f "$TMP_DIR/mb_api_test.mbw" ]; then
    pass "Test index created for API tests"
else
    fail "Failed to create test index for API"
    exit 1
fi

# Test ex-one.c (single read alignment)
echo "[API 3] mbmap-one: single read alignment"
"$API_DIR/mbmap-one" $TMP_DIR/mb_api_test "$DATA_DIR/chrM-read_1.fa.gz" > $TMP_DIR/mb_api_one_out.txt 2>/dev/null
if [ -s "$TMP_DIR/mb_api_one_out.txt" ]; then
    pass "mbmap-one produced output"
    line_count=$(wc -l < $TMP_DIR/mb_api_one_out.txt)
    if [ "$line_count" -gt 0 ]; then
        pass "mbmap-one produced $line_count alignment lines"
    else
        fail "mbmap-one produced empty alignment output"
    fi
else
    fail "mbmap-one produced no output"
fi

# Test ex-batch.c (batch alignment)
echo "[API 4] mbmap-batch: batch alignment"
"$API_DIR/mbmap-batch" $TMP_DIR/mb_api_test "$DATA_DIR/chrM-read_1.fa.gz" > $TMP_DIR/mb_api_batch_out.txt 2>/dev/null
if [ -s "$TMP_DIR/mb_api_batch_out.txt" ]; then
    pass "mbmap-batch produced output"
    line_count=$(wc -l < $TMP_DIR/mb_api_batch_out.txt)
    if [ "$line_count" -gt 0 ]; then
        pass "mbmap-batch produced $line_count alignment lines"
    else
        fail "mbmap-batch produced empty alignment output"
    fi
else
    fail "mbmap-batch produced no output"
fi

# Test that single-read and batch outputs are similar
echo "[API 5] Consistency: single vs batch alignment"
# Both should produce alignments for the same reads
single_count=$(wc -l < $TMP_DIR/mb_api_one_out.txt)
batch_count=$(wc -l < $TMP_DIR/mb_api_batch_out.txt)
if [ "$single_count" -gt 0 ] && [ "$batch_count" -gt 0 ]; then
    # Allow some difference due to batching behavior
    diff_pct=$(( (single_count > batch_count ? single_count - batch_count : batch_count - single_count) * 100 / (single_count > batch_count ? single_count : batch_count) ))
    if [ "$diff_pct" -lt 50 ]; then
        pass "Single and batch produce comparable results (diff: ${diff_pct}%)"
    else
        fail "Single and batch results differ too much (diff: ${diff_pct}%)"
    fi
else
    fail "Could not compare single vs batch (one produced no output)"
fi

# Test API with no arguments
echo "[API 6] API examples show usage with no arguments"
rc=0; "$API_DIR/mbmap-one" 2>/dev/null || rc=$?
if [ "$rc" -ne 0 ]; then
    pass "mbmap-one exits non-zero with no args"
else
    fail "mbmap-one should exit non-zero with no args"
fi

rc=0; "$API_DIR/mbmap-batch" 2>/dev/null || rc=$?
if [ "$rc" -ne 0 ]; then
    pass "mbmap-batch exits non-zero with no args"
else
    fail "mbmap-batch should exit non-zero with no args"
fi

# Test API with missing index
echo "[API 7] API handles missing index gracefully"
rc=0; "$API_DIR/mbmap-one" $TMP_DIR/nonexistent_index $TMP_DIR/mb_test_noop.fa 2>/dev/null || rc=$?
if [ "$rc" -ne 0 ]; then
    pass "mbmap-one exits non-zero with missing index"
else
    fail "mbmap-one should handle missing index"
fi

# Cleanup
echo ""
echo "Cleaning up..."
rm -f $TMP_DIR/mb_api_test.l2b $TMP_DIR/mb_api_test.mbw $TMP_DIR/mb_api_one_out.txt $TMP_DIR/mb_api_batch_out.txt 2>/dev/null || true
rm -f "$API_DIR/mbmap-one" "$API_DIR/mbmap-batch" 2>/dev/null || true

echo ""
echo "============================================"
echo " API Test Summary"
echo "============================================"
echo -e " \033[0;32m$PASS\033[0m passed"
echo -e " \033[0;31m$FAIL\033[0m failed"
echo "============================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
