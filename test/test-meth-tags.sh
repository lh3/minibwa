#!/bin/sh
# End-to-end smoke test for --meth SAM tag emission.
#
# Builds a small synthetic reference with two CpG-rich blocks and a PE
# read pair simulating an OT-source fragment with methylated CpGs:
# - R1 reads the top strand directly (BS preserves Cs at methylated CpGs).
# - R2 reads the bottom strand (revcomp displayed in SAM).
#
# Verifies XR/XG (Bismark OT encoding) and XM (Z at every CpG cytosine).
#
# Run from the repo root after `make`.

set -e

MINIBWA="${MINIBWA:-./minibwa}"
TMP="${TMPDIR:-/tmp}/minibwa-meth-tags.$$"
mkdir -p "$TMP"
trap "rm -rf $TMP" EXIT

cat > "$TMP/ref.fa" <<'EOF'
>chr1
AAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCG
EOF

cat > "$TMP/r1.fq" <<'EOF'
@read_meth/1
AAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > "$TMP/r2.fq" <<'EOF'
@read_meth/2
CGCGCGCGTTTTTTTTTTCGCGCGCGTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

"$MINIBWA" index --meth "$TMP/ref.fa" 2>/dev/null
"$MINIBWA" map --meth -a "$TMP/ref.fa" "$TMP/r1.fq" "$TMP/r2.fq" 2>/dev/null > "$TMP/out.sam"

# minibwa emits SAM records in segment order (R1 then R2 per query), so
# the first non-header record is R1 and the second is R2.
R1=$(grep -v '^@' "$TMP/out.sam" | sed -n '1p')
R2=$(grep -v '^@' "$TMP/out.sam" | sed -n '2p')

# OT fragment, R1: XR=CT, XG=CT.
echo "$R1" | grep -q 'XR:Z:CT' || { echo "FAIL: R1 XR != CT" >&2; echo "$R1"; exit 1; }
echo "$R1" | grep -q 'XG:Z:CT' || { echo "FAIL: R1 XG != CT" >&2; echo "$R1"; exit 1; }

# OT fragment, R2: XR=GA, XG=CT (Bismark CTOT-style encoding for R2 of OT).
echo "$R2" | grep -q 'XR:Z:GA' || { echo "FAIL: R2 XR != GA" >&2; echo "$R2"; exit 1; }
echo "$R2" | grep -q 'XG:Z:CT' || { echo "FAIL: R2 XG != CT" >&2; echo "$R2"; exit 1; }

# XM: Z at top-strand C of every CpG (positions 10,12,14,16, 28,30,32,34).
EXPECTED_XM='..........Z.Z.Z.Z...........Z.Z.Z.Z'
for rec in "$R1" "$R2"; do
	XM=$(echo "$rec" | grep -oE 'XM:Z:[^[:space:]]+' | sed 's/XM:Z://')
	# trailing-position behavior is reference-end dependent; check the prefix.
	prefix=$(echo "$XM" | cut -c1-${#EXPECTED_XM})
	if [ "$prefix" != "$EXPECTED_XM" ]; then
		echo "FAIL: XM prefix mismatch" >&2
		echo "  expected: $EXPECTED_XM" >&2
		echo "  got:      $prefix" >&2
		exit 1
	fi
done

# ---------------------------------------------------------------------------
# Reverse-strand soft-clip framing.
#
# A read whose 3' tail runs off the reference maps to the reverse strand with
# a *leading* soft-clip in SAM (SAM SEQ is the reverse complement of the
# read). XM is written in SAM SEQ orientation, so the soft-clipped columns
# must be '.'. A regression that seeds the XM walk at r->qs (correct only for
# forward reads) instead of qlen - r->qe for reverse reads shifts the calls
# into the clipped region.
TOP='GATTACAGGCATGCCACCGTGCCCGGCTAATTTACGTGGAGCTCTAGCT'
BOTTOM=$(printf '%s' "$TOP" | rev | tr ACGTacgt TGCAtgca)
TAIL='TTAAGGCCTTAA'
READ="${BOTTOM}${TAIL}"
QUAL=$(printf '%s' "$READ" | tr '[:print:]' 'I')

cat > "$TMP/refrc.fa" <<EOF
>chr1
$TOP
EOF

cat > "$TMP/rc.fq" <<EOF
@rev_clip
$READ
+
$QUAL
EOF

"$MINIBWA" index --meth "$TMP/refrc.fa" 2>/dev/null
"$MINIBWA" map --meth -a "$TMP/refrc.fa" "$TMP/rc.fq" 2>/dev/null > "$TMP/rc.sam"

REC=$(grep -v '^@' "$TMP/rc.sam" | sed -n '1p')
FLAG=$(echo "$REC" | cut -f2)
CIGAR=$(echo "$REC" | cut -f6)
XM=$(echo "$REC" | grep -oE 'XM:Z:[^[:space:]]+' | sed 's/XM:Z://')

# Must map to the reverse strand (FLAG 0x10) with a leading soft-clip.
if [ $(( FLAG & 16 )) -eq 0 ]; then
	echo "FAIL: reverse-clip read not reverse-mapped (FLAG=$FLAG)" >&2
	echo "$REC"; exit 1
fi
case "$CIGAR" in
	[0-9]*S[0-9]*M*) : ;;
	*) echo "FAIL: reverse-clip read has no leading soft-clip (CIGAR=$CIGAR)" >&2
	   echo "$REC"; exit 1 ;;
esac

# Every XM column inside the leading soft-clip must be '.'.
CLIP=${CIGAR%%S*}
CLIPPED=$(printf '%s' "$XM" | cut -c1-"$CLIP")
case "$CLIPPED" in
	*[!.]*) echo "FAIL: XM has a methylation call inside the soft-clip" >&2
	        echo "  clip=$CLIP CIGAR=$CIGAR" >&2
	        echo "  XM=$XM" >&2
	        exit 1 ;;
esac

# The read is the genomic reverse complement with no bisulfite conversion, so
# every cytosine is "methylated" and all XM calls must be UPPERCASE. A walk
# seeded at r->qs (correct only for forward reads) reads shifted read bases for
# reverse reads and emits spurious lowercase (unmethylated) / 'u' calls.
case "$XM" in
	*[zxhu]*) echo "FAIL: reverse-clip XM has lowercase calls (mis-framed reverse walk)" >&2
	          echo "  CIGAR=$CIGAR XM=$XM" >&2
	          exit 1 ;;
esac

echo "PASS: --meth XR/XG/XM tags emit correctly for a methylated OT fragment"
echo "PASS: --meth XM is framed correctly for reverse-strand soft-clipped reads"
