#!/bin/sh
# Functional test for --meth bisulfite alignment.
#
# Bisulfite-converts the bundled chrM reads (R1 C->T, R2 G->A) to simulate a
# fully-unmethylated directional library, then maps them with --meth. The
# asymmetric scoring matrix must tolerate the conversions, so converted reads
# align near full length with low edit distance. As a control, mapping the same
# converted reads WITHOUT --meth degrades badly (conversions read as
# mismatches), so meth mode must be clearly better.
#
# Run from the repo root after `make`.

set -e

MINIBWA="${MINIBWA:-./minibwa}"
TMP="${TMPDIR:-/tmp}/minibwa-meth-align.$$"
mkdir -p "$TMP"
trap "rm -rf $TMP" EXIT

gzip -dc test/chrM-human.fa.gz > "$TMP/ref.fa"
# Bisulfite-convert sequence lines only (FASTA = alternating header/seq lines).
gzip -dc test/chrM-read_1.fa.gz | awk 'NR%2==0{gsub(/[Cc]/,"T")}1' > "$TMP/bs_1.fa"
gzip -dc test/chrM-read_2.fa.gz | awk 'NR%2==0{gsub(/[Gg]/,"A")}1' > "$TMP/bs_2.fa"

"$MINIBWA" index "$TMP/ref.fa" 2>/dev/null
"$MINIBWA" index --meth "$TMP/ref.fa" 2>/dev/null

"$MINIBWA" map -a --meth "$TMP/ref.fa" "$TMP/bs_1.fa" "$TMP/bs_2.fa" 2>/dev/null > "$TMP/meth.sam"
"$MINIBWA" map -a        "$TMP/ref.fa" "$TMP/bs_1.fa" "$TMP/bs_2.fa" 2>/dev/null > "$TMP/norm.sam"

# Mean NM and primary-alignment count over a SAM file.
# Prints "<count> <mean_nm>".
nm_stats () {
	grep -v '^@' "$1" | awk '{
		flag = $2
		if (int(flag/256)%2 || int(flag/2048)%2) next   # skip secondary/supplementary
		if (int(flag/4)%2) next                          # skip unmapped
		for (i = 12; i <= NF; ++i)
			if ($i ~ /^NM:i:/) { split($i, a, ":"); s += a[3]; ++n }
	} END { printf "%d %.3f", n, (n ? s/n : 999) }'
}

set -- $(nm_stats "$TMP/meth.sam"); METH_N=$1; METH_NM=$2
set -- $(nm_stats "$TMP/norm.sam"); NORM_N=$1; NORM_NM=$2

echo "  --meth : $METH_N primary alns, mean NM $METH_NM"
echo "  normal : $NORM_N primary alns, mean NM $NORM_NM"

# At least 95% of the 2000 reads must map in meth mode.
if [ "$METH_N" -lt 1900 ]; then
	echo "FAIL: --meth mapped only $METH_N/2000 reads" >&2; exit 1
fi
# Converted reads must align cleanly under the asymmetric matrix.
if ! awk -v v="$METH_NM" 'BEGIN{exit !(v < 2.0)}'; then
	echo "FAIL: --meth mean NM ($METH_NM) too high; bisulfite tolerance broken" >&2; exit 1
fi
# ...and clearly beat plain mapping of the same converted reads.
if ! awk -v m="$METH_NM" -v n="$NORM_NM" 'BEGIN{exit !(n > m*3)}'; then
	echo "FAIL: --meth ($METH_NM) not clearly better than normal ($NORM_NM)" >&2; exit 1
fi

echo "PASS: --meth aligns bisulfite-converted reads with low edit distance"
