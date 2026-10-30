#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bam.h"
#include "kommon.h"

/*
 * BGZF is a multi-member gzip stream, so gzread() inflates a BAM file without
 * help; only the binary record format is decoded here. This avoids a dependency
 * on htslib at the cost of single-threaded decompression.
 */

#define BAM_BUF_SIZE  0x20000
#define BAM_MIN_BLOCK 32 // the fixed-length part of a record

struct mb_bam_file_s {
	gzFile fp;
	unsigned char *buf; // decompressed input
	int32_t begin, end;
	int32_t is_eof, is_err;
	char *text;         // plain-text header
	unsigned char *rec; // current record, without the leading block_size
	int32_t l_rec, m_rec;
	int32_t has_rec;    // rec[] holds a record not returned to the caller yet
	int32_t first_flag; // FLAG of the first record; -1 if there is none
};

static inline int32_t bam_le32(const unsigned char *p)
{
	return (int32_t)((uint32_t)p[0] | (uint32_t)p[1]<<8 | (uint32_t)p[2]<<16 | (uint32_t)p[3]<<24);
}

static inline int32_t bam_le16(const unsigned char *p)
{
	return (int32_t)((uint32_t)p[0] | (uint32_t)p[1]<<8);
}

/***********
 * Raw I/O *
 ***********/

static int32_t bam_refill(mb_bam_file_t *b)
{
	int ret;
	if (b->is_eof) return 0;
	b->begin = b->end = 0;
	ret = gzread(b->fp, b->buf, BAM_BUF_SIZE);
	if (ret <= 0) {
		b->is_eof = 1;
		if (ret < 0) b->is_err = 1;
		return 0;
	}
	b->end = ret;
	return ret;
}

// read len bytes; return the number of bytes actually read
static int32_t bam_read(mb_bam_file_t *b, void *dst, int32_t len)
{
	unsigned char *p = (unsigned char*)dst;
	int32_t n = 0;
	while (n < len) {
		int32_t c;
		if (b->begin >= b->end && bam_refill(b) == 0) break;
		c = b->end - b->begin;
		if (c > len - n) c = len - n;
		if (p) memcpy(p + n, b->buf + b->begin, c);
		b->begin += c, n += c;
	}
	return n;
}

// read the next primary record into b->rec; return its length, 0 at EOF and -1 on error
static int32_t bam_next_rec(mb_bam_file_t *b)
{
	unsigned char tmp[4];
	for (;;) {
		int32_t bs, n;
		n = bam_read(b, tmp, 4);
		if (n == 0) return 0;
		if (n != 4) { b->is_err = 1; return -1; }
		bs = bam_le32(tmp);
		if (bs < BAM_MIN_BLOCK) { b->is_err = 1; return -1; }
		if (bs > b->m_rec) { // don't overwrite b->rec: a corrupt block_size may ask for more than we can get
			unsigned char *rec;
			if ((rec = kom_realloc(unsigned char, b->rec, bs)) == 0) { b->is_err = 1; return -1; }
			b->rec = rec, b->m_rec = bs;
		}
		if (bam_read(b, b->rec, bs) != bs) { b->is_err = 1; return -1; }
		if (bam_le16(b->rec + 14) & 0x900) continue; // skip secondary and supplementary records
		return bs;
	}
}

/**********
 * Header *
 **********/

mb_bam_file_t *mb_bam_open(gzFile fp)
{
	mb_bam_file_t *b;
	unsigned char tmp[4];
	int32_t l_text, n_ref, i;

	b = kom_calloc(mb_bam_file_t, 1);
	b->fp = fp;
	b->buf = kom_malloc(unsigned char, BAM_BUF_SIZE);
	b->first_flag = -1;
	if (bam_read(b, tmp, 4) != 4) goto bam_open_fail;
	l_text = bam_le32(tmp);
	if (l_text < 0) goto bam_open_fail;
	b->text = kom_malloc(char, (int64_t)l_text + 1); // 64-bit: l_text may be corrupt
	if (b->text == 0) goto bam_open_fail;
	if (bam_read(b, b->text, l_text) != l_text) goto bam_open_fail;
	b->text[l_text] = 0;
	if (bam_read(b, tmp, 4) != 4) goto bam_open_fail;
	n_ref = bam_le32(tmp);
	if (n_ref < 0) goto bam_open_fail;
	for (i = 0; i < n_ref; ++i) { // skipped; contigs come from the minibwa index
		int32_t l_name;
		if (bam_read(b, tmp, 4) != 4) goto bam_open_fail;
		l_name = bam_le32(tmp);
		if (l_name < 0) goto bam_open_fail;
		if (bam_read(b, 0, l_name + 4) != l_name + 4) goto bam_open_fail; // name and l_ref
	}
	// buffer the first record for mb_bam_is_pe()
	b->l_rec = bam_next_rec(b);
	if (b->l_rec < 0) goto bam_open_fail;
	if (b->l_rec > 0) {
		b->has_rec = 1;
		b->first_flag = bam_le16(b->rec + 14);
	}
	return b;

bam_open_fail:
	if (kom_verbose >= 1)
		fprintf(stderr, "[ERROR] failed to parse the BAM header\n");
	mb_bam_close(b);
	return 0;
}

void mb_bam_close(mb_bam_file_t *b)
{
	if (b == 0) return;
	free(b->rec);
	free(b->text);
	free(b->buf);
	free(b);
}

const char *mb_bam_hdr_text(const mb_bam_file_t *b) { return b->text; }

int mb_bam_is_pe(const mb_bam_file_t *b) { return b->first_flag >= 0 && (b->first_flag & 0x1); }

int mb_bam_eof(const mb_bam_file_t *b) { return b->is_eof && !b->has_rec; }

/******************
 * Auxiliary tags *
 ******************/

/*
 * Aux tags not carried over from a BAM input. A tag is dropped when either:
 *
 *   1. minibwa writes it itself. SAM ("The alignment section: optional fields")
 *      requires that "within each alignment line, no TAG may appear more than
 *      once", so passing the input copy through would emit an invalid record.
 *      minibwa writes NM/AS/ms/md on every aligned record (write_tags() in
 *      format.c), plus MC, MQ, SA, XA, CG and n2, and cs, ds or MD under -b.
 *      sam_write_aux() in format.c enforces this regardless of what is listed
 *      here, so a gap in this group costs a tag, not a malformed record.
 *   2. It describes an alignment somebody else made, and says nothing true
 *      about the one minibwa is about to produce.
 *
 * RG is not listed: it is read-intrinsic, so it is carried over unless -R gives
 * a replacement, which format.c writes instead.
 *
 * Where each group is defined. SAM reserves tags starting with X, Y or Z, and
 * any tag containing a lowercase letter, for local use, so only the first group
 * is formally specified; the rest come from the tool that emits them:
 *
 *   SAMtags        NM MD MC MQ SA CG CC CP H0 H1 H2 OA OC OP AS SM AM
 *   bwa            X0 X1 XN XM XO XG XT XA XS   (bwa.1, "SAM ALIGNMENT TAGS")
 *   bowtie2        XS XN XM XO XG YS YF YT      (bowtie2 MANUAL.markdown)
 *   bismark        XM XR XG XA XB YS            (bismark; XM:Z is the meth call
 *                  string, XR/XG the read and genome conversion state)
 *   HISAT2         ZS                           (bismark, "HISAT2 uses ZS:i:")
 *   bwa-meth       YC = conversion applied, YD = strand   (bwameth.py)
 *   minimap2       cs ms s1 s2 cm tp            (minimap2.1)
 *   minibwa        ds md n2                     (minibwa.1, format.c, cs.c)
 *
 * Names are matched without the type, which is what we want even though the
 * groups above collide: bismark's XM:Z/XG:Z/XA:Z mean something quite different
 * from bwa's XM:i/XG:i/XA:Z, and YS is an int for bowtie2 but a string for
 * bismark. Every one of those readings is alignment-derived, so all of them
 * should go.
 */
static const char *bam_aux_drop[] = {
	"NM", "MD", "AS", "XS", "SA", "MC", "MQ", "CG", "XA", "CC", "CP", "SM", "AM",
	"H0", "H1", "H2", "X0", "X1", "XN", "XM", "XO", "XG", "XT", "OA", "OC", "OP",
	"XR", "XB", "YS", "YF", "YT", "ZS", "YC", "YD",
	"cs", "ds", "ms", "md", "s1", "s2", "cm", "tp", "n2", 0
};

static int bam_aux_is_dropped(const char *tag)
{
	int i;
	for (i = 0; bam_aux_drop[i]; ++i)
		if (tag[0] == bam_aux_drop[i][0] && tag[1] == bam_aux_drop[i][1])
			return 1;
	return 0;
}

static void bam_aux_num(kstring_t *v, const unsigned char *p, int type)
{
	if (type == 'c') kom_sprintf_lite(v, "%d", (int)(int8_t)p[0]);
	else if (type == 'C') kom_sprintf_lite(v, "%d", (int)p[0]);
	else if (type == 's') kom_sprintf_lite(v, "%d", (int)(int16_t)bam_le16(p));
	else if (type == 'S') kom_sprintf_lite(v, "%d", bam_le16(p));
	else if (type == 'i') kom_sprintf_lite(v, "%d", bam_le32(p));
	else if (type == 'I') kom_sprintf_lite(v, "%u", (uint32_t)bam_le32(p));
	else if (type == 'f') {
		char tmp[32];
		union { uint32_t i; float f; } u;
		u.i = (uint32_t)bam_le32(p);
		snprintf(tmp, sizeof(tmp), "%g", u.f);
		kom_sprintf_lite(v, "%s", tmp);
	}
}

static inline int bam_aux_num_size(int type)
{
	if (type == 'c' || type == 'C') return 1;
	if (type == 's' || type == 'S') return 2;
	if (type == 'i' || type == 'I' || type == 'f') return 4;
	return 0;
}

// convert the aux section to tab-led SAM fields; return NULL if nothing is kept
static char *bam_aux2str(const unsigned char *p, const unsigned char *end)
{
	kstring_t s = {0,0,0}, v = {0,0,0};
	while (p + 3 <= end) {
		char tag[3], sam_type = 0;
		int type, keep, sz;
		tag[0] = p[0], tag[1] = p[1], tag[2] = 0;
		type = p[2];
		p += 3;
		keep = !bam_aux_is_dropped(tag);
		v.l = 0;
		if (v.s) v.s[0] = 0;
		if (type == 'A') {
			if (p + 1 > end) break;
			kom_sprintf_lite(&v, "%c", *p);
			p += 1, sam_type = 'A';
		} else if ((sz = bam_aux_num_size(type)) > 0) {
			if (p + sz > end) break;
			bam_aux_num(&v, p, type);
			p += sz, sam_type = type == 'f'? 'f' : 'i';
		} else if (type == 'Z' || type == 'H') {
			const unsigned char *q = p;
			while (q < end && *q) ++q;
			if (q == end) break; // unterminated
			kom_sprintf_lite(&v, "%s", (const char*)p);
			p = q + 1, sam_type = type;
		} else if (type == 'B') {
			int32_t n, i, esz;
			if (p + 5 > end) break;
			esz = bam_aux_num_size(p[0]);
			n = bam_le32(p + 1);
			if (esz == 0 || n < 0 || p + 5 + (int64_t)n * esz > end) break;
			kom_sprintf_lite(&v, "%c", p[0]);
			for (i = 0; i < n; ++i) {
				kom_sprintf_lite(&v, ",");
				bam_aux_num(&v, p + 5 + (int64_t)i * esz, p[0]);
			}
			p += 5 + (int64_t)n * esz, sam_type = 'B';
		} else break; // unknown type; the rest can't be parsed
		if (keep) kom_sprintf_lite(&s, "\t%s:%c:%s", tag, sam_type, v.s? v.s : "");
	}
	free(v.s);
	return s.s;
}

/**********
 * Record *
 **********/

static int bam_rec2bseq(const unsigned char *p, int32_t bs, mb_bseq1_t *s, int with_qual, int with_aux)
{
	static const char nt16_table[] = "=ACMGRSVTWYHKDBN";
	int32_t l_qname = p[8], n_cigar = bam_le16(p + 12), flag = bam_le16(p + 14), l_seq = bam_le32(p + 16);
	const unsigned char *seq, *qual;
	int64_t off;
	int32_t i;

	off = (int64_t)BAM_MIN_BLOCK + l_qname + 4 * (int64_t)n_cigar;
	if (l_qname < 1 || l_seq < 0 || off + ((int64_t)l_seq + 1) / 2 + l_seq > bs) return -1; // 64-bit: l_seq may be corrupt
	seq = p + off, qual = seq + (l_seq + 1) / 2;

	memset(s, 0, sizeof(mb_bseq1_t));
	s->l_seq = l_seq;
	s->flag = flag;
	s->name = kom_malloc(char, l_qname);
	memcpy(s->name, p + BAM_MIN_BLOCK, l_qname); // NUL-terminated in BAM
	s->name[l_qname - 1] = 0;
	if (s->name[0] == 0 && kom_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->seq = kom_malloc(char, l_seq + 1);
	for (i = 0; i < l_seq; ++i)
		s->seq[i] = nt16_table[seq[i>>1] >> ((~i&1)<<2) & 0xf];
	s->seq[l_seq] = 0;
	if (with_qual && l_seq > 0 && qual[0] != 0xff) { // 0xff marks a missing QUAL
		s->qual = kom_malloc(char, l_seq + 1);
		for (i = 0; i < l_seq; ++i)
			s->qual[i] = qual[i] + 33;
		s->qual[l_seq] = 0;
	}
	if (flag & 0x10) { // the reverse strand is stored reverse complemented
		kom_revcomp(l_seq, s->seq);
		if (s->qual) kom_reverse(char, l_seq, s->qual);
	}
	if (with_aux) s->aux = bam_aux2str(qual + l_seq, p + bs);
	return 0;
}

int mb_bam_read1(mb_bam_file_t *b, mb_bseq1_t *s, int with_qual, int with_aux)
{
	if (b->has_rec) b->has_rec = 0;
	else {
		b->l_rec = bam_next_rec(b);
		if (b->l_rec <= 0) return b->l_rec == 0? -1 : -2;
	}
	if (bam_rec2bseq(b->rec, b->l_rec, s, with_qual, with_aux) < 0) { b->is_err = 1; return -2; }
	return s->l_seq;
}
