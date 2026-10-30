#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
#include "bam.h"
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

#define kvec_t(type) struct { size_t n, m; type *a; }

#define kv_resize(type, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_push(type, v, x) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		(v).a[(v).n++] = (x); \
	} while (0)

#define kv_pushp(type, v, p) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		*(p) = &(v).a[(v).n++]; \
	} while (0)

#define CHECK_PAIR_THRES 1000000

struct mb_bseq_file_s {
	gzFile fp;
	kseq_t *ks;       // for FASTA/FASTQ input
	mb_bam_file_t *b; // for BAM input
	mb_bseq1_t s;
};

mb_bseq_file_t *mb_bseq_open(const char *fn)
{
	mb_bseq_file_t *fp;
	unsigned char magic[4];
	int n_magic;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	gzbuffer(f, 0x20000);
	fp = (mb_bseq_file_t*)calloc(1, sizeof(mb_bseq_file_t));
	fp->fp = f;
	n_magic = gzread(f, magic, 4); // BGZF is a gzip stream, so this reaches the BAM magic
	if (n_magic == 4 && memcmp(magic, "BAM\1", 4) == 0) {
		fp->b = mb_bam_open(f);
		if (fp->b == 0) { gzclose(f); free(fp); return 0; }
	} else {
		if (n_magic == 4 && memcmp(magic, "@HD\t", 4) == 0 && kom_verbose >= 2) // only warn: a read named "HD" looks the same
			fprintf(stderr, "[WARNING]\033[1;31m '%s' looks like SAM. Convert it to BAM first.\033[0m\n", fn? fn : "-");
		fp->ks = kseq_init(fp->fp);
		if (n_magic > 0) { // push the peeked bytes back to the kstream buffer
			memcpy(fp->ks->f->buf, magic, n_magic);
			fp->ks->f->begin = 0, fp->ks->f->end = n_magic;
		}
	}
	return fp;
}

void mb_bseq_close(mb_bseq_file_t *fp)
{
	mb_bseq_free1(&fp->s);
	if (fp->b) mb_bam_close(fp->b);
	if (fp->ks) kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

void mb_bseq_free1(mb_bseq1_t *s)
{
	free(s->name); free(s->seq); free(s->qual); free(s->comment); free(s->aux);
	memset(s, 0, sizeof(mb_bseq1_t));
}

const char *mb_bseq_hdr_text(const mb_bseq_file_t *fp)
{
	return fp->b? mb_bam_hdr_text(fp->b) : 0;
}

int mb_bseq_is_pe(const mb_bseq_file_t *fp)
{
	return fp->b? mb_bam_is_pe(fp->b) : 0;
}

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mb_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	if (ks->name.l == 0)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->aux = 0, s->flag = 0;
	s->l_seq = ks->seq.l;
}

// read one sequence; return its length, -1 at EOF and <-1 on a parse error
static int bseq_read1(mb_bseq_file_t *fp, mb_bseq1_t *s, int with_qual, int with_comment, int with_aux)
{
	int ret;
	if (fp->b) return mb_bam_read1(fp->b, s, with_qual, with_aux);
	if ((ret = kseq_read(fp->ks)) < 0) return ret;
	kseq2bseq(fp->ks, s, with_qual, with_comment);
	return s->l_seq;
}

mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int with_aux, int frag_mode, int min_cnt, int64_t max_chunk_size, int *n_)
{
	int64_t size = 0;
	int ret = 0;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	mb_bseq1_t t;
	memset(&t, 0, sizeof(mb_bseq1_t));
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mb_bseq1_t, a, 256);
		kv_push(mb_bseq1_t, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mb_bseq1_t));
	}
	if (max_chunk_size < chunk_size)
		max_chunk_size = chunk_size;
	while ((ret = bseq_read1(fp, &t, with_qual, with_comment, with_aux)) >= 0) {
		int32_t to_stop = 0;
		assert(t.l_seq <= INT32_MAX);
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		size += t.l_seq;
		kv_push(mb_bseq1_t, a, t);
		memset(&t, 0, sizeof(mb_bseq1_t));
		if (chunk_size <= 0 || max_chunk_size <= 0) to_stop = 1;
		else if (size >= max_chunk_size) to_stop = 1;
		else if (size >= chunk_size && a.n >= min_cnt) to_stop = 1;
		if (to_stop) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = bseq_read1(fp, &fp->s, with_qual, with_comment, with_aux)) >= 0) {
					size += fp->s.l_seq;
					if (mb_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mb_bseq1_t, a, fp->s);
						memset(&fp->s, 0, sizeof(mb_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	if (ret < -1) {
		const char *what = fp->b? "BAM" : "FASTA/FASTQ";
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the %s record next to '%s'. Continue anyway.\033[0m\n", what, a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first %s record. Continue anyway.\033[0m\n", what);
	}
	if (a.n == 0) { free(a.a); a.a = 0; } // NULL tells the caller the input has ended
	*n_ = a.n;
	return a.a;
}

mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int with_aux, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	mb_bseq1_t *reads;
	*n_ = 0;
	if (n_fp < 1) return 0;
	reads = (mb_bseq1_t*)calloc(n_fp, sizeof(mb_bseq1_t));
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (bseq_read1(fp[i], &reads[i], with_qual, with_comment, with_aux) >= 0)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				fprintf(stderr, "[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			for (i = 0; i < n_fp; ++i) mb_bseq_free1(&reads[i]);
			break; // some file reaches the end
		}
		for (i = 1; i < n_fp; ++i) // mates are paired by position, so the order must agree
			if (!mb_qname_same(reads[0].name, reads[i].name)) {
				fprintf(stderr, "[ERROR] query files are out of sync: read %ld of file 1 is '%s' but file %d has '%s'.\n",
						(long)(a.n / n_fp + 1), reads[0].name, i + 1, reads[i].name);
				fprintf(stderr, "        mates are paired by position, so the files must list reads in the same order\n");
				fprintf(stderr, "        sort them the same way, or interleave the reads into one input\n");
				exit(1);
			}
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		for (i = 0; i < n_fp; ++i) {
			kv_push(mb_bseq1_t, a, reads[i]);
			size += reads[i].l_seq;
			memset(&reads[i], 0, sizeof(mb_bseq1_t));
		}
		if (size >= chunk_size) break;
	}
	free(reads);
	*n_ = a.n;
	return a.a;
}

int mb_bseq_eof(mb_bseq_file_t *fp)
{
	if (fp->b) return (mb_bam_eof(fp->b) && fp->s.seq == 0);
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
