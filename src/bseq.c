#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
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
	kseq_t *ks;
	mb_bseq1_t s;
};

mb_bseq_file_t *mb_bseq_open(const char *fn)
{
	mb_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mb_bseq_file_t*)calloc(1, sizeof(mb_bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void mb_bseq_close(mb_bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
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
	s->l_seq = ks->seq.l;
}

mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int min_cnt, int64_t max_chunk_size, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mb_bseq1_t, a, 256);
		kv_push(mb_bseq1_t, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mb_bseq1_t));
	}
	if (max_chunk_size < chunk_size)
		max_chunk_size = chunk_size;
	while ((ret = kseq_read(ks)) >= 0) {
		int32_t to_stop = 0;
		mb_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		kv_pushp(mb_bseq1_t, a, &s);
		kseq2bseq(ks, s, with_qual, with_comment);
		size += ks->seq.l;
		if (chunk_size <= 0 || max_chunk_size <= 0) to_stop = 1;
		else if (size >= max_chunk_size) to_stop = 1;
		else if (size >= chunk_size && a.n >= min_cnt) to_stop = 1;
		if (to_stop) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = kseq_read(ks)) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual, with_comment);
					size += ks->seq.l;
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
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}

mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (kseq_read(fp[i]->ks) >= 0)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				fprintf(stderr, "[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			break; // some file reaches the end
		}
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mb_bseq1_t *s;
			kv_pushp(mb_bseq1_t, a, &s);
			kseq2bseq(fp[i]->ks, s, with_qual, with_comment);
			size += s->l_seq;
		}
		if (size >= chunk_size) break;
	}
	*n_ = a.n;
	return a.a;
}

int mb_bseq_eof(mb_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
