#include <zlib.h>
#include <stdio.h>
#include "mbpriv.h"
#include "ketopt.h"
#include "kommon.h"
#include "kseq.h"
#include "kalloc.h"
KSEQ_INIT(gzFile, gzread);

typedef struct {
	char *name;
	int32_t l_seq;
	uint8_t *seq;
	mb_sai_v v;
} batch_seq1_t;

static void batch_smem(const mb_idx_t *idx, int32_t n, batch_seq1_t *t, int32_t min_len, int32_t min_occ)
{
	int32_t i;
	mb_smem_entry_t *s;
	s = kom_calloc(mb_smem_entry_t, n);
	for (i = 0; i < n; ++i) {
		batch_seq1_t *ti = &t[i];
		mb_smem_entry_t *si = &s[i];
		si->min_len = min_len, si->min_occ = min_occ;
		si->st = 0, si->en = ti->l_seq;
		si->q = ti->seq;
		si->v = &ti->v;
		si->v->n = 0;
	}
	mb_bwt_smem_batch(0, idx->bwt, n, s);
	free(s);
}

static void write_intv(const mb_idx_t *idx, const mb_sai_t *p, int32_t max_size_out, uint64_t *sa, kstring_t *out)
{
	int64_t len;
	kom_sprintf_lite(out, "EM\t%ld\t%ld\t%ld", p->info>>32, p->info&0xffffffffull, p->size);
	len = (p->info&0xffffffffull) - (p->info>>32);
	if (p->size <= max_size_out) {
		int64_t j, n_sa = p->size;
		for (j = 0; j < p->size; ++j)
			sa[j] = p->x[0] + j;
		mb_bwt_sa_batch(0, idx->bwt, p->size, sa);
		for (j = 0; j < n_sa; ++j) {
			int rev;
			int64_t cid, cst;
			cid = l2b_intv2cid(idx->l2b, sa[j], sa[j] + len, &cst, &rev);
			if (cid < 0) kom_sprintf_lite(out, "\t.");
			else kom_sprintf_lite(out, "\t%s:%c%ld", idx->l2b->ctg[cid].name, "+-"[rev], cst + 1);
		}
	} else kom_sprintf_lite(out, "\t*");
	kom_sprintf_lite(out, "\n");
}

static void process_batch(const mb_idx_t *idx, int32_t n, batch_seq1_t *t, int32_t min_len, int32_t min_occ, int32_t max_size_out, uint64_t *sa, kstring_t *out)
{
	int32_t i, j;
	batch_smem(idx, n, t, min_len, min_occ);
	for (i = 0; i < n; ++i) {
		batch_seq1_t *ti = &t[i];
		out->l = 0;
		kom_sprintf_lite(out, "SQ\t%s\t%d\n", ti->name, ti->l_seq);
		for (j = 0; j < ti->v.n; ++j)
			write_intv(idx, &ti->v.a[j], max_size_out, sa, out);
		free(ti->name);
		free(ti->seq);
		kom_sprintf_lite(out, "//\n");
		fputs(out->s, stdout);
	}
}

static int usage_fastmap(FILE *fp, int min_len, int min_occ, int max_size_out, int max_seq)
{
	fprintf(fp, "Usage: minibwa fastmap [options] <idx-prefix> <in.fq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -l INT     min seed length [%d]\n", min_len);
	fprintf(fp, "  -s INT     min interval size [%d]\n", min_occ);
	fprintf(fp, "  -w INT     max interval size to output coordinates [%d]\n", max_size_out);
	fprintf(fp, "  -b INT     batch size [%d]\n", max_seq);
	fprintf(fp, "  --help     print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_fastmap(int argc, char *argv[])
{
	mb_idx_t *idx;
	mb_bwt_t *bwt;
	ketopt_t o = KETOPT_INIT;
	gzFile fp;
	kseq_t *ks;
	int c, min_len = 19, min_occ = 1, max_size_out = 20, max_seq = 1;
	uint64_t *sa, m_a = 0;
	mb_sai_t *a = 0;
	kstring_t out = {0};
	int32_t n_seq = 0;
	batch_seq1_t *seq = 0;
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};

	while ((c = ketopt(&o, argc, argv, 1, "l:s:w:b:", long_opts)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 's') min_occ = atoi(o.arg);
		else if (c == 'w') max_size_out = atoi(o.arg);
		else if (c == 'b') max_seq = atoi(o.arg);
		else if (c == 901) return usage_fastmap(stdout, min_len, min_occ, max_size_out, max_seq);
	}

	if (argc - o.ind < 2) return usage_fastmap(stderr, min_len, min_occ, max_size_out, max_seq);

	idx = mb_idx_load(argv[o.ind], 0);
	bwt = idx->bwt;
	kom_assert(bwt, "failed to open the BWT file.");
	fp = strcmp(argv[o.ind+1], "-")? gzopen(argv[o.ind+1], "rb") : gzdopen(0, "rb");
	ks = kseq_init(fp);
	sa = kom_calloc(uint64_t, max_size_out);
	if (max_seq > 1) seq = kom_calloc(batch_seq1_t, max_seq);
	while (kseq_read(ks) >= 0) {
		int64_t i;
		for (i = 0; i < ks->seq.l; ++i)
			ks->seq.s[i] = kom_nt4_table[(uint8_t)ks->seq.s[i]];
		if (max_seq > 1) {
			batch_seq1_t *s;
			if (n_seq == max_seq) {
				process_batch(idx, n_seq, seq, min_len, min_occ, max_size_out, sa, &out);
				n_seq = 0;
			}
			s = &seq[n_seq++];
			s->l_seq = ks->seq.l;
			s->name = kom_strdup(ks->name.s);
			s->seq = kom_malloc(uint8_t, s->l_seq);
			memcpy(s->seq, ks->seq.s, s->l_seq);
		} else {
			int64_t x = 0, n_a = 0;
			mb_sai_t p;
			out.l = 0;
			kom_sprintf_lite(&out, "SQ\t%s\t%ld\n", ks->name.s, ks->seq.l);
			do {
				x = mb_bwt_smem(bwt, ks->seq.l, (uint8_t*)ks->seq.s, x, min_len, min_occ, &p);
				if (p.size > 0) {
					kom_grow(mb_sai_t, a, n_a, m_a);
					a[n_a++] = p;
				}
			} while (x < ks->seq.l);
			for (i = 0; i < n_a; ++i)
				write_intv(idx, &a[i], max_size_out, sa, &out);
			kom_sprintf_lite(&out, "//\n");
			fputs(out.s, stdout);
		}
	}
	if (max_seq > 1)
		process_batch(idx, n_seq, seq, min_len, min_occ, max_size_out, sa, &out); // TODO: then free seq[]
	free(sa);
	kseq_destroy(ks);
	gzclose(fp);
	mb_bwt_destroy(bwt);
	return 0;
}
