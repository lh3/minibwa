#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kalloc.h"
#include "kommon.h"
#include "meth.h"

/* Strip leading 'f'/'r' c2t prefix; collapse paired f<x>/r<x> contigs onto a
 * single output name; non-prefixed contigs pass through with direction 0. */
mb_meth_cmap_t *mb_meth_cmap_build(const l2b_t *l2b)
{
	mb_meth_cmap_t *m;
	int i, j;

	if (l2b == 0) return 0;
	m = kom_calloc(mb_meth_cmap_t, 1);
	m->n_internal = (int)l2b->n_ctg;
	if (m->n_internal > 0) {
		m->out_tid      = kom_calloc(int,     m->n_internal);
		m->direction    = kom_calloc(char,    m->n_internal);
		m->output_names = kom_calloc(char *,  m->n_internal);
		m->output_lens  = kom_calloc(int64_t, m->n_internal);
	}
	for (i = 0; i < m->n_internal; ++i) {
		const char *name = l2b->ctg[i].name;
		const char *stripped = (name[0] == 'f' || name[0] == 'r')? name + 1 : name;
		int existing = -1;
		m->direction[i] = (name[0] == 'f' || name[0] == 'r')? name[0] : 0;
		for (j = 0; j < m->n_output; ++j)
			if (strcmp(m->output_names[j], stripped) == 0) { existing = j; break; }
		if (existing >= 0) m->out_tid[i] = existing;
		else {
			int idx = m->n_output;
			m->output_names[idx] = kom_strdup(stripped);
			m->output_lens[idx] = (int64_t)l2b->ctg[i].len;
			m->out_tid[i] = idx;
			++m->n_output;
		}
	}
	return m;
}

void mb_meth_cmap_free(mb_meth_cmap_t *m)
{
	int i;
	if (m == 0) return;
	free(m->out_tid);
	free(m->direction);
	if (m->output_names) {
		for (i = 0; i < m->n_output; ++i) free(m->output_names[i]);
		free(m->output_names);
	}
	free(m->output_lens);
	free(m);
}

void mb_meth_ingest(int n_seq, mb_bseq1_t *seq, int frag_mode, int split_files)
{
	int i, prev_is_r1 = 0;
	for (i = 0; i < n_seq; ++i) {
		mb_bseq1_t *t = &seq[i];
		const char *prior;
		size_t prior_len, cap;
		char *comment, from, from_lo, to;
		const char *yc;
		int is_r2, l = (int)t->l_seq, k, off;

		if (!frag_mode)        is_r2 = 0;
		else if (split_files)  is_r2 = (i & 1);                    /* split-file PE: bseq_read_frag interleaves */
		else if (prev_is_r1 && i > 0 && mb_qname_same(seq[i-1].name, t->name))
		                       is_r2 = 1, prev_is_r1 = 0;
		else                   is_r2 = 0, prev_is_r1 = 1;

		yc = is_r2? "GA" : "CT";
		from = is_r2? 'G' : 'C';
		from_lo = is_r2? 'g' : 'c';
		to = is_r2? 'A' : 'T';

		prior = t->comment;
		prior_len = prior? strlen(prior) : 0;
		cap = (size_t)l + 32 + (prior_len? prior_len + 1 : 0);
		comment = kom_malloc(char, cap);
		off = snprintf(comment, cap, "YS:Z:");
		memcpy(comment + off, t->seq, (size_t)l);
		off += l;
		off += snprintf(comment + off, cap - off, "\tYC:Z:%s", yc);
		if (prior_len)
			snprintf(comment + off, cap - off, "\t%s", prior);
		free(t->comment);
		t->comment = comment;
		for (k = 0; k < l; ++k) {
			char c = t->seq[k];
			if (c == from || c == from_lo) t->seq[k] = to;
		}
	}
}

/* Longest M/=/X run on a CIGAR. */
static int meth_cigar_longest_m(const uint32_t *cigar, int n)
{
	int i, longest = 0;
	for (i = 0; i < n; ++i) {
		int op = cigar[i] & 0xf;
		if (op == 0 || op == 7 || op == 8) {
			int len = (int)(cigar[i] >> 4);
			if (len > longest) longest = len;
		}
	}
	return longest;
}

#define MB_METH_MIN_LONGEST_M_PCT 44 /* matches bwameth.py */

static int meth_qc_one(const mb_opt_t *opt, const mb_meth_cmap_t *cmap,
                       const mb_bseq1_t *t, mb_hit_t *r)
{
	int failed = 0;
	char dir = 0;
	if (r == 0 || r->p == 0) return 0;
	if (cmap && r->tid >= 0 && r->tid < cmap->n_internal)
		dir = cmap->direction[(int)r->tid];
	if (dir == 0) return 0; /* not on a doubled contig */
	if (opt->meth_set_as_failed != 0 && opt->meth_set_as_failed == dir)
		failed = 1;
	if (!opt->meth_no_chim && r->p->n_cigar > 0 && t->l_seq > 0) {
		int lm = meth_cigar_longest_m(r->p->cigar, r->p->n_cigar);
		if (100 * lm < MB_METH_MIN_LONGEST_M_PCT * t->l_seq) {
			failed = 1;
			if (r->mapq > 1) r->mapq = 1;
		}
	}
	if (failed) r->qc_fail = 1, r->proper_pair = 0;
	return failed;
}

void mb_meth_apply_qc(const mb_opt_t *opt, const mb_meth_cmap_t *cmap,
                      const mb_bseq1_t *seq, int32_t n_seg,
                      const int32_t *n_hit, mb_hit_t *const *hit)
{
	int any_fail = 0, s, j;
	if (opt == 0 || hit == 0 || n_hit == 0) return;
	for (s = 0; s < n_seg; ++s)
		for (j = 0; j < n_hit[s]; ++j) {
			mb_hit_t *r = &hit[s][j];
			if (r->parent == r->id && meth_qc_one(opt, cmap, &seq[s], r))
				any_fail = 1;
		}
	if (!any_fail) return;
	for (s = 0; s < n_seg; ++s)
		for (j = 0; j < n_hit[s]; ++j)
			hit[s][j].qc_fail = 1, hit[s][j].proper_pair = 0;
}
