#include <zlib.h>
#include <stdio.h>
#include "mbpriv.h"
#include "ketopt.h"
#include "kommon.h"
#include "kseq.h"
#include "kalloc.h"
KSEQ_INIT(gzFile, gzread);

int main_fastmap(int argc, char *argv[])
{
	mb_idx_t *idx;
	mb_bwt_t *bwt;
	ketopt_t o = KETOPT_INIT;
	gzFile fp;
	kseq_t *ks;
	int c, min_len = 19, min_occ = 1, max_occ = 1, max_size_out = 20;
	int use_sa1 = 0;
	uint64_t *sa, m_a = 0;
	mb_sai_t *a = 0;
	kstring_t out = {0};


	while ((c = ketopt(&o, argc, argv, 1, "l:s:w:1", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 's') min_occ = atoi(o.arg);
		else if (c == 'c') max_occ = atoi(o.arg);
		else if (c == 'w') max_size_out = atoi(o.arg);
		else if (c == '1') use_sa1 = 1;
	}

	if (max_occ < min_occ) max_occ = min_occ;
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa fastmap [options] <idx-prefix> <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT     min seed length [%d]\n", min_len);
		fprintf(stderr, "  -s INT     min interval size [%d]\n", min_occ);
		fprintf(stderr, "  -c INT     max interval size [%d]\n", max_occ);
		fprintf(stderr, "  -w INT     max interval size to output coordinates [%d]\n", max_size_out);
		fprintf(stderr, "  -1         use unbatched sa\n");
		return 1;
	}

	idx = mb_idx_load(argv[o.ind]);
	bwt = idx->bwt;
	kom_assert(bwt, "failed to open the BWT file.");
	fp = strcmp(argv[o.ind+1], "-")? gzopen(argv[o.ind+1], "rb") : gzdopen(0, "rb");
	ks = kseq_init(fp);
	sa = kom_calloc(uint64_t, max_size_out);
	while (kseq_read(ks) >= 0) {
		int64_t x = 0, i, n_a = 0;
		mb_sai_t p;
		out.l = 0;
		for (i = 0; i < ks->seq.l; ++i)
			ks->seq.s[i] = kom_nt4_table[(uint8_t)ks->seq.s[i]];
		if (0) { // temporary code for debugging only; will be removed later
			mb_sai_v u = {0,0,0};
			mb_anchor_v v = {0,0,0};
			mb_anchor_t *a;
			int32_t n_hit;
			uint64_t *w;
			mb_mopt_t mo, *opt = &mo;

			mb_mopt_init(opt);
			mb_preset(opt, "sr");
			mb_seed_intv(0, bwt, ks->seq.l, (uint8_t*)ks->seq.s, 19, 10, &u);
			mb_anchor(0, idx, &u, ks->seq.l, 500, &v);
			free(u.a);
			a = mb_lchain_dp(0, opt->max_gap, opt->max_gap, opt->bw, opt->max_chain_skip, opt->max_chain_iter,
							 opt->min_chain_score, opt->chn_pen_gap, opt->chn_pen_skip, v.n, v.a, &n_hit, &w);
			v.a = 0; v.n = v.m = 0; // ownership transferred to a
			free(a);
			continue;
		}
		kom_sprintf_lite(&out, "SQ\t%s\t%ld\n", ks->name.s, ks->seq.l);
		do {
			x = mb_bwt_smem(bwt, ks->seq.l, (uint8_t*)ks->seq.s, x, min_len, min_occ, max_occ, &p);
			if (p.size > 0) {
				kom_grow(mb_sai_t, a, n_a, m_a);
				a[n_a++] = p;
			}
		} while (x < ks->seq.l);
		for (i = 0; i < n_a; ++i) {
			int64_t len;
			kom_sprintf_lite(&out, "EM\t%ld\t%ld\t%ld", a[i].info>>32, a[i].info&0xffffffffull, a[i].size);
			len = (a[i].info&0xffffffffull) - (a[i].info>>32);
			if (a[i].size <= max_size_out) {
				int64_t j, n_sa = a[i].size;
				if (use_sa1) {
					for (j = 0; j < a[i].size; ++j)
						sa[j] = mb_bwt_sa(bwt, a[i].x[0] + j);
				} else {
					for (j = 0; j < a[i].size; ++j)
						sa[j] = a[i].x[0] + j;
					mb_bwt_sa_batch(0, bwt, a[i].size, sa);
				}
				for (j = 0; j < n_sa; ++j) {
					int rev;
					int64_t cid, cst;
					cid = l2b_intv2cid(idx->l2b, sa[j], sa[j] + len, &cst, &rev);
					if (cid < 0) kom_sprintf_lite(&out, "\t.");
					else kom_sprintf_lite(&out, "\t%s:%c%ld", idx->l2b->ctg[cid].name, "+-"[rev], cst + 1);
				}
			} else kom_sprintf_lite(&out, "\t*");
			kom_sprintf_lite(&out, "\n");
		}
		kom_sprintf_lite(&out, "//\n");
		fputs(out.s, stdout);
	}
	free(sa);
	kseq_destroy(ks);
	gzclose(fp);
	mb_bwt_destroy(bwt);
	return 0;
}
