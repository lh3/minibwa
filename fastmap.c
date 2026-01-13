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
	mb_mopt_t opt;
	ketopt_t o = KETOPT_INIT;
	gzFile fp;
	kseq_t *ks;
	int c, min_len = 19, min_occ = 1, max_occ = 1, max_size_out = 20;
	int max_sub_occ = 10, max_anchor_occ = 500;
	int use_sa1 = 0, test_seed = 0, test_anchor = 0, test_chain = 0;
	const char *preset = 0;
	uint64_t *sa, m_a = 0;
	mb_sai_t *a = 0;
	kstring_t out = {0};
	mb_sai_v u = {0,0,0};
	mb_anchor_v v = {0,0,0};

	mb_mopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "l:s:w:1c:SACx:", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 's') min_occ = atoi(o.arg);
		else if (c == 'c') max_occ = atoi(o.arg);
		else if (c == 'w') max_size_out = atoi(o.arg);
		else if (c == '1') use_sa1 = 1;
		else if (c == 'S') test_seed = 1;
		else if (c == 'A') test_seed = test_anchor = 1;
		else if (c == 'C') test_seed = test_anchor = test_chain = 1;
		else if (c == 'x') preset = o.arg;
	}
	if (preset && mb_set_preset(preset, &opt) < 0) {
		fprintf(stderr, "[ERROR] unknown preset '%s'\n", preset);
		return 1;
	}
	if (max_occ < min_occ) max_occ = min_occ;
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa fastmap [options] <idx-prefix> <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -x STR     preset: sr (short reads) []\n");
		fprintf(stderr, "  -l INT     min seed length [%d]\n", min_len);
		fprintf(stderr, "  -s INT     min interval size [%d]\n", min_occ);
		fprintf(stderr, "  -c INT     max interval size [%d]\n", max_occ);
		fprintf(stderr, "  -w INT     max interval size to output coordinates [%d]\n", max_size_out);
		fprintf(stderr, "  -S         test the seeding algorithm\n");
		fprintf(stderr, "  -A         test the anchoring algorithm\n");
		fprintf(stderr, "  -C         test the chaining algorithm\n");
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
		if (test_anchor || test_seed) {
			mb_seed_intv(0, bwt, ks->seq.l, (uint8_t*)ks->seq.s, min_len, max_sub_occ, &u);
			if (test_anchor) {
				int32_t j;
				mb_anchor(0, idx, &u, ks->seq.l, max_anchor_occ, &v);
				if (test_chain) {
					int n_u;
					uint64_t *cu;
					mb_anchor_t *ca;
					mb_hit_t *hit;
					ca = mb_lchain_dp(opt.max_gap, opt.max_gap, opt.bw, opt.max_chain_skip, opt.max_chain_iter,
									  opt.min_chain_score, opt.chn_pen_gap, opt.chn_pen_skip,
									  v.n, v.a, &n_u, &cu, 0);
					v.a = 0; v.n = v.m = 0; // ownership transferred to ca
					hit = mb_gen_hit(0, 0, ks->seq.l, idx->l2b, n_u, cu, ca);
					mb_set_parent(0, opt.mask_level, opt.mask_len, n_u, hit, opt.sub_diff, 0);
					for (j = 0; j < n_u; ++j) {
						mb_hit_t *h = &hit[j];
						if (h->tid < 0) continue;
						kom_sprintf_lite(&out, "%s\t%ld\t%d\t%d\t%c\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\ttp:A:%c\tcm:i:%d\n",
							ks->name.s, (long)ks->seq.l, h->qs, h->qe, "+-"[h->rev],
							idx->l2b->ctg[h->tid].name, (long)idx->l2b->ctg[h->tid].len,
							(long)h->ts, (long)h->te, h->mlen, h->blen, 255, h->parent == h->id? 'P' : 'S', h->cnt);
					}
					free(ca); free(cu); free(hit);
				} else {
					for (j = 0; j < v.n; ++j) {
						mb_anchor_t *q = &v.a[j];
						int32_t qs, qe;
						if (!(q->sid&1)) {
							qs = q->qpos + 1 - q->len;
							qe = q->qpos + 1;
						} else {
							qs = ks->seq.l - (q->qpos + 1);
							qe = ks->seq.l - (q->qpos + 1 - q->len);
						}
						kom_sprintf_lite(&out, "AR\t%s\t%d\t%d\t%s\t%ld\t%ld\t%c\n",
							ks->name.s, qs, qe, idx->l2b->ctg[q->sid>>1].name, q->tpos + 1 - q->len, q->tpos + 1, "+-"[q->sid&1]);
					}
				}
			} else {
				for (i = 0; i < u.n; ++i) {
					uint32_t st = u.a[i].info >> 32;
					uint32_t en = (uint32_t)u.a[i].info;
					kom_sprintf_lite(&out, "%s\t%u\t%u\t%ld\n", ks->name.s, st, en, (long)u.a[i].size);
				}
			}
		} else {
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
		}
		fputs(out.s, stdout);
	}
	free(u.a); free(v.a);
	free(sa);
	kseq_destroy(ks);
	gzclose(fp);
	mb_bwt_destroy(bwt);
	return 0;
}
