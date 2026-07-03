#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksw2.h"

static inline int mb_insert_dir(const mb_hit_t *h0, const mb_hit_t *h1, int64_t *dist)
{
	int64_t p0, p1;
	p0 = h0->rev? h0->te : h0->ts;
	p1 = h1->rev? h1->te : h1->ts;
	*dist = p0 > p1? p0 - p1 : p1 - p0;
	return ((int32_t)h0->rev << 1 | (int32_t)h1->rev) ^ (p0 < p1? 0 : 3);
}

static inline double mb_pair_score(const mb_hit_t *h0, const mb_hit_t *h1, const mb_pestat_t pes[4], int32_t match_sc)
{
	const double MB_SQRT1_2 = 0.707106781186547524401;
	int32_t dir;
	int64_t dist;
	dir = mb_insert_dir(h0, h1, &dist);
	if (!pes[dir].failed && dist >= pes[dir].lo && dist <= pes[dir].hi) {
		double ns = (dist - pes[dir].avg) / pes[dir].std; // normalized score
		return h0->p->dp_max + h1->p->dp_max + .721 * log(2. * erfc(fabs(ns) * MB_SQRT1_2)) * match_sc; // .721 = 1/log(4)
	}
	return -1.0;
}

static const mb_hit_t *mb_select_unique_se(int32_t n_hit, const mb_hit_t *hit)
{
	int32_t j, n_pri = 0, mapq = 0, k = -1;
	for (j = 0; j < n_hit; ++j)
		if (hit[j].id == hit[j].parent) 
			++n_pri, mapq = hit[j].mapq, k = j;
	return n_pri == 1 && mapq >= 10? &hit[k] : 0;
}

static int32_t mb_hit_sum_score(void *km, int32_t n_hit, const mb_hit_t *hit)
{
	static const int32_t max_cnt_heap = 8;
	int32_t i, k, n_pri, sc, qe;
	uint64_t *a, aa[max_cnt_heap];
	for (i = 0, n_pri = sc = 0; i < n_hit; ++i) // precalculate size
		if (hit[i].id == hit[i].parent && hit[i].p)
			++n_pri, sc = hit[i].p->dp_max; // NB: requiring base alignment
	if (n_pri == 0) return 0;
	if (n_pri == 1) return sc;
	a = n_pri < max_cnt_heap? aa : Kmalloc(km, uint64_t, n_pri);
	for (i = k = 0; i < n_hit; ++i)
		if (hit[i].id == hit[i].parent && hit[i].p)
			a[k++] = (uint64_t)hit[i].qs << 32 | i;
	radix_sort_mb64(a, a + k);
	for (i = 0, sc = qe = 0; i < k; ++i) {
		const mb_hit_t *h = &hit[(uint32_t)a[i]];
		if (h->qe <= qe) continue;
		if (h->qs < qe)
			sc += (int32_t)((double)(h->qe - qe) / (h->qe - h->qs) * h->p->dp_max + .499);
		else sc += h->p->dp_max;
		qe = h->qe;
	}
	if (a != aa) kfree(km, a);
	return sc;
}

void mb_pestat(void *km, const mb_opt_t *opt, int32_t n_frag, const int32_t *seg_off, const int32_t *seg_cnt, const int32_t *n_hit, mb_hit_t *const *hit, mb_pestat_t pes[4])
{
	const int MIN_DIR_CNT = 20;
	const double MIN_DIR_RATIO = 0.05, OUTLIER_BOUND = 2.0, MAPPING_BOUND = 3.0, MAX_STDDEV = 4.0;
	int32_t i, d, max;
	struct { int32_t n, m; uint64_t *a; } is[4], *q;
	memset(is, 0, sizeof(is[0]) * 4);
	memset(pes, 0, sizeof(pes[0]) * 4);
	for (i = 0; i < n_frag; ++i) {
		const mb_hit_t *r[2];
		int32_t off, dir;
		int64_t dist;
		if (seg_cnt[i] != 2) continue;
		off = seg_off[i];
		r[0] = mb_select_unique_se(n_hit[off + 0], hit[off + 0]);
		r[1] = mb_select_unique_se(n_hit[off + 1], hit[off + 1]);
		if (r[0] == 0 || r[1] == 0) continue;
		if (r[0]->tid != r[1]->tid) continue; // not on the same contig
		dir = mb_insert_dir(r[0], r[1], &dist);
		if (dist < opt->max_pe_ins) {
			if (is[dir].n == is[dir].m)
				Kgrow(km, uint64_t, is[dir].a, is[dir].n, is[dir].m);
			is[dir].a[is[dir].n++] = dist;
		}
	}
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%d, %d, %d, %d)\n", __func__, is[0].n, is[1].n, is[2].n, is[3].n);
	for (d = 0, max = 0; d < 4; ++d)
		max = max > is[d].n? max : is[d].n;
	for (d = 0; d < 4; ++d) {
		mb_pestat_t *r = &pes[d];
		q = &is[d];
		int p25, p50, p75, x;
		if (q->n < MIN_DIR_CNT || q->n < max * MIN_DIR_RATIO) {
			r->failed = 1;
			kfree(km, q->a);
			continue;
		}
		radix_sort_mb64(q->a, q->a + q->n);
		p25 = q->a[(int)(.25 * q->n + .499)];
		p50 = q->a[(int)(.50 * q->n + .499)];
		p75 = q->a[(int)(.75 * q->n + .499)];
		r->lo  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->lo < 1) r->lo = 1;
		r->hi = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->lo && q->a[i] <= r->hi)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->lo && q->a[i] <= r->hi)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] (25, 50, 75) percentile: (%d, %d, %d); mean and std.dev: (%.2f, %.2f)\n",
				__func__, "FR"[d>>1&1], "FR"[d&1], p25, p50, p75, r->avg, r->std);
		r->lo  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->hi = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->lo  > r->avg - MAX_STDDEV * r->std) r->lo  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->hi < r->avg + MAX_STDDEV * r->std) r->hi = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->lo < 1) r->lo = 1;
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] low and high boundaries for proper pairs: (%d, %d)\n", __func__, "FR"[d>>1&1], "FR"[d&1], r->lo, r->hi);
		kfree(km, q->a);
	}
}

typedef struct {
	int32_t score, sub_sc, n_sub, n_pp;
	int32_t i[2];
} mb_pairaux_t;

static void mb_pair_hits(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], mb_pairaux_t *ret)
{
	int32_t r, i, k, n_pa, y[4], n_pp = 0, m_pp = 0;
	mb128_t *pa, *pp = 0; // pp: proper pairs

	ret->i[0] = ret->i[1] = ret->score = ret->sub_sc = -1, ret->n_sub = ret->n_pp = 0;
	if (n_hit[0] == 0 || n_hit[1] == 0) return;
	pa = Kcalloc(km, mb128_t, n_hit[0] + n_hit[1]);
	for (r = n_pa = 0; r < 2; ++r) {
		for (i = 0; i < n_hit[r]; ++i) {
			mb128_t *p = &pa[n_pa++];
			mb_hit_t *h = &hit[r][i];
			h->proper_pair = 0;
			p->x = l2b->ctg[h->tid].off + (h->rev? h->te : h->ts);
			p->y = (uint64_t)i << 2 | (uint64_t)h->rev << 1 | r;
		}
	}
	radix_sort_mb128x(pa, pa + n_pa);

	y[0] = y[1] = y[2] = y[3] = -1;
	for (i = 0; i < n_pa; ++i) {
		mb128_t *pi = &pa[i];
		mb_hit_t *hi = &hit[pi->y&1][pi->y>>2];
		for (r = 0; r < 2; ++r) {
			int which, dir = r << 1 | (pi->y>>1&1);
			if (pes[dir].failed) continue; // invalid orientation
			which = r << 1 | ((pi->y&1) ^ 1);
			if (y[which] < 0) continue; // no previous hit
			for (k = y[which]; k >= 0; --k) {
				mb128_t *pk = &pa[k], *q;
				mb_hit_t *hk = &hit[pk->y&1][pk->y>>2];
				int64_t dist;
				double s;
				if ((pk->y&3) != which) continue;
				if (hi->tid != hk->tid) break;
				dist = pi->x - pk->x;
				if (dist > pes[dir].hi) break;
				if (dist < pes[dir].lo) continue;
				hk->proper_pair = hi->proper_pair = 1; // paired
				s = mb_pair_score(hk, hi, pes, opt->a);
				if (s < 0.) s = 0.;
				if (n_pp == m_pp) Kgrow(km, mb128_t, pp, n_pp, m_pp);
				q = &pp[n_pp++];
				q->y = (pk->y&1) == 0? (uint64_t)(pk->y>>2) << 32 | (pi->y>>2) : (uint64_t)(pi->y>>2) << 32 | (pk->y>>2); // upper bits: index to read[0]
				q->x = (uint64_t)(s + .499) << 32 | (hk->hash ^ hi->hash);
			}
		}
		y[pi->y&3] = i;
	}

	ret->n_pp = n_pp;
	if (n_pp > 0) {
		uint64_t max = 0, max2 = 0;
		int32_t tmp = opt->a + opt->b > opt->q + opt->e? opt->a + opt->b : opt->q + opt->e;
		for (i = 0; i < n_pp; ++i) { // find max and max2
			mb128_t *q = &pp[i];
			if (q->x > max) max2 = max, max = q->x, ret->i[0] = q->y>>32, ret->i[1] = (uint32_t)q->y;
			else if (q->x > max2) max2 = q->x;
		}
		assert(ret->i[0] < n_hit[0] && ret->i[1] < n_hit[1]);
		ret->score = max>>32, ret->sub_sc = max2>>32;
		for (i = 0; i < n_pp; ++i)
			if (pp[i].x>>32 >= ret->score - tmp)
				ret->n_sub++;
	}
	kfree(km, pp);
	kfree(km, pa);
}

static int32_t mb_ungap(void *km, int32_t qlen, const uint8_t *qseq, int32_t tlen, const uint8_t *tseq, int32_t kmer, l2b_meth_t mt, int32_t *max_i, int32_t *n_good, int32_t *n_kmer)
{ // a linear algorithm to find ungapped alignment
	static uint8_t c2t[4] = { 0, 3, 2, 3 };
	static uint8_t g2a[4] = { 0, 1, 0, 3 };
	int32_t i, l, cap = 1 << kmer*2, mask = cap - 1, max, *a;
	uint16_t *h, x;
	*max_i = -1, *n_good = *n_kmer = 0;
	if (qlen >= UINT16_MAX) return 0;
	a = Kcalloc(km, int32_t, tlen);
	h = Kcalloc(km, uint16_t, cap);
	for (i = l = 0, x = 0; i < qlen; ++i) {
		if (qseq[i] < 4) {
			uint8_t c = mt == L2B_METH_C2T? c2t[qseq[i]] : mt == L2B_METH_G2A? g2a[qseq[i]] : qseq[i];
			x = (x << 2 | c) & mask;
			if (++l >= kmer) {
				if (h[x] == 0) ++*n_kmer; // n_kmer is the number of distinct k-mers on the query sequence
				h[x] = i;
			}
		} else x = 0, l = 0;
	}
	for (i = l = 0, x = 0; i < tlen; ++i) {
		if (tseq[i] < 4) {
			uint8_t c = mt == L2B_METH_C2T? c2t[tseq[i]] : mt == L2B_METH_G2A? g2a[tseq[i]] : tseq[i];
			x = (x << 2 | c) & mask;
			if (++l >= kmer) {
				if (h[x] > 0 && i >= h[x])
					++a[i - h[x]]; // inspired by the Hough transform for line finding
			}
		} else x = 0, l = 0;
	}
	for (i = 0, max = 0; i < tlen; ++i)
		if (max < a[i])
			max = a[i], *max_i = i;
	for (i = 0, *n_good = 0; i < tlen; ++i)
		if (a[i] > max>>1) ++*n_good;
	kfree(km, h);
	kfree(km, a);
	return max;
}

static void mb_matesw_align(void *km, const mb_opt_t *opt, int32_t qlen, uint8_t *qseq, int32_t tlen, uint8_t *tseq, mb_hit_t *h, int32_t min_sc, l2b_meth_t mt, ksw_extz_t *ez)
{
	int8_t mat[25];
	int32_t max_sc = qlen < tlen? qlen : tlen;
	int32_t b_mm = (opt->b + opt->a - 1) / opt->a;
	int32_t b_ts = (opt->b_ts + opt->a - 1) / opt->a;
	int32_t b_ambi = (opt->b_ambi + opt->a - 1) / opt->a;
	int32_t gapo = (opt->q + opt->a - 1) / opt->a;
	int32_t gape = (opt->e + opt->a - 1) / opt->a;
	int32_t sz, xtra;
	void *qp;
	ksw_llrst_t rst;

	memset(h, 0, sizeof(*h));
	if (max_sc >= 32767) return;
	ksw_gen_nt4_mat(mat, 1, b_mm, b_ts, b_ambi, (int)mt);
	sz = max_sc < 255 - b_mm? 1 : 2;
	qp = ksw_ll_qinit(km, sz, qlen, qseq, 5, mat);
	xtra = KSW_LL_SUBO | opt->min_len;
	if (sz == 1)
		rst = ksw_ll_u8_core(qp, tlen, tseq, gapo, gape, xtra);
	else
		rst = ksw_ll_i16_core(qp, tlen, tseq, gapo, gape, xtra);
	if (kom_dbg_flag & MB_DBG_ALN_PE) {
		int i;
		fprintf(stderr, "===> qlen=%d; tlen=%d; score=%d; qe=%d; te=%d <===\n", qlen, tlen, rst.score, rst.qe + 1, rst.te + 1);
		for (i = 0; i < qlen; ++i) fputc("ACGTN"[qseq[i]], stderr);
		fputc('\n', stderr);
		for (i = 0; i < tlen; ++i) fputc("ACGTN"[tseq[i]], stderr);
		fputc('\n', stderr);
	}
	kfree(km, qp);
	if (rst.score >= opt->min_dp_max && rst.score >= min_sc) { // min_sc is already divided by opt->a
		int32_t te = rst.te + 1, qe = rst.qe + 1, ksw_flag = KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR;
		mb_seq_rev(qe, qseq);
		mb_seq_rev(te, tseq);
		ksw_gen_nt4_mat(mat, opt->a, opt->b, opt->b_ts, opt->b_ambi, (int)mt);
		if (mt != L2B_METH_NONE) ksw_flag |= KSW_EZ_GENERIC_SC;
		ksw_extz2_sse(km, qe, qseq, te, tseq, 5, mat, opt->q, opt->e, opt->bw, opt->zdrop, opt->end_bonus, ksw_flag, ez);
		mb_seq_rev(qe, qseq);
		mb_seq_rev(te, tseq);
		if (ez->n_cigar > 0 && ez->max >= opt->min_dp_max * opt->a) {
			mb_append_cigar(h, ez->n_cigar, ez->cigar);
			h->rescued = 1;
			h->qe = qe, h->te = te;
			h->ts = te - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
			h->qs = qe - (ez->reach_end? qe : ez->max_q + 1);
			h->p->dp_max = h->p->dp_score = ez->max;
			h->p->dp_max2 = (int32_t)((double)ez->max / rst.score * rst.score2 + .499);
			if (h->p->dp_max2 < 0) h->p->dp_max2 = 0;
			h->score = h->score0 = rst.score;
			h->subsc = rst.score2;
			h->cnt = 0, h->as = -1;
			h->parent = MB_PARENT_UNSET;
			if (kom_dbg_flag & MB_DBG_ALN_PE) {
				int i;
				fprintf(stderr, "max=%d; ts=%ld; qs=%d; reach_end=%d; cigar=", ez->max, (long)h->ts, h->qs, ez->reach_end);
				for (i = 0; i < ez->n_cigar; ++i) fprintf(stderr, "%d%c", ez->cigar[i]>>4, MB_CIGAR_STR[ez->cigar[i]&0xf]);
				fputc('\n', stderr);
			}
			mb_update_extra(km, h, &qseq[h->qs], &tseq[h->ts], mat, opt->q, opt->e, opt->flag, 0);
		}
	}
}

typedef struct {
	int32_t n, m;
	mb_hit_t *a;
} mb_hit_v;

static const mb_hit_t *mb_matesw_core(void *km, const mb_opt_t *opt, const l2b_t *l2b, const mb_pestat_t pes[4], const mb_hit_t *h0, int32_t r0,
	int32_t len, uint8_t *seq[2], l2b_meth_t mt0, mb_hit_v *h1, int32_t min_sc, ksw_extz_t *ez)
{
	int32_t dir, skip[4];
	int64_t pos5;
	const mb_hit_t *ret = 0;
	// find permitted orientation
	for (dir = 0; dir < 4; ++dir) skip[dir] = !!pes[dir].failed;
	if (skip[0] + skip[1] + skip[2] + skip[3] == 4) return 0; // no need to perform SW
	// perform SW
	pos5 = h0->rev? h0->te : h0->ts;
	for (dir = 0; dir < 4; ++dir) {
		int is_rev, is_larger;
		int64_t ts, te;
		if (skip[dir]) continue;
		is_rev = (dir>>1 != (dir&1)) ^ h0->rev; // whether to reverse complement the mate
		if (dir>>1 != (dir&1)) is_larger = dir>>1 ^ is_rev; // whether the mate has larger coordinate (FR or RF)
		else is_larger = dir>>1 ^ r0 ^ h0->rev; // FF or RR
		ts = (is_larger? pos5 + pes[dir].lo : pos5 - pes[dir].hi) - (!is_rev? 0 : len);
		te = (is_larger? pos5 + pes[dir].hi : pos5 - pes[dir].lo) + (!is_rev? len : 0);
		if (ts < 0) ts = 0;
		if (te > l2b->ctg[h0->tid].len) te = l2b->ctg[h0->tid].len;
		if (te - ts > len) {
			int64_t ts2 = ts, te2 = te;
			int32_t max_ug, max_i, n_good, n_kmer;
			uint8_t *ref;
			l2b_meth_t mt = mt0;
			mb_hit_t ht;
			ht.p = 0;
			ref = Kmalloc(km, uint8_t, te - ts);
			if (is_rev) mt = l2b_meth_rev(mt0);
			l2b_getseq(l2b, h0->tid, ts, te, ref);
			max_ug = mb_ungap(km, len, seq[is_rev], te - ts, ref, 7, mt, &max_i, &n_good, &n_kmer);
			if (max_ug >= 10 && max_ug >= len>>1 && n_good == 1) {
				ts2 = ts + max_i - len / 2;
				te2 = ts2 + len * 2;
				if (ts2 < ts) ts2 = ts;
				if (te2 > te) te2 = te;
			}
			if (max_ug >= 10 || max_ug >= n_kmer * 0.33)
				mb_matesw_align(km, opt, len, seq[is_rev], te2 - ts2, &ref[ts2 - ts], &ht, min_sc, mt, ez);
			if (ht.p) { // a good hit found
				ht.tid = h0->tid;
				ht.ts += ts2, ht.te += ts2;
				ht.rev = is_rev;
				if (is_rev) {
					int32_t qt = ht.qs;
					ht.qs = len - ht.qe;
					ht.qe = len - qt;
				}
				//fprintf(stderr, "X\t%s\tsc0=%d\tnew_sc=%d\tnew_ts=%ld\tn_cigar=%d\tcigar[0]=%d\n", l2b->ctg[ht.tid].name, h0->p->dp_max, ht.p->dp_max, (long)ht.ts, ht.p->n_cigar, ht.p->cigar[0]>>4);
				if (h1->n == h1->m) kom_grow(mb_hit_t, h1->a, h1->n, h1->m);
				h1->a[h1->n++] = ht;
				ret = &h1->a[h1->n - 1];
			}
			kfree(km, ref);
		}
	}
	return ret;
}

static int32_t mb_matesw(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], const mb_pairaux_t *paux0,
	int32_t qlen[2], char *const qseq[2], int32_t is_meth)
{
	int32_t i, r, n_add, n_res, max[2], max2[2], skip[2], min_sc[2];
	mb_hit_v ha[2];
	ksw_extz_t ez;
	mb128_t *a;
	uint8_t *qs[2][2];
	if (opt->max_rescue == 0) return 0;
	// precalculate the number of rescue candidates
	for (r = 0, n_res = 0; r < 2; ++r) {
		const mb_hit_t *h0 = hit[r];
		int32_t m, n = n_hit[r];
		for (i = 0, m = 0; i < n && m < opt->max_rescue; ++i)
			if (h0[i].proper_pair == 0 && h0[i].p->dp_max >= h0[0].p->dp_max - opt->pen_unpair * opt->a)
				++m;
		n_res += m;
	}
	if (n_res == 0) return 0;
	// collect rescue candidates; MUST match the loop above
	a = Kcalloc(km, mb128_t, n_res);
	for (r = 0, n_res = 0; r < 2; ++r) {
		const mb_hit_t *h0 = hit[r];
		int32_t m, n = n_hit[r];
		for (i = 0, m = 0; i < n && m < opt->max_rescue; ++i) {
			if (h0[i].proper_pair == 0 && h0[i].p->dp_max >= h0[0].p->dp_max - opt->pen_unpair * opt->a) {
				mb128_t *p = &a[n_res++];
				p->x = (uint64_t)h0[i].p->dp_max << 32 | h0[i].hash;
				p->y = i << 1 | r;
				++m;
			}
		}
	}
	radix_sort_mb128x(a, a + n_res);

	// prepare sequences for alignment
	qs[0][0] = Kcalloc(km, uint8_t, (qlen[0] + qlen[1]) * 2);
	qs[0][1] = qs[0][0] + qlen[0];
	qs[1][0] = qs[0][1] + qlen[0];
	qs[1][1] = qs[1][0] + qlen[1];
	for (r = 0; r < 2; ++r) {
		ha[r].n = ha[r].m = n_hit[r];
		ha[r].a = hit[r];
		for (i = 0; i < qlen[r]; ++i) {
			int32_t c = kom_nt4_table[(uint8_t)qseq[r][i]];
			qs[r][0][i] = c;
			qs[r][1][qlen[r] - 1 - i] = c < 4? 3 - c : 4;
		}
	}
	min_sc[0] = mb_hit_sum_score(km, n_hit[0], hit[0]) / opt->a - opt->pen_unpair;
	min_sc[1] = mb_hit_sum_score(km, n_hit[1], hit[1]) / opt->a - opt->pen_unpair;

	// do alignment
	max[0] = paux0->score, max2[0] = paux0->sub_sc, skip[0] = 0;
	max[1] = paux0->score, max2[1] = paux0->sub_sc, skip[1] = 0;
	memset(&ez, 0, sizeof(ez));
	ez.m_cigar = 16;
	ez.cigar = Kmalloc(km, uint32_t, ez.m_cigar);
	for (i = n_res - 1; i >= 0; --i) {
		int32_t sc, r = a[i].y&1, j = a[i].y>>1;
		const mb_hit_t *h0 = &ha[r].a[j], *h1;
		l2b_meth_t mt = !is_meth? L2B_METH_NONE : r == 0? L2B_METH_G2A : L2B_METH_C2T; // rescuing mate: r=0 rescues R2(G2A), r=1 rescues R1(C2T)
		if (skip[r]) continue;
		h1 = mb_matesw_core(km, opt, l2b, pes, h0, r, qlen[!r], qs[!r], mt, &ha[!r], min_sc[!r], &ez);
		if (h1) { // rescue successful
			sc = mb_pair_score(h0, h1, pes, opt->a);
			if (sc > max[r]) max2[r] = max[r], max[r] = sc;
			else if (sc > max2[r]) max2[r] = sc;
			if (max[r] == max2[r]) skip[r] = 1;
		}
	}

	kfree(km, ez.cigar);
	kfree(km, qs[0][0]);
	kfree(km, a);
	n_add = (ha[0].n - n_hit[0]) + (ha[1].n - n_hit[1]);
	for (r = 0; r < 2; ++r)
		n_hit[r] = ha[r].n, hit[r] = ha[r].a;
	return n_add;
}

void mb_pair(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], int32_t qlen[2], char *const qseq[2])
{
	int32_t r, i, dp_max_se[2], score_se, do_matesw, is_meth = !!(opt->flag & MB_F_METH);
	mb_pairaux_t paux;
	mb_hit_t *h[2];
	int32_t seed_ratio[2], min_seed_ratio;

	if (n_hit[0] == 0 && n_hit[1] == 0) return;
	seed_ratio[0] = n_hit[0] > 0? hit[0][0].seed_ratio : 255;
	seed_ratio[1] = n_hit[1] > 0? hit[1][0].seed_ratio : 255;
	min_seed_ratio = seed_ratio[0] < seed_ratio[1]? seed_ratio[0] : seed_ratio[1];
	mb_pair_hits(km, opt, l2b, n_hit, hit, pes, &paux);
	do_matesw = paux.n_pp > 0 && paux.score == paux.sub_sc? 0 : 1; // skip mate rescue if we see two equally best pairs
	if (do_matesw && opt->max_rescue > 0) {
		int32_t sub_diff = opt->a + opt->b > opt->q + opt->e? opt->a + opt->b : opt->q + opt->e;
		if (mb_matesw(km, opt, l2b, n_hit, hit, pes, &paux, qlen, qseq, is_meth) > 0) {
			for (r = 0; r < 2; ++r) {
				for (i = 0; i < n_hit[r]; ++i) {
					mb_hit_t *h = &hit[r][i];
					if (!h->rescued)
						h->n_sub = h->subsc = h->p->dp_max2 = 0;
				}
				mb_hit_sort(km, &n_hit[r], hit[r]);
				mb_set_parent(km, opt->mask_level, opt->mask_len, n_hit[r], hit[r], sub_diff, 0);
				mb_set_mapq(km, qlen[r], n_hit[r], hit[r], opt->min_chain_score, opt->a, mb_is_sr_mode(opt, qlen[r]), opt->max_sr_len);
			}
			mb_pair_hits(km, opt, l2b, n_hit, hit, pes, &paux); // pair again if new hits rescued
		}
	}
	if (paux.n_pp == 0) goto end_pairing;

	for (r = 0; r < 2; ++r) {
		for (dp_max_se[r] = 0, i = 0; i < n_hit[r]; ++i)
			if (dp_max_se[r] < hit[r][i].p->dp_max)
				dp_max_se[r] = hit[r][i].p->dp_max;
	}
	score_se = dp_max_se[0] + dp_max_se[1];

	h[0] = &hit[0][paux.i[0]];
	h[1] = &hit[1][paux.i[1]];
	if (paux.score >= score_se - opt->pen_unpair * opt->a) {
		int32_t mapq_pe, score2 = paux.sub_sc;
		double identity;
		mb_sync_high_cov(n_hit[0], hit[0]);
		mb_sync_high_cov(n_hit[1], hit[1]);
		identity = (double)(h[0]->mlen + h[1]->mlen) / (h[0]->blen + h[1]->blen);
		if ((h[0]->id != h[0]->parent || h[1]->id != h[1]->parent) && score2 < score_se - opt->pen_unpair * opt->a)
			score2 = score_se - opt->pen_unpair * opt->a;
		mapq_pe = (int)(6.02 * identity * identity * (paux.score - score2) / opt->a - 4.343 * log(paux.n_sub + 1) + .499);
		mapq_pe = (int)(mapq_pe * (1. - .5 * (h[0]->frac_high / 255. + h[1]->frac_high / 255.)) + .499);
		if (min_seed_ratio < 50) mapq_pe *= (double)min_seed_ratio * min_seed_ratio / 2500.0;
		if (mapq_pe > 60) mapq_pe = 60;
		if (mapq_pe == 0 && paux.score > score2) mapq_pe = 1;
		for (r = 0; r < 2; ++r) {
			if (h[r]->mapq < mapq_pe)
				h[r]->mapq = (int32_t)(.2 * h[r]->mapq + .8 * mapq_pe + .499);
			if (h[r]->id != h[r]->parent) { // then lift the paired hit to primary and update parent
				mb_hit_t *p = &hit[r][h[r]->parent];
				for (i = 0; i < n_hit[r]; ++i)
					if (hit[r][i].parent == p->id) // h[r]->parent is always set to h[r]->id
						hit[r][i].parent = h[r]->id;
				p->mapq = 0;
			}
			for (i = 0; i < n_hit[r]; ++i) { // handle other chimeric hits
				mb_hit_t *p = &hit[r][i], *q = h[r];
				if (q != p && p->id == p->parent) { // p is a chimeric hit that is not h[r]
					int32_t j, ol = p->qe <= q->qs || p->qs >= q->qe? 0 : (p->qe < q->qe? p->qe : q->qe) - (p->qs > q->qs? p->qs : q->qs);
					if (ol > opt->mask_level * (p->qe - p->qs)) { // if p overlaps with h[r] a lot, make it a secondary hit
						for (j = 0; j < n_hit[r]; ++j) // FIXME: quadratic time complexity, but almost never an issue on real data
							if (hit[r][j].parent == p->id)
								hit[r][j].parent = q->id;
						p->mapq = 0;
					}
				}
			}
		}
	}
end_pairing:
	mb_set_sam_pri(n_hit[0], hit[0], !!(opt->flag & MB_F_PRIMARY5));
	mb_set_sam_pri(n_hit[1], hit[1], !!(opt->flag & MB_F_PRIMARY5));
}
