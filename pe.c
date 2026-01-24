#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mbpriv.h"
#include "kalloc.h"

#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0

static const mb_hit_t *mb_select_unique_se(int32_t n_hit, const mb_hit_t *hit)
{
	int32_t j, n_pri = 0, mapq = 0, k = -1;
	for (j = 0; j < n_hit; ++j)
		if (hit[j].id == hit[j].parent) 
			++n_pri, mapq = hit[j].mapq, k = j;
	return n_pri == 1 && mapq >= 10? &hit[k] : 0;
}

void mb_pestat(void *km, const mb_opt_t *opt, int32_t n_frag, const int32_t *seg_off, const int32_t *seg_cnt, const int32_t *n_hit, mb_hit_t *const *hit, mb_pestat_t pes[4])
{
	int32_t i, d, max;
	struct { int32_t n, m; uint64_t *a; } is[4], *q;
	memset(is, 0, sizeof(is[0]) * 4);
	memset(pes, 0, sizeof(pes[0]) * 4);
	for (i = 0; i < n_frag; ++i) {
		const mb_hit_t *r[2];
		int32_t off, dir;
		int64_t dist, pos[2];
		if (seg_cnt[i] != 2) continue;
		off = seg_off[i];
		r[0] = mb_select_unique_se(n_hit[off + 0], hit[off + 0]);
		r[1] = mb_select_unique_se(n_hit[off + 1], hit[off + 1]);
		if (r[0] == 0 || r[1] == 0) continue;
		if (r[0]->tid != r[1]->tid) continue; // not on the same contig
		pos[0] = r[0]->rev? r[0]->te : r[0]->ts; // 5'-end of the read
		pos[1] = r[1]->rev? r[1]->te : r[1]->ts;
		dist = pos[0] > pos[1]? pos[0] - pos[1] : pos[1] - pos[0];
		dir = ((int32_t)r[0]->rev << 1 | (int32_t)r[1]->rev) ^ (pos[0] < pos[1]? 0 : 3);
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
		r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->low < 1) r->low = 1;
		r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] (25, 50, 75) percentile: (%d, %d, %d); mean and std.dev: (%.2f, %.2f)\n",
				__func__, "FR"[d>>1&1], "FR"[d&1], p25, p50, p75, r->avg, r->std);
		r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->high < r->avg + MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->low < 1) r->low = 1;
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] low and high boundaries for proper pairs: (%d, %d)\n", __func__, "FR"[d>>1&1], "FR"[d&1], r->low, r->high);
		kfree(km, q->a);
	}
}

#define MB_SQRT1_2 0.707106781186547524401

typedef struct {
	int32_t score, sub_sc, n_sub, n_pp;
	int32_t i[2];
} mb_pairaux_t;

void mb_pair(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], mb_pairaux_t *ret)
{
	int32_t r, i, k, n_pa, y[4], n_pp = 0, m_pp = 0;
	mb128_t *pa, *pp = 0; // pp: proper pairs

	ret->i[0] = ret->i[1] = ret->score = ret->sub_sc = -1, ret->n_sub = 0;
	if (n_hit[0] == 0 || n_hit[1] == 0) return;
	pa = Kcalloc(km, mb128_t, n_hit[0] + n_hit[1]);
	for (r = n_pa = 0; r < 2; ++r) {
		for (i = 0; i < n_hit[r]; ++i) {
			mb128_t *p = &pa[n_pa++];
			const mb_hit_t *h = &hit[r][i];
			p->x = l2b->ctg[h->tid].off + (h->rev? h->te : h->ts);
			p->y = (uint64_t)i << 2 | (uint64_t)h->rev << 1 | r;
		}
	}
	radix_sort_mb128x(pa, pa + n_pa);

	y[0] = y[1] = y[2] = y[3] = -1;
	for (i = 0; i < n_pa; ++i) {
		mb128_t *pi = &pa[i];
		const mb_hit_t *hi = &hit[pi->y&1][pi->y>>2];
		for (r = 0; r < 2; ++r) {
			int which, dir = r << 1 | (pi->y>>1&1);
			if (pes[dir].failed) continue; // invalid orientation
			which = r << 1 | ((pi->y&1) ^ 1);
			if (y[which] < 0) continue; // no previous hit
			for (k = y[which]; k >= 0; --k) {
				mb128_t *pk = &pa[k], *q;
				const mb_hit_t *hk = &hit[pk->y&1][pk->y>>2];
				int64_t dist;
				double ns, s;
				if ((pk->y&3) != which) continue;
				if (hi->tid != hk->tid) break;
				dist = pi->x - pk->x;
				if (dist > pes[dir].high) break;
				if (dist < pes[dir].low) continue;
				ns = (dist - pes[dir].avg) / pes[dir].std; // normalized score
				s = hk->p->dp_max + hi->p->dp_max + .721 * log(2. * erfc(fabs(ns) * MB_SQRT1_2)) * opt->a; // .721 = 1/log(4)
				if (s < 0.) s = 0.;
				if (n_pp == m_pp) Kgrow(km, mb128_t, pp, n_pp, m_pp);
				q = &pp[n_pp++];
				q->y = (pk->y&1) == 0? (uint64_t)k << 32 | i : (uint64_t)i << 32 | k; // upper bits: index to read[0]; lower: to read[1]
				q->x = (uint64_t)(s + .499) << 32 | ((hk->hash ^ hi->hash) & 0xffffffffULL);
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
			if (q->x > max) max = q->x, ret->i[0] = q->y>>32, ret->i[1] = (uint32_t)q->y;
			else if (q->x > max2) max2 = q->x;
		}
		ret->score = max>>32, ret->sub_sc = max2>>32;
		for (i = 0; i < n_pp; ++i)
			if (pp[i].x>>32 <= ret->score - tmp)
				ret->n_sub++;
	}
	kfree(km, pp);
	kfree(km, pa);
}
