#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mbpriv.h"
#include "kalloc.h"

#define MIN_RATIO     0.8
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
		int64_t dist;
		if (seg_cnt[i] != 2) continue;
		off = seg_off[i];
		r[0] = mb_select_unique_se(n_hit[off + 0], hit[off + 0]);
		r[1] = mb_select_unique_se(n_hit[off + 1], hit[off + 1]);
		if (r[0] == 0 || r[1] == 0) continue;
		if (r[0]->tid != r[1]->tid) continue; // not on the same contig
		dist = r[0]->ts > r[1]->ts? r[0]->ts - r[1]->ts : r[1]->ts - r[0]->ts;
		dir = ((int32_t)r[0]->rev << 1 | (int32_t)r[1]->rev) ^ (r[0]->ts < r[1]->ts? 0 : 3);
		if (dist < opt->max_pe_ins) {
			if (is[dir].n == is[dir].m)
				Kgrow(km, uint64_t, is[dir].a, is[dir].n, is[dir].m);
			is[dir].a[is[dir].n++] = dist;
		}
	}
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%d, %d, %d, %d)\n", __func__, is[0].n, is[1].n, is[2].n, is[3].n);
	for (d = 0; d < 4; ++d) {
		mb_pestat_t *r = &pes[d];
		q = &is[d];
		int p25, p50, p75, x;
		if (q->n < MIN_DIR_CNT) {
			fprintf(stderr, "[M::%s] skip orientation %c%c as there are not enough pairs\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
			r->failed = 1;
			kfree(km, q->a);
			continue;
		} else fprintf(stderr, "[M::%s] analyzing insert size distribution for orientation %c%c...\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
		radix_sort_mb64(q->a, q->a + q->n);
		p25 = q->a[(int)(.25 * q->n + .499)];
		p50 = q->a[(int)(.50 * q->n + .499)];
		p75 = q->a[(int)(.75 * q->n + .499)];
		r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->low < 1) r->low = 1;
		r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
		fprintf(stderr, "[M::%s] (25, 50, 75) percentile: (%d, %d, %d)\n", __func__, p25, p50, p75);
		fprintf(stderr, "[M::%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, r->low, r->high);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
		fprintf(stderr, "[M::%s] mean and std.dev: (%.2f, %.2f)\n", __func__, r->avg, r->std);
		r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->high < r->avg + MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->low < 1) r->low = 1;
		fprintf(stderr, "[M::%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, r->low, r->high);
		kfree(km, q->a);
	}
	for (d = 0, max = 0; d < 4; ++d)
		max = max > is[d].n? max : is[d].n;
	for (d = 0; d < 4; ++d) {
		if (pes[d].failed == 0 && is[d].n < max * MIN_DIR_RATIO) {
			pes[d].failed = 1;
			fprintf(stderr, "[M::%s] skip orientation %c%c\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
		}
	}
}
