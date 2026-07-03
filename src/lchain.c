#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksort.h"

/*
 * Anchor format in minibwa (mb_anchor_t):
 *   tid2: target (reference) sequence ID
 *   len:  length of the seed (q_span)
 *   qpos: query coordinate of the last base
 *   tpos: target coordinate of the last base
 *   qocc: query occurrence
 *   tocc: target occurrence
 *   flt:  filtered flag
 *
 * For chaining, we need:
 *   - Target position (tid2, tpos) for gap calculation on reference
 *   - Query position (qpos) for gap calculation on query
 *   - Seed length (len) as q_span for scoring
 */

static int64_t mb_chain_bk_end(int32_t max_drop, const mb128_t *z, const int32_t *f, const int64_t *p, int32_t *t, int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	if (i < 0 || t[i] != 0) return i;
	do {
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (s > max_s) max_s = s, max_i = i;
		else if (max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0);
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) // reset modified t[]
		t[i] = 0;
	return max_i;
}

static uint64_t *mb_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_sc, int32_t max_drop, int32_t *n_u_, int32_t *n_v_)
{
	mb128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	z = Kmalloc(km, mb128_t, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_mb128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mb_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				++n_v, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0)
				++n_u;
			else n_v = n_v0;
		}
	}
	u = Kmalloc(km, uint64_t, n_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mb_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				v[n_v++] = i, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0)
				u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
			else n_v = n_v0;
		}
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

static mb_anchor_t *compact_a(void *km, const l2b_t *l2b, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mb_anchor_t *a)
{
	mb_anchor_t *b;
	mb128_t *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	b = Kmalloc(km, mb_anchor_t, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	w = Kmalloc(km, mb128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[b[k].sid>>1];
		w[i].x = b[k].tpos + ctg->off * 2 + ctg->len * (b[k].sid&1);
		w[i].y = (uint64_t)k<<32 | i;
		k += (int32_t)u[i];
	}
	radix_sort_mb128x(w, w + n_u);
	u2 = Kmalloc(km, uint64_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mb_anchor_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mb_anchor_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, u2); kfree(km, w); kfree(km, a);
	return b;
}

// Compute chaining score between two anchors
static inline int32_t comput_sc(const mb_anchor_t *ai, const mb_anchor_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap)
{
	int64_t dq = ai->qpos - aj->qpos, dr, dd, dg, sc;
	if (dq <= 0 || dq > max_dist_y + ai->len) return INT32_MIN;
	// Check if on same target; use tid for comparison
	if (ai->sid != aj->sid) return INT32_MIN;
	dr = ai->tpos - aj->tpos;
	if (dr <= 0 || dq > max_dist_x + ai->len) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (dd > bw) return INT32_MIN;
	dg = dr < dq? dr : dq;
	// Use successor span for variable-length anchors (SMEMs)
	sc = ai->len < dg? ai->len : dg;
	if (dd) { // with a gap
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd;
		log_pen = dd >= 1? mb_log2(dd + 1) : 0.0f;
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

/* Input:
 *   a[].sid: target sequence ID (tid<<1|rev)
 *   a[].tpos: target coordinate of the last base
 *   a[].qpos: query coordinate of the last base
 *   a[].len:  seed length (q_span)
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
mb_anchor_t *mb_lchain_dp(void *km, const l2b_t *l2b, int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_sc, float chn_pen_gap,
						  int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
	int64_t *p, i, j, max_ii;
	uint64_t *u;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw) max_dist_y = bw;
	v = Kmalloc(km, int32_t, n);
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);

	// fill the score and backtrack arrays
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].len, n_skip = 0;
		for (j = i - 1; j >= 0 && j >= i - max_iter; --j) {
			int32_t sc;
			if (a[i].tpos - a[j].tpos >= max_dist_x + a[i].len) break;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap);
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;
		if (max_ii < 0 || a[i].tpos - a[max_ii].tpos > max_dist_x + a[i].len) {
			int32_t max = INT32_MIN;
			max_ii = -1; // the index of the best score. Rescue in case the best is missed due to the max_skip heuristic
			for (j = i - 1; j >= end_j && j >= 0; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].tpos - a[max_ii].tpos <= max_dist_x + a[i].len && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) mmax_f = max_f;
	}

	u = mb_chain_backtrack(km, n, f, p, v, t, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, t); kfree(km, f); kfree(km, p);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, l2b, n_u, u, n_v, v, a);
}
