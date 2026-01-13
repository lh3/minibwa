#include <stdlib.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"
#include "ksort.h"

#define key_128x(a) ((a).x)
KRADIX_SORT_INIT(mb128x, mb128_t, key_128x, 8)

#define key_64(a) (a)
KRADIX_SORT_INIT(mb64, uint64_t, key_64, 8)

/*****************
 * Thread buffer *
 *****************/

struct mb_tbuf_s {
	void *km;
};

mb_tbuf_t *mb_tbuf_init(void)
{
	mb_tbuf_t *b;
	b = kom_calloc(mb_tbuf_t, 1);
	if (!(kom_dbg_flag & MB_DBG_NO_KALLOC))
		b->km = km_init();
	return b;
}

void mb_tbuf_destroy(mb_tbuf_t *b)
{
	if (b->km) km_destroy(b->km);
	free(b);
}

/************************
 * Basic hit operations *
 ************************/

static inline uint64_t hash64(uint64_t x)
{
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

int32_t mb_cal_mblen(int32_t n, const mb_anchor_t *a, int32_t *blen_)
{
	int32_t i;
	int64_t mlen, blen;
	*blen_ = 0;
	if (n <= 0) return 0;
	mlen = blen = a[0].len;
	for (i = 1; i < n; ++i) {
		int span = a[i].len;
		int tl = (int32_t)a[i].tpos - (int32_t)a[i-1].tpos;
		int ql = (int32_t)a[i].qpos - (int32_t)a[i-1].qpos;
		blen += tl > ql? tl : ql;
		mlen += tl > span && ql > span? span : tl < ql? tl : ql;
	}
	*blen_ = blen;
	return mlen;
}

static void mb_cal_fuzzy_len(mb_hit_t *r, const mb_anchor_t *a)
{
	r->mlen = mb_cal_mblen(r->cnt, &a[r->as], &r->blen);
}

static inline void mb_hit_set_coor(mb_hit_t *r, int32_t qlen, const l2b_t *l2b, const mb_anchor_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
	int32_t k = r->as;
	const mb_anchor_t *ak0 = &a[k];
	const mb_anchor_t *ak1 = &a[k + r->cnt - 1];

	r->tid = ak0->sid>>1, r->rev = ak0->sid&1;
	r->ts = ak0->tpos + 1 - ak0->len;
	r->te = ak1->tpos + 1;
	if (!r->rev) { // forward strand
		r->qs = ak0->qpos + 1 - ak0->len;
		r->qe = ak1->qpos + 1;
	} else { // reverse strand
		r->qs = qlen - (ak1->qpos + 1);
		r->qe = qlen - (ak0->qpos + 1 - ak0->len);
	}
	mb_cal_fuzzy_len(r, a);
}

mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const l2b_t *l2b, int n_u, uint64_t *u, mb_anchor_t *a)
{ // convert chains to hits
	mb128_t *z, tmp;
	mb_hit_t *r;
	int i, k;

	if (n_u <= 0) return 0;

	// sort by score
	z = Kmalloc(km, mb128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		uint32_t h;
		h = (uint32_t)hash64((hash64(a[k].tpos) + hash64(a[k].qpos)) ^ hash);
		z[i].x = u[i] ^ h; // u[i] -- higher 32 bits: chain score; lower 32 bits: number of anchors
		z[i].y = (uint64_t)k << 32 | (int32_t)u[i];
		k += (int32_t)u[i];
	}
	radix_sort_mb128x(z, z + n_u);
	for (i = 0; i < n_u>>1; ++i) // reverse, s.t. larger score first
		tmp = z[i], z[i] = z[n_u-1-i], z[n_u-1-i] = tmp;

	// populate r[]
	r = (mb_hit_t*)calloc(n_u, sizeof(mb_hit_t));
	for (i = 0; i < n_u; ++i) {
		mb_hit_t *ri = &r[i];
		ri->id = i;
		ri->parent = MB_PARENT_UNSET;
		ri->score = ri->score0 = z[i].x >> 32;
		ri->hash = (uint32_t)z[i].x;
		ri->cnt = (int32_t)z[i].y;
		ri->as = z[i].y >> 32;
		mb_hit_set_coor(ri, qlen, l2b, a);
	}
	kfree(km, z);
	return r;
}

void mb_split_hit(mb_hit_t *r, mb_hit_t *r2, int n, int qlen, mb_anchor_t *a, const l2b_t *l2b)
{
	if (n <= 0 || n >= r->cnt) return;
	*r2 = *r;
	r2->id = -1;
	r2->sam_pri = 0;
	r2->p = 0;
	r2->split_inv = 0;
	r2->cnt = r->cnt - n;
	r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
	r2->as = r->as + n;
	if (r->parent == r->id) r2->parent = MB_PARENT_TMP_PRI;
	mb_hit_set_coor(r2, qlen, l2b, a);
	r->cnt -= r2->cnt;
	r->score -= r2->score;
	mb_hit_set_coor(r, qlen, l2b, a);
	r->split |= 1, r2->split |= 2;
}

void mb_sync_hits(void *km, int n_regs, mb_hit_t *regs)
{
	int *tmp, i, max_id = -1, n_tmp, n_pri = 0;
	if (n_regs <= 0) return;
	for (i = 0; i < n_regs; ++i)
		max_id = max_id > regs[i].id? max_id : regs[i].id;
	n_tmp = max_id + 1;
	tmp = (int*)kmalloc(km, n_tmp * sizeof(int));
	for (i = 0; i < n_tmp; ++i) tmp[i] = -1;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id >= 0) tmp[regs[i].id] = i;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		r->id = i;
		if (r->parent == MB_PARENT_TMP_PRI)
			r->parent = i;
		else if (r->parent >= 0 && tmp[r->parent] >= 0)
			r->parent = tmp[r->parent];
		else r->parent = MB_PARENT_UNSET;
	}
	kfree(km, tmp);
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id == regs[i].parent) {
			++n_pri;
			regs[i].sam_pri = (n_pri == 1);
		} else regs[i].sam_pri = 0;
}

/**********************************
 * Set primary and secondary hits *
 **********************************/

void mb_set_parent(void *km, float mask_level, int mask_len, int n, mb_hit_t *r, int sub_diff, int hard_mask_level)
{
	int i, j, k, *w;
	uint64_t *cov;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	cov = (uint64_t*)kmalloc(km, n * sizeof(uint64_t));
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		mb_hit_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe, n_cov = 0, uncov_len = 0;
		if (hard_mask_level) goto skip_uncov;
		for (j = 0; j < k; ++j) {
			mb_hit_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			if (ej <= si || sj >= ei) continue;
			if (sj < si) sj = si;
			if (ej > ei) ej = ei;
			cov[n_cov++] = (uint64_t)sj<<32 | ej;
		}
		if (n_cov == 0) {
			goto set_parent_test;
		} else {
			int j, x = si;
			radix_sort_mb64(cov, cov + n_cov);
			for (j = 0; j < n_cov; ++j) {
				if ((int)(cov[j]>>32) > x) uncov_len += (cov[j]>>32) - x;
				x = (int32_t)cov[j] > x? (int32_t)cov[j] : x;
			}
			if (ei > x) uncov_len += ei - x;
		}
skip_uncov:
		for (j = 0; j < k; ++j) {
			mb_hit_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe, min, max, ol;
			if (ej <= si || sj >= ei) continue;
			min = ej - sj < ei - si? ej - sj : ei - si;
			max = ej - sj > ei - si? ej - sj : ei - si;
			ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si);
			if ((float)ol / min - (float)uncov_len / max > mask_level && uncov_len <= mask_len) {
				int cnt_sub = 0, sci = ri->score;
				ri->parent = rp->parent;
				rp->subsc = rp->subsc > sci? rp->subsc : sci;
				if (ri->cnt >= rp->cnt) cnt_sub = 1;
				if (cnt_sub) ++rp->n_sub;
				break;
			}
		}
set_parent_test:
		if (j == k) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, cov);
	kfree(km, w);
}

void mb_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mb_hit_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		for (i = k = 0; i < n; ++i) {
			int p = r[i].parent;
			if (p == i || r[i].inv) {
				r[k++] = r[i];
			} else if ((r[i].score >= r[p].score * pri_ratio || r[i].score + min_diff >= r[p].score) && n_2nd < best_n) {
				if (!(r[i].qs == r[p].qs && r[i].qe == r[p].qe && r[i].tid == r[p].tid && r[i].ts == r[p].ts && r[i].te == r[p].te))
					r[k++] = r[i], ++n_2nd;
			}
		}
		if (k != n) mb_sync_hits(km, k, r);
		*n_ = k;
	}
}

mb_hit_t *mb_map(const mb_idx_t *idx, int64_t qlen, const char *seq0, int32_t *n_hit, mb_tbuf_t *b, const mb_mopt_t *opt, const char *qname)
{
	int64_t i;
	void *km = b->km;
	uint8_t *seq;

	seq = Kcalloc(km, uint8_t, qlen);
	for (i = 0; i < qlen; ++i)
		seq[i] = kom_nt4_table[(uint8_t)seq0[i]];
	kfree(km, seq);
	return 0;
}
