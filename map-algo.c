#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

mb_tbuf_t *mb_tbuf_init(int no_kalloc)
{
	mb_tbuf_t *b;
	b = kom_calloc(mb_tbuf_t, 1);
	if (!no_kalloc) b->km = km_init();
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
		h = (uint32_t)mb_hash64((mb_hash64(a[k].tpos) + mb_hash64(a[k].qpos)) ^ hash);
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
				if (rp->p && ri->p && (rp->tid != ri->tid || rp->ts != ri->ts || rp->te != ri->te || ol != min)) { // the last condition excludes identical hits after DP
					sci = ri->p->dp_max;
					rp->p->dp_max2 = rp->p->dp_max2 > sci? rp->p->dp_max2 : sci;
					if (rp->p->dp_max - ri->p->dp_max <= sub_diff) cnt_sub = 1;
				}
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

void mb_hit_sort(void *km, int *n_regs, mb_hit_t *r)
{
	int32_t i, n_aux, n = *n_regs;
	mb128_t *aux;
	mb_hit_t *t;

	if (n <= 1) return;
	aux = (mb128_t*)kmalloc(km, (size_t)n * 16);
	t = (mb_hit_t*)kmalloc(km, (size_t)n * sizeof(mb_hit_t));
	for (i = n_aux = 0; i < n; ++i) {
		if (r[i].inv || r[i].cnt > 0) { // squeeze out elements with cnt==0 (soft deleted)
			int score;
			if (r[i].p) score = r[i].p->dp_max;
			else score = r[i].score;
			aux[n_aux].x = (uint64_t)score << 32 | r[i].hash;
			aux[n_aux++].y = i;
		} else if (r[i].p) {
			free(r[i].p);
			r[i].p = 0;
		}
	}
	radix_sort_mb128x(aux, aux + n_aux);
	for (i = n_aux - 1; i >= 0; --i)
		t[n_aux - 1 - i] = r[aux[i].y];
	memcpy(r, t, sizeof(mb_hit_t) * n_aux);
	*n_regs = n_aux;
	kfree(km, aux);
	kfree(km, t);
}

void mb_filter_hits(const mb_opt_t *opt, int qlen, int *n_regs, mb_hit_t *regs)
{
	int i, k;
	for (i = k = 0; i < *n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		int flt = r->flt;
		if (r->p) {
			if (r->mlen < opt->min_chain_score) flt = 1;
			else if (r->p->dp_max < opt->min_dp_max) flt = 1;
			if (flt) free(r->p);
		}
		if (!flt) {
			if (k < i) regs[k++] = regs[i];
			else ++k;
		}
	}
	*n_regs = k;
}

int mb_squeeze_a(void *km, int n_regs, mb_hit_t *regs, mb_anchor_t *a)
{
	int i, as = 0;
	uint64_t *aux;
	aux = (uint64_t*)kmalloc(km, (size_t)n_regs * 8);
	for (i = 0; i < n_regs; ++i)
		aux[i] = (uint64_t)regs[i].as << 32 | i;
	radix_sort_mb64(aux, aux + n_regs);
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[(int32_t)aux[i]];
		if (r->as != as) {
			memmove(&a[as], &a[r->as], (size_t)r->cnt * sizeof(mb_anchor_t));
			r->as = as;
		}
		as += r->cnt;
	}
	kfree(km, aux);
	return as;
}

/*******************
 * Mapping quality *
 *******************/

static void mb_set_inv_mapq(void *km, int n_regs, mb_hit_t *regs)
{
	int i, n_aux;
	mb128_t *aux;
	if (n_regs < 3) return;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].inv) break;
	if (i == n_regs) return; // no inversion hits

	aux = Kmalloc(km, mb128_t, n_regs);
	for (i = n_aux = 0; i < n_regs; ++i)
		if (regs[i].parent == i || regs[i].parent < 0)
			aux[n_aux].y = i, aux[n_aux++].x = (uint64_t)regs[i].tid << 32 | regs[i].ts;
	radix_sort_mb128x(aux, aux + n_aux);

	for (i = 1; i < n_aux - 1; ++i) {
		mb_hit_t *inv = &regs[aux[i].y];
		if (inv->inv) {
			mb_hit_t *l = &regs[aux[i-1].y];
			mb_hit_t *r = &regs[aux[i+1].y];
			inv->mapq = l->mapq < r->mapq? l->mapq : r->mapq;
		}
	}
	kfree(km, aux);
}

void mb_set_mapq(void *km, int n_regs, mb_hit_t *regs, int min_chain_sc, int match_sc, int is_sr)
{
	static const float q_coef = 40.0f;
	int i;
	if (n_regs == 0) return;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (r->inv) {
			r->mapq = 0;
		} else if (r->parent == r->id) {
			int mapq, mapq_alt, subsc;
			float pen_s1 = r->score > 100? 1.0f : 0.01f * r->score;
			subsc = r->subsc > min_chain_sc? r->subsc : min_chain_sc;
			if (r->p && r->p->dp_max2 > 0 && r->p->dp_max > 0) {
				float x, identity = (float)r->mlen / r->blen;
				x = (float)r->p->dp_max2 * subsc / r->p->dp_max / r->score0;
				mapq = (int)(pen_s1 * identity * q_coef * (1.0f - x * x) * logf((float)r->p->dp_max / match_sc));
				// mapq_alt is close to the bwa-mem formula for short reads but minimap2 uses it only for long reads. What is happening?!
				mapq_alt = (int)(6.02f * identity * identity * (r->p->dp_max - r->p->dp_max2) / match_sc + .499f);
				mapq = mapq < mapq_alt? mapq : mapq_alt; // in case the long-read heuristic fails
			} else {
				float x = (float)subsc / r->score0;
				if (r->p) {
					float identity = (float)r->mlen / r->blen;
					mapq = (int)(pen_s1 * identity * q_coef * (1.0f - x) * logf((float)r->p->dp_max / match_sc));
				} else {
					mapq = (int)(pen_s1 * q_coef * (1.0f - x) * logf(r->score));
				}
			}
			mapq -= (int)(4.343f * logf(r->n_sub + 1) + .499f);
			mapq = mapq > 0? mapq : 0;
			r->mapq = mapq < 60? mapq : 60;
			if (r->p && r->p->dp_max > r->p->dp_max2 && r->mapq == 0) r->mapq = 1;
		} else r->mapq = 0;
	}
	mb_set_inv_mapq(km, n_regs, regs);
}

/************************
 * Core mapping routine *
 ************************/

static void mb_dbg_seed(int64_t n, const mb_sai_t *u, const char *qname)
{
	int64_t i;
	for (i = 0; i < n; ++i) {
		const mb_sai_t *p = &u[i];
		fprintf(stderr, "SD\t%s\t%d\t%d\t%ld\n", qname? qname : "*", (int32_t)(p->info>>32), (int32_t)p->info, (long)p->size);
	}
}

static void mb_dbg_anchor(const mb_idx_t *idx, int qlen, int64_t n, const mb_anchor_t *a, const char *qname)
{
	int64_t i;
	for (i = 0; i < n; ++i) {
		const mb_anchor_t *ai = &a[i];
		int rid = ai->sid >> 1;
		int rev = ai->sid & 1;
		int32_t qs = rev? qlen - 1 - ai->qpos : ai->qpos + 1 - ai->len;
		int64_t ts = ai->tpos + 1 - ai->len;
		fprintf(stderr, "AC\t%s\t%d\t%c\t%s\t%ld\t%d\n", qname? qname : "*", qs, "+-"[rev], idx->l2b->ctg[rid].name, (long)ts, ai->len);
	}
}

mb_hit_t *mb_map(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const char *seq0, int32_t *n_hit_, mb_tbuf_t *b, const char *qname)
{
	uint8_t *seq;
	uint32_t hash;
	int32_t i, n_hit;
	uint64_t *w;
	float chn_pen_gap, chn_pen_skip;
	mb_sai_v u = {0,0,0};
	mb_anchor_v v = {0,0,0};
	mb_anchor_t *a;
	mb_hit_t *hit;

	*n_hit_ = 0;
	seq = Kmalloc(b->km, uint8_t, qlen);
	for (i = 0; i < qlen; ++i)
		seq[i] = kom_nt4_table[(uint8_t)seq0[i]];

	hash  = qname? mb_hash_str(qname) : 0;
	hash ^= mb_hash64(qlen) + mb_hash64(opt->seed);
	hash  = mb_hash64(hash);

	chn_pen_gap = opt->chain_gap_scale * .01 * opt->min_len;
	chn_pen_skip = opt->chain_skip_scale * .01 * opt->min_len;
	mb_seed_intv(b->km, idx->bwt, qlen, seq, opt->min_len, opt->max_sub_occ, &u);
	if (kom_dbg_flag & MB_DBG_SEED) mb_dbg_seed(u.n, u.a, qname);
	mb_anchor(b->km, idx, &u, qlen, opt->max_occ, &v);
	if (kom_dbg_flag & MB_DBG_ANCHOR) mb_dbg_anchor(idx, qlen, v.n, v.a, qname);
	a = mb_lchain_dp(b->km, opt->max_gap, opt->max_gap, opt->bw, opt->max_chain_skip, opt->max_chain_iter,
					 opt->min_chain_score, chn_pen_gap, chn_pen_skip, v.n, v.a, &n_hit, &w);
	v.a = 0; v.n = v.m = 0; // ownership transferred to a
	kfree(b->km, u.a); // no longer needed

	hit = mb_gen_hit(b->km, hash, qlen, idx->l2b, n_hit, w, a);
	kfree(b->km, w);
	mb_set_parent(b->km, opt->mask_level, opt->mask_len, n_hit, hit, opt->a * 2 + opt->b, 0);
	mb_select_sub(b->km, opt->pri_ratio, opt->min_len * 2, opt->best_n, &n_hit, hit);
	if (!(opt->flag & MB_F_NO_ALN)) {
		hit = mb_align_skeleton(b->km, opt, idx, qlen, seq0, &n_hit, hit, a);
		mb_set_parent(b->km, opt->mask_level, opt->mask_len, n_hit, hit, opt->a * 2 + opt->b, 0);
		mb_select_sub(b->km, opt->pri_ratio, opt->min_len * 2, opt->best_n, &n_hit, hit);
	}
	mb_set_mapq(b->km, n_hit, hit, opt->min_chain_score, opt->a, !(opt->flag & MB_F_LONG));
	kfree(b->km, a);
	kfree(b->km, seq);
	*n_hit_ = n_hit;
	#if 0
	km_stat_t kmst;
	km_stat(b->km, &kmst);
	assert(kmst.n_blocks == kmst.n_cores);
	#endif
	return hit;
}
