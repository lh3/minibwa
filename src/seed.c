#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define key_sai0(a) ((a).x[0])
KRADIX_SORT_INIT(mb_sai0, mb_sai_t, key_sai0, 8)

#define key_sais(a) ((a).size)
KRADIX_SORT_INIT(mb_sais, mb_sai_t, key_sais, 8)

#define key_saii(a) ((a).info)
KRADIX_SORT_INIT(mb_saii, mb_sai_t, key_saii, 8)

#define key_anchor(a) ((a).tpos)
KRADIX_SORT_INIT(mb_anchor, mb_anchor_t, key_anchor, 8)

/***********
 * Seeding *
 ***********/

void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{
	int64_t x = 0, i, n_a0;
	mb_sai_t p;

	v->n = 0;
	do { // pass 1: standard SMEMs
		x = mb_bwt_smem(bwt, len, seq, x, min_len, 1, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, v->a, v->n, v->m);
			v->a[v->n++] = p;
		}
	} while (x < len);

	n_a0 = v->n;
	for (i = 0; i < n_a0; ++i) { // pass 2: sub-SMEMs
		int32_t sub_min_len;
		uint32_t st = v->a[i].info>>32, en = (uint32_t)v->a[i].info;
		if (en - st < min_len * 2 || v->a[i].size > max_sub_occ)
			continue;
		x = st;
		sub_min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
		do { // if two SMEMs have large overlaps, we may find the same sub intervals in both. A rare case not worth optimizing
			x = mb_bwt_smem(bwt, en, seq, x, sub_min_len, v->a[i].size + 1, &p);
			if (p.size > v->a[i].size) {
				Kgrow(km, mb_sai_t, v->a, v->n, v->m);
				v->a[v->n++] = p;
			}
		} while (x < en);
	}
}

void mb_seed_intv_batch(void *km, const mb_bwt_t *bwt, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{ // identical to mb_seed_intv() though the order of intervals is often different
	const int max_batch_size = 50;
	mb_smem_entry_t *s;
	int32_t i, j, n_s, *nv;

	// first pass: standard SMEMs
	s = Kcalloc(km, mb_smem_entry_t, max_batch_size);
	nv = Kcalloc(km, int32_t, n_seq);
	for (i = 0; i < n_seq; ++i) v[i].n = 0;
	for (i = 0; i < n_seq; i += max_batch_size) {
		int32_t en = i + max_batch_size < n_seq? i + max_batch_size : n_seq;
		for (j = i; j < en; ++j) {
			mb_smem_entry_t *t = &s[j - i];
			t->min_len = min_len;
			t->min_occ = 1;
			t->st = 0, t->en = len[j];
			t->q = seq[j];
			t->v = &v[j];
		}
		mb_bwt_smem_batch(km, bwt, en - i, s);
	}

	// second pass; sub-SMEMs
	for (i = 0; i < n_seq; ++i) nv[i] = v[i].n;
	for (i = n_s = 0; i < n_seq; ++i) {
		for (j = 0; j < nv[i]; ++j) {
			mb_smem_entry_t *t;
			uint32_t st = v[i].a[j].info>>32, en = (uint32_t)v[i].a[j].info;
			if (en - st < min_len * 2 || v[i].a[j].size > max_sub_occ)
				continue;
			t = &s[n_s++];
			t->min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
			t->min_occ = v[i].a[j].size + 1;
			t->st = st, t->en = en;
			t->q = seq[i];
			t->v = &v[i];
			if (n_s == max_batch_size) {
				mb_bwt_smem_batch(km, bwt, n_s, s);
				n_s = 0;
			}
		}
	}
	if (n_s > 0)
		mb_bwt_smem_batch(km, bwt, n_s, s);
	kfree(km, nv);
	kfree(km, s);
}

/*****************************
 * Seed/anchor deduplication *
 *****************************/

static void mb_seed_sort_dedup(mb_sai_v *u)
{
	int64_t i, i0, j;
	if (u->n <= 1) return;
	// sort by ::{x[0],size,info}
	radix_sort_mb_sai0(u->a, u->a + u->n); // sort by ::x[0]
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0]) {
			if (i - i0 > 1) {
				int64_t k, k0, n1 = i - i0;
				mb_sai_t *a1 = &u->a[i0];
				radix_sort_mb_sais(&u->a[i0], &u->a[i]); // sort by ::size
				kom_reverse(mb_sai_t, n1, a1);
				for (k = i0 + 1, k0 = i0; k <= i; ++k) {
					if (k == i || u->a[k0].size != u->a[k].size) {
						if (k - k0 > 1)
							radix_sort_mb_saii(&u->a[k0], &u->a[k]); // sort by ::info
						k0 = k;
					}
				}
			}
			i0 = i;
		}
	}
	// dedup
	for (i = 1, j = 0; i < u->n; ++i)
		if (!(u->a[i].x[0] == u->a[j].x[0] && u->a[i].size == u->a[j].size && u->a[i].info == u->a[j].info))
			u->a[++j] = u->a[i];
	u->n = j + 1;
}

/* Remove duplicated anchors. The two-round seeding algorithm may lead an
 * anchor precisely contained in a longer anchor. This routine filters out the
 * shorter anchor. This wouldn't happen to minimap2. */
static void mb_anchor_dedup(mb_anchor_v *v) // NB: assuming sorted by tpos
{
	const int max_back = 100; // to avoid quadratic behavior in the worst case
	int64_t i, j, k;
	for (i = 1; i < v->n; ++i) {
		mb_anchor_t *ai = &v->a[i];
		int64_t tsj, tsi = ai->tpos + 1 - ai->len;
		int32_t qsj, qsi = ai->qpos + 1 - ai->len;
		for (j = i - 1, k = 0; j >= 0 && k < max_back; --j, ++k) {
			mb_anchor_t *aj = &v->a[j];
			if (aj->sid != ai->sid) break;
			if (aj->tpos < tsi) break;
			tsj = aj->tpos + 1 - aj->len;
			qsj = aj->qpos + 1 - aj->len;
			if (tsj >= tsi) { // then j is contained in i
				if (tsj - tsi == qsj - qsi) aj->flt = 1;
			} else if (ai->tpos == aj->tpos) { // then i is contained in j
				if (ai->qpos == aj->qpos) ai->flt = 1;
			}
		}
	}
	for (i = j = 0; i < v->n; ++i)
		if (!v->a[i].flt) v->a[j++] = v->a[i];
	v->n = j;
}

/************************
 * Get contig positions *
 ************************/

typedef struct { int64_t st, en; } anchor_aux_t;
typedef struct { int64_t a, i; } sa_aux_t;

static void process_batch(void *km, const mb_idx_t *idx, const anchor_aux_t *aux, int32_t m, const sa_aux_t *b, uint64_t *a, int32_t qlen, l2b_meth_t mt, const mb_sai_v *u, mb_anchor_v *v)
{
	int64_t j, k;
	for (k = 0; k < m; ++k) a[k] = b[k].a;
	mb_bwt_sa_batch(km, idx->bwt, m, a);
	for (k = 0; k < m; ++k) {
		const anchor_aux_t *p = &aux[b[k].i];
		for (j = p->st; j < p->en; ++j) {
			int32_t qs = u->a[j].info>>32, qe = (int32_t)u->a[j].info;
			int32_t rev, len = qe - qs;
			int64_t tid, cst;
			l2b_meth_t mt_anchor;
			const l2b_ctg_t *ctg;
			mb_anchor_t *q;
			if (mt != L2B_METH_NONE) {
				tid = l2b_intv2cid_meth(idx->l2b, a[k], a[k] + len, &mt_anchor, &cst, &rev);
				if (tid < 0) continue;
				// R1(C2T): keep c2t_f(copy0) and g2a_r(copy2); R2(G2A): keep g2a_f(copy1) and c2t_r(copy3)
				if ((mt_anchor == mt) != (rev == 0)) continue; // filter
			} else {
				tid = l2b_intv2cid(idx->l2b, a[k], a[k] + len, &cst, &rev);
				if (tid < 0) continue;
			}
			rev = !!rev;
			ctg = &idx->l2b->ctg[tid];
			Kgrow(km, mb_anchor_t, v->a, v->n, v->m);
			q = &v->a[v->n++];
			memset(q, 0, sizeof(*q));
			q->sid = tid << 1 | rev;
			q->len = len;
			q->qpos = rev? qlen - 1 - qs : qs + len - 1;
			q->tpos = ctg->off * 2 + ctg->len * rev + cst + len - 1; // for sorting; will be adjusted later
		}
	}
}

static void mb_anchor_split_meth(void *km, const l2b_t *l2b, int32_t min_len, int32_t qlen, const uint8_t *qseq0, l2b_meth_t mt0, mb_anchor_v *v)
{
	int64_t i, m_a = v->n * 2, n_a = 0;
	int32_t max_len = 0;
	mb_anchor_t *a;
	uint8_t *tseq, *qseq2[2];

	for (i = 0; i < v->n; ++i)
		if (max_len < v->a[i].len)
			max_len = v->a[i].len;
	tseq = Kmalloc(km, uint8_t, max_len + qlen * 2);
	qseq2[0] = tseq + max_len;
	qseq2[1] = qseq2[0] + qlen;
	memcpy(qseq2[0], qseq0, qlen);
	for (i = 0; i < qlen; ++i)
		qseq2[1][qlen - 1 - i] = qseq0[i] > 3? 4 : 3 - qseq0[i];

	a = Kmalloc(km, mb_anchor_t, m_a);
	for (i = 0; i < v->n; ++i) {
		mb_anchor_t *p, *q = &v->a[i];
		const l2b_ctg_t *ctg = &l2b->ctg[q->sid>>1];
		int32_t rev = q->sid&1;
		int64_t tpos = q->tpos - (ctg->off * 2 + ctg->len * rev); // NB: requiring concatenated ::tpos!!
		int64_t ts = tpos + 1 - q->len;
		int32_t qs = q->qpos + 1 - q->len, j, j0;
		uint8_t t_allow, q_allow;
		l2b_meth_t mt;
		const uint8_t *qseq = qseq2[rev] + qs;
		l2b_getseq(l2b, q->sid>>1, ts, ts + q->len, tseq);
		mt = q->sid&1? l2b_meth_rev(mt0) : mt0;
		t_allow = mt == L2B_METH_C2T? 1 : 2;
		q_allow = mt == L2B_METH_C2T? 3 : 0;
		for (j0 = j = 0; j <= q->len; ++j) {
			if (j == q->len || tseq[j] == 4 || qseq[j] == 4 || (tseq[j] != qseq[j] && !(tseq[j] == t_allow && qseq[j] == q_allow))) {
				if (j - j0 >= min_len) {
					Kgrow(km, mb_anchor_t, a, n_a, m_a);
					p = &a[n_a++];
					*p = *q;
					p->len = j - j0;
					p->qpos = q->qpos - (q->len - j);
					p->tpos = q->tpos - (q->len - j);
				}
				j0 = j + 1;
			}
		}
	}
	kfree(km, tseq);
	Kgrow(km, mb_anchor_t, v->a, n_a, v->m);
	memcpy(v->a, a, n_a * sizeof(mb_anchor_t));
	v->n = n_a;
	kfree(km, a);
}

/* Converting seed intervals to anchors. This function batches small SA
 * intervals and calls mb_bwt_sa_batch() in process_batch(). With prefetch, the
 * strategy noticeably improves the performance. */
double mb_anchor(void *km, const mb_idx_t *idx, mb_sai_v *u, int32_t min_len, int32_t qlen, const uint8_t *qseq, l2b_meth_t mt, int32_t max_occ, mb_anchor_v *v)
{
	const int batch_size = 20;
	int32_t n_aux, m, m_a;
	int64_t i, i0, j, k;
	uint64_t *a;
	double seed_ratio = 1.0;
	sa_aux_t *b;
	anchor_aux_t *aux;

	v->n = 0;
	if (u->n == 0) return seed_ratio; // no anchors
	mb_seed_sort_dedup(u);

	for (i = 0, k = 0; i < u->n; ++i) // pre-calculate the size of v->a
		k += u->a[i].size < max_occ? u->a[i].size : max_occ;
	Kgrow(km, mb_anchor_t, v->a, k - 1, v->m); // preallocate

	for (i = 1, i0 = 0, n_aux = 0; i <= u->n; ++i) // pre-compute n_aux
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size)
			++n_aux, i0 = i;
	aux = Kmalloc(km, anchor_aux_t, n_aux);
	for (i = 1, i0 = 0, n_aux = 0; i <= u->n; ++i) // populate aux[]
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size)
			aux[n_aux].st = i0, aux[n_aux++].en = i, i0 = i;

	m_a = max_occ > batch_size? max_occ : batch_size; // max size of a[] and b[]
	a = Kmalloc(km, uint64_t, m_a);
	b = Kmalloc(km, sa_aux_t, m_a);
	for (i = 0, m = 0; i < n_aux; ++i) {
		const anchor_aux_t *p = &aux[i];
		const mb_sai_t *q = &u->a[p->st];
		if (q->size + m > batch_size) {
			process_batch(km, idx, aux, m, b, a, qlen, mt, u, v);
			m = 0;
		}
		if (q->size <= max_occ) { // get SA for all of them
			for (j = 0; j < q->size; ++j)
				b[m].a = q->x[0] + j, b[m++].i = i;
		} else { // sample up to max_occ
			int32_t n = 0;
			for (j = 0; j < q->size && n < max_occ; ++n) {
				int32_t step = (q->size - j) / (max_occ - n);
				if (step < 1) step = 1;
				b[m].a = q->x[0] + j, b[m++].i = i;
				j += step;
			}
		}
		assert(m <= m_a); // shouldn't happen!
	}
	process_batch(km, idx, aux, m, b, a, qlen, mt, u, v);
	kfree(km, b);
	kfree(km, a);
	kfree(km, aux);

	if (mt != L2B_METH_NONE && v->n > 0) {
		int64_t t0, t1;
		for (i = 0, t0 = 0; i < v->n; ++i) t0 += v->a[i].len;
		mb_anchor_split_meth(km, idx->l2b, min_len, qlen, qseq, mt, v);
		for (i = 0, t1 = 0; i < v->n; ++i) t1 += v->a[i].len;
		seed_ratio = (double)t1 / t0;
	}

	radix_sort_mb_anchor(v->a, v->a + v->n);
	for (i = 0; i < v->n; ++i) { // adjust mb_anchor_t::tpos
		mb_anchor_t *q = &v->a[i];
		const l2b_ctg_t *ctg = &idx->l2b->ctg[q->sid>>1];
		q->tpos -= ctg->off * 2 + ctg->len * (q->sid&1);
	}
	mb_anchor_dedup(v);
	return seed_ratio;
}

void mb_anchor_sort(const l2b_t *l2b, int64_t n_a, mb_anchor_t *a)
{
	int64_t i;
	if (n_a <= 1) return;
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos += ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
	radix_sort_mb_anchor(a, a + n_a);
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos -= ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
}
