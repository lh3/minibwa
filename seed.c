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

#define key_anchor(a) ((a).tpos)
KRADIX_SORT_INIT(mb_anchor, mb_anchor_t, key_anchor, 8)

static void mb_bwt_extend_eq(const mb_bwt_t *bwt, int32_t len, const uint8_t *q, mb_sai_t *p)
{
	uint32_t st = p->info>>32, en = (uint32_t)p->info;
	int32_t j, c;
	mb_sai_t ok[4];
	for (j = en; j < len; ++j) { // extend forwardly
		if (q[j] > 3) break;
		c = 3 - q[j];
		mb_bwt_extend(bwt, p, ok, 0);
		if (ok[c].size < p->size) break;
		*p = ok[c];
	}
	en = j;
	for (j = st - 1; j >= 0; --j) { // extend backwardly
		if (q[j] > 3) break;
		c = q[j];
		mb_bwt_extend(bwt, p, ok, 1);
		if (ok[c].size < p->size) break;
		*p = ok[c];
	}
	st = j + 1;
	p->info = (uint64_t)st<<32 | en;
}

static int32_t mb_bwt_seed_greedy(const mb_bwt_t *bwt, int32_t len, const uint8_t *q, int32_t x, int32_t min_len, int32_t max_intv, mb_sai_t *p)
{ // inspired by PMID:21209072 but more crude
	int32_t i, c;
	mb_sai_t ik, ok[4];
	memset(p, 0, sizeof(*p));
	if (q[x] > 3) return x + 1;
	mb_bwt_set_intv(bwt, q[x], &ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			mb_bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].size < max_intv && i - x >= min_len) {
				*p = ok[c];
				p->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	*p = ik;
	p->info = (uint64_t)x<<32 | len;
	return len;
}

void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{
	int64_t x = 0, i, n_a0;
	mb_sai_t p;

	v->n = 0;
	do {
		x = mb_bwt_smem(bwt, len, seq, x, min_len, 1, 1, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, v->a, v->n, v->m);
			v->a[v->n++] = p;
		}
	} while (x < len);

	n_a0 = v->n;
	for (i = 0; i < n_a0; ++i) {
		uint32_t st = v->a[i].info>>32, en = (uint32_t)v->a[i].info;
		if (en - st < min_len * 2 || v->a[i].size > max_sub_occ)
			continue;
		x = st;
		#if 1
		int32_t sub_min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
		do {
			x = mb_bwt_smem(bwt, len, seq, x, sub_min_len, v->a[i].size + 1, v->a[i].size + 1, &p);
			if (p.size > v->a[i].size) {
				Kgrow(km, mb_sai_t, v->a, v->n, v->m);
				v->a[v->n++] = p;
			}
		} while (x < en);
		#else
		do {
			x = mb_bwt_seed_greedy(bwt, en, seq, x, min_len, max_sub_occ * 2, &p);
			if (p.size > v->a[i].size) {
				int32_t to_add = 1;
				if (v->n > 0 && p.size == v->a[v->n-1].size && (uint32_t)p.info <= (uint32_t)v->a[v->n-1].info)
					to_add = 0;
				if (to_add) {
					mb_bwt_extend_eq(bwt, len, seq, &p);
					Kgrow(km, mb_sai_t, v->a, v->n, v->m);
					v->a[v->n++] = p;
				}
			}
		} while (x < en);
		#endif
	}
}

/************************
 * Get contig positions *
 ************************/

typedef struct { int64_t st, en; } anchor_aux_t;
typedef struct { int64_t a, i; } sa_aux_t;

static void process_batch(void *km, const mb_idx_t *idx, const anchor_aux_t *aux, int32_t m, const sa_aux_t *b, uint64_t *a, int32_t qlen, const mb_sai_v *u, mb_anchor_v *v)
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
			const l2b_ctg_t *ctg;
			mb_anchor_t *q;
			tid = l2b_intv2cid(idx->l2b, a[k], a[k] + len, &cst, &rev);
			rev = !!rev; // make sure rev is 0 or 1
			if (tid < 0) continue; // bridging boundaries
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

void mb_anchor(void *km, const mb_idx_t *idx, const mb_sai_v *u, int32_t qlen, int32_t max_occ, mb_anchor_v *v)
{
	const int batch_size = 20;
	int32_t n_aux, m, m_a;
	int64_t i, i0, j, k;
	uint64_t *a;
	sa_aux_t *b;
	anchor_aux_t *aux;

	v->n = 0;
	if (u->n == 0) return; // no anchors

	// sort by ::x[0] and then by ::size
	if (u->n > 1) radix_sort_mb_sai0(u->a, u->a + u->n);
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0]) {
			if (i - i0 > 1)
				radix_sort_mb_sais(&u->a[i0], &u->a[i]);
			i0 = i;
		}
	}

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
			process_batch(km, idx, aux, m, b, a, qlen, u, v);
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
	process_batch(km, idx, aux, m, b, a, qlen, u, v);
	kfree(km, b);
	kfree(km, a);
	kfree(km, aux);

	radix_sort_mb_anchor(v->a, v->a + v->n);
	for (i = 0; i < v->n; ++i) { // adjust mb_anchor_t::tpos
		mb_anchor_t *q = &v->a[i];
		const l2b_ctg_t *ctg = &idx->l2b->ctg[q->sid>>1];
		q->tpos -= ctg->off * 2 + ctg->len * (q->sid&1);
	}
}
