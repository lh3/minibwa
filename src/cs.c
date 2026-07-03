#include <assert.h>
#include <stdio.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"

static char *alloc_tmp(void *km, const mb_hit_t *r)
{
	int32_t i, min_tmp_len = 31;
	for (i = 0; i < (int)r->p->n_cigar; ++i) {
		int op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL)
			min_tmp_len = min_tmp_len > len + 5? min_tmp_len : len + 5;
	}
	return Kmalloc(km, char, min_tmp_len + 1);
}

static void write_indel_ds(void *km, kstring_t *str, int64_t len, const uint8_t *seq, int64_t ll, int64_t lr) // write an indel to ds; adapted from minigraph
{
	int64_t i;
	if (ll + lr >= len) {
		km_sprintf_lite(km, str, "[");
		for (i = 0; i < len; ++i)
			km_sprintf_lite(km, str, "%c", "acgtn"[seq[i]]);
		km_sprintf_lite(km, str, "]");
	} else {
		int64_t k = 0;
		if (ll > 0) {
			km_sprintf_lite(km, str, "[");
			for (i = 0; i < ll; ++i)
				km_sprintf_lite(km, str, "%c", "acgtn"[seq[k+i]]);
			km_sprintf_lite(km, str, "]");
			k += ll;
		}
		for (i = 0; i < len - lr - ll; ++i)
			km_sprintf_lite(km, str, "%c", "acgtn"[seq[k+i]]);
		k += len - lr - ll;
		if (lr > 0) {
			km_sprintf_lite(km, str, "[");
			for (i = 0; i < lr; ++i)
				km_sprintf_lite(km, str, "%c", "acgtn"[seq[k+i]]);
			km_sprintf_lite(km, str, "]");
		}
	}
}

void mb_write_cs_ds(void *km, kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mb_hit_t *r, int is_ds)
{
	int i, q_off, t_off, q_len = 0, t_len = 0;
	char *tmp;
	km_sprintf_lite(km, s, "%cs:Z:", is_ds? 'd' : 'c');
	for (i = 0; i < (int)r->p->n_cigar; ++i) {
		int op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		if (op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) {
			q_len += len, t_len += len;
		} else if (op == MB_CIGAR_INS) {
			q_len += len;
		} else if (op == MB_CIGAR_DEL || op == MB_CIGAR_N_SKIP) {
			t_len += len;
		}
	}
	tmp = alloc_tmp(km, r);
	for (i = q_off = t_off = 0; i < (int)r->p->n_cigar; ++i) {
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert((op >= MB_CIGAR_MATCH && op <= MB_CIGAR_N_SKIP) || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH);
		if (op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) {
			int l_tmp = 0;
			for (j = 0; j < len; ++j) {
				if (qseq[q_off + j] != tseq[t_off + j]) {
					if (l_tmp > 0) {
						km_sprintf_lite(km, s, ":%d", l_tmp);
						l_tmp = 0;
					}
					km_sprintf_lite(km, s, "*%c%c", "acgtn"[tseq[t_off + j]], "acgtn"[qseq[q_off + j]]);
				} else ++l_tmp;
			}
			if (l_tmp > 0)
				km_sprintf_lite(km, s, ":%d", l_tmp);
			q_off += len, t_off += len;
		} else if (op == MB_CIGAR_INS) {
			if (is_ds) {
				int z, ll, lr, y = q_off;
				for (z = 1; z <= len; ++z)
					if (y - z < 0 || qseq[y + len - z] != qseq[y - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (y + len + z >= q_len || qseq[y + len + z] != qseq[y + z])
						break;
				ll = z;
				km_sprintf_lite(km, s, "+");
				write_indel_ds(km, s, len, &qseq[y], ll, lr);
			} else {
				for (j = 0, tmp[len] = 0; j < len; ++j)
					tmp[j] = "acgtn"[qseq[q_off + j]];
				km_sprintf_lite(km, s, "+%s", tmp);
			}
			q_off += len;
		} else if (op == MB_CIGAR_DEL) {
			if (is_ds) {
				int z, ll, lr, x = t_off;
				for (z = 1; z <= len; ++z)
					if (x - z < 0 || tseq[x + len - z] != tseq[x - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (x + len + z >= t_len || tseq[x + z] != tseq[x + len + z])
						break;
				ll = z;
				km_sprintf_lite(km, s, "-");
				write_indel_ds(km, s, len, &tseq[x], ll, lr);
			} else {
				for (j = 0, tmp[len] = 0; j < len; ++j)
					tmp[j] = "acgtn"[tseq[t_off + j]];
				km_sprintf_lite(km, s, "-%s", tmp);
			}
			t_off += len;
		} else { // intron
			assert(len >= 2);
			km_sprintf_lite(km, s, "~%c%c%d%c%c", "acgtn"[tseq[t_off]], "acgtn"[tseq[t_off+1]],
				len, "acgtn"[tseq[t_off+len-2]], "acgtn"[tseq[t_off+len-1]]);
			t_off += len;
		}
	}
	kfree(km, tmp);
	assert(t_off == r->te - r->ts && q_off == r->qe - r->qs);
}

void mb_write_MD(void *km, kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mb_hit_t *r)
{
	int i, q_off, t_off, l_MD = 0;
	char *tmp;
	km_sprintf_lite(km, s, "MD:Z:");
	tmp = alloc_tmp(km, r);
	for (i = q_off = t_off = 0; i < (int)r->p->n_cigar; ++i) {
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert((op >= MB_CIGAR_MATCH && op <= MB_CIGAR_N_SKIP) || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH);
		if (op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) {
			for (j = 0; j < len; ++j) {
				if (qseq[q_off + j] != tseq[t_off + j]) {
					km_sprintf_lite(km, s, "%d%c", l_MD, "ACGTN"[tseq[t_off + j]]);
					l_MD = 0;
				} else ++l_MD;
			}
			q_off += len, t_off += len;
		} else if (op == MB_CIGAR_INS) {
			q_off += len;
		} else if (op == MB_CIGAR_DEL) {
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "ACGTN"[tseq[t_off + j]];
			km_sprintf_lite(km, s, "%d^%s", l_MD, tmp);
			l_MD = 0;
			t_off += len;
		} else if (op == MB_CIGAR_N_SKIP) {
			t_off += len;
		}
	}
	if (l_MD > 0) km_sprintf_lite(km, s, "%d", l_MD);
	kfree(km, tmp);
	assert(t_off == r->te - r->ts && q_off == r->qe - r->qs);
}
