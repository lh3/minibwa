#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mbpriv.h"
#include "kommon.h"
#include "kalloc.h"
#include "ksw2.h"

static inline void update_max_zdrop(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int *max_zdrop, int pos[2][2])
{
	if (score < *max) {
		int li = i - *max_i;
		int lj = j - *max_j;
		int diff = li > lj? li - lj : lj - li;
		int z = *max - score - diff * e;
		if (z > *max_zdrop) {
			*max_zdrop = z;
			pos[0][0] = *max_i, pos[0][1] = i;
			pos[1][0] = *max_j, pos[1][1] = j;
		}
	} else *max = score, *max_i = i, *max_j = j;
}

static int mm_test_zdrop(void *km, const mb_opt_t *opt, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat, int32_t is_sr)
{
	uint32_t k;
	int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
	int pos[2][2] = {{-1, -1}, {-1, -1}}, q_len, t_len;

	// find the score and the region where score drops most along diagonal
	for (k = 0, score = 0; k < n_cigar; ++k) {
		uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == MB_CIGAR_MATCH) {
			for (l = 0; l < len; ++l) {
				score += mat[tseq[i + l] * 5 + qseq[j + l]];
				update_max_zdrop(score, i+l, j+l, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
			}
			i += len, j += len;
		} else if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL || op == MB_CIGAR_N_SKIP) {
			score -= opt->q + opt->e * len;
			if (op == MB_CIGAR_INS) j += len;
			else i += len;
			update_max_zdrop(score, i, j, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
		}
	}

	// test if there is an inversion in the most dropped region
	q_len = pos[1][1] - pos[1][0], t_len = pos[0][1] - pos[0][0];
	if (!is_sr && max_zdrop > opt->zdrop_inv && q_len < opt->max_gap && t_len < opt->max_gap) {
		uint8_t *qseq2;
		void *qp;
		int q_off, t_off;
		qseq2 = Kmalloc(km, uint8_t, q_len);
		for (i = 0; i < q_len; ++i) {
			int c = qseq[pos[1][1] - i - 1];
			qseq2[i] = c >= 4? 4 : 3 - c;
		}
		qp = ksw_ll_qinit(km, 2, q_len, qseq2, 5, mat);
		score = ksw_ll_i16(qp, t_len, tseq + pos[0][0], opt->q, opt->e, &q_off, &t_off);
		kfree(km, qseq2);
		kfree(km, qp);
		if (score >= opt->min_chain_score * opt->a && score >= opt->min_dp_max * opt->a)
			return 2; // there is a potential inversion
	}
	return max_zdrop > opt->zdrop? 1 : 0;
}

static void mb_fix_cigar(mb_hit_t *r, const uint8_t *qseq, const uint8_t *tseq, int *qshift, int *tshift)
{
	mb_extra_t *p = r->p;
	int32_t toff = 0, qoff = 0, to_shrink = 0;
	uint32_t k;
	*qshift = *tshift = 0;
	if (p->n_cigar <= 1) return;
	for (k = 0; k < p->n_cigar; ++k) { // indel left alignment
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (len == 0) to_shrink = 1;
		if (op == MB_CIGAR_MATCH) {
			toff += len, qoff += len;
		} else if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL) {
			if (k > 0 && k < p->n_cigar - 1 && (p->cigar[k-1]&0xf) == 0 && (p->cigar[k+1]&0xf) == 0) {
				int l, prev_len = p->cigar[k-1] >> 4;
				if (op == MB_CIGAR_INS) {
					for (l = 0; l < prev_len; ++l)
						if (qseq[qoff - 1 - l] != qseq[qoff + len - 1 - l])
							break;
				} else {
					for (l = 0; l < prev_len; ++l)
						if (tseq[toff - 1 - l] != tseq[toff + len - 1 - l])
							break;
				}
				if (l > 0)
					p->cigar[k-1] -= l<<4, p->cigar[k+1] += l<<4, qoff -= l, toff -= l;
				if (l == prev_len) to_shrink = 1;
			}
			if (op == MB_CIGAR_INS) qoff += len;
			else toff += len;
		} else if (op == MB_CIGAR_N_SKIP) {
			toff += len;
		}
	}
	assert(qoff == r->qe - r->qs && toff == r->te - r->ts);
	for (k = 0; k < p->n_cigar - 2; ++k) { // fix CIGAR like 5I6D7I
		if ((p->cigar[k]&0xf) > 0 && (p->cigar[k]&0xf) + (p->cigar[k+1]&0xf) == 3) {
			uint32_t l, s[3] = {0,0,0};
			for (l = k; l < p->n_cigar; ++l) { // count number of adjacent I and D
				uint32_t op = p->cigar[l]&0xf;
				if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL || p->cigar[l]>>4 == 0)
					s[op] += p->cigar[l] >> 4;
				else break;
			}
			if (s[1] > 0 && s[2] > 0 && l - k > 2) { // turn to a single I and a single D
				p->cigar[k]   = s[1]<<4|MB_CIGAR_INS;
				p->cigar[k+1] = s[2]<<4|MB_CIGAR_DEL;
				for (k += 2; k < l; ++k)
					p->cigar[k] &= 0xf;
				to_shrink = 1;
			}
			k = l;
		}
	}
	if (to_shrink) { // squeeze out zero-length operations
		int32_t l = 0;
		for (k = 0; k < p->n_cigar; ++k) // squeeze out zero-length operations
			if (p->cigar[k]>>4 != 0)
				p->cigar[l++] = p->cigar[k];
		p->n_cigar = l;
		for (k = l = 0; k < p->n_cigar; ++k) // merge two adjacent operations if they are the same
			if (k == p->n_cigar - 1 || (p->cigar[k]&0xf) != (p->cigar[k+1]&0xf))
				p->cigar[l++] = p->cigar[k];
			else p->cigar[k+1] += p->cigar[k]>>4<<4; // add length to the next CIGAR operator
		p->n_cigar = l;
	}
	if ((p->cigar[0]&0xf) == MB_CIGAR_INS || (p->cigar[0]&0xf) == MB_CIGAR_DEL) { // get rid of leading I or D
		int32_t l = p->cigar[0] >> 4;
		if ((p->cigar[0]&0xf) == MB_CIGAR_INS) {
			if (r->rev) r->qe -= l;
			else r->qs += l;
			*qshift = l;
		} else r->ts += l, *tshift = l;
		--p->n_cigar;
		memmove(p->cigar, p->cigar + 1, p->n_cigar * 4);
	}
}

static void mm_update_cigar_eqx(mb_hit_t *r, const uint8_t *qseq, const uint8_t *tseq) // written by @armintoepfer
{
	uint32_t n_EQX = 0, n_X = 0;
	uint32_t k, l, m, cap, toff = 0, qoff = 0, n_M = 0;
	mb_extra_t *p;
	if (r->p == 0) return;
	for (k = 0; k < r->p->n_cigar; ++k) {
		uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
		if (op == MB_CIGAR_MATCH) {
			while (len > 0) {
				for (l = 0; l < len && qseq[qoff + l] == tseq[toff + l] && qseq[qoff + l] < 4; ++l) {} // run of "="; N is a mismatch (cf. e41830b for NM)
				if (l > 0) { ++n_EQX; len -= l; toff += l; qoff += l; }
				for (l = 0; l < len && !(qseq[qoff + l] == tseq[toff + l] && qseq[qoff + l] < 4); ++l) {} // run of "X" (includes N<=>N)
				if (l > 0) { ++n_EQX; ++n_X; len -= l; toff += l; qoff += l; }
			}
			++n_M;
		} else if (op == MB_CIGAR_INS) {
			qoff += len;
		} else if (op == MB_CIGAR_DEL) {
			toff += len;
		} else if (op == MB_CIGAR_N_SKIP) {
			toff += len;
		}
	}
	// update in-place only if every M op is a single pure "=" run (no "X",
	// hence no mismatch or N); otherwise the emission pass below is needed
	if (n_X == 0) {
		for (k = 0; k < r->p->n_cigar; ++k) {
			uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
			if (op == MB_CIGAR_MATCH) r->p->cigar[k] = len << 4 | MB_CIGAR_EQ_MATCH;
		}
		return;
	}
	// allocate new storage
	cap = r->p->n_cigar + (n_EQX - n_M) + sizeof(mb_extra_t);
	kom_roundup32(cap);
	p = (mb_extra_t*)calloc(cap, 4);
	memcpy(p, r->p, sizeof(mb_extra_t));
	p->cap = cap;
	// update cigar while copying
	toff = qoff = m = 0;
	for (k = 0; k < r->p->n_cigar; ++k) {
		uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
		if (op == MB_CIGAR_MATCH) {
			while (len > 0) {
				// match ("="); N is a mismatch, so require both bases < 4
				for (l = 0; l < len && qseq[qoff + l] == tseq[toff + l] && qseq[qoff + l] < 4; ++l) {}
				if (l > 0) p->cigar[m++] = l << 4 | MB_CIGAR_EQ_MATCH;
				len -= l;
				toff += l, qoff += l;
				// mismatch ("X"); includes N<=>N
				for (l = 0; l < len && !(qseq[qoff + l] == tseq[toff + l] && qseq[qoff + l] < 4); ++l) {}
				if (l > 0) p->cigar[m++] = l << 4 | MB_CIGAR_X_MISMATCH;
				len -= l;
				toff += l, qoff += l;
			}
			continue;
		} else if (op == MB_CIGAR_INS) {
			qoff += len;
		} else if (op == MB_CIGAR_DEL) {
			toff += len;
		} else if (op == MB_CIGAR_N_SKIP) {
			toff += len;
		}
		p->cigar[m++] = r->p->cigar[k];
	}
	p->n_cigar = m;
	free(r->p);
	r->p = p;
}

void mb_update_extra(void *km, mb_hit_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e, uint64_t opt_flag, int log_gap)
{
	uint32_t k, l;
	int32_t qshift, tshift, toff = 0, qoff = 0, len4;
	double s = 0.0, max = 0.0;
	mb_extra_t *p = r->p;
	if (p == 0) return;
	mb_fix_cigar(r, qseq, tseq, &qshift, &tshift);
	qseq += qshift, tseq += tshift; // qseq and tseq may be shifted due to the removal of leading I/D
	r->blen = r->mlen = 0;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == MB_CIGAR_MATCH) {
			int n_ambi = 0, n_diff = 0;
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l];
				if (ct > 3 || cq > 3) ++n_ambi;
				else if (ct != cq) n_diff += (mat[ct * 5 + cq] < 0);
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
			}
			r->blen += len - n_ambi, r->mlen += len - (n_ambi + n_diff), p->n_ambi += n_ambi;
			toff += len, qoff += len;
		} else if (op == MB_CIGAR_INS) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			if (log_gap) s -= q + (double)e * mb_log2(1.0 + len);
			else s -= q + e;
			if (s < 0) s = 0;
			qoff += len;
		} else if (op == MB_CIGAR_DEL) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			if (log_gap) s -= q + (double)e * mb_log2(1.0 + len);
			else s -= q + e;
			if (s < 0) s = 0;
			toff += len;
		}
	}
	p->dp_max0 = p->dp_max = (int32_t)(max + .499);
	assert(qoff == r->qe - r->qs && toff == r->te - r->ts);
	if (opt_flag & MB_F_EQX) mm_update_cigar_eqx(r, qseq, tseq); // NB: it has to be called here as changes to qseq and tseq are not returned
	if (opt_flag & (MB_F_WRITE_DS|MB_F_WRITE_CS|MB_F_WRITE_MD)) {
		kstring_t str = {0,0,0};
		str.m = 256;
		str.s = kmalloc(km, str.m);
		if (opt_flag & (MB_F_WRITE_DS|MB_F_WRITE_CS))
			mb_write_cs_ds(km, &str, tseq, qseq, r, !!(opt_flag & MB_F_WRITE_DS));
		else
			mb_write_MD(km, &str, tseq, qseq, r);
		r->p->cs = 1;
		len4 = r->p->n_cigar + sizeof(mb_extra_t)/4 + (str.l + 1 + 3) / 4;
		if (len4 > r->p->cap) {
			r->p->cap = len4;
			r->p = (mb_extra_t*)realloc(r->p, r->p->cap * 4);
		}
		memcpy(&r->p->cigar[r->p->n_cigar], str.s, str.l + 1);
		kfree(km, str.s);
	}
}

static void mb_enlarge_cigar(mb_hit_t *r, uint32_t n_cigar) // TODO: this calls the libc realloc()
{
	if (n_cigar == 0) return;
	if (r->p == 0) {
		uint32_t cap = n_cigar + sizeof(mb_extra_t)/4;
		kom_roundup32(cap);
		r->p = (mb_extra_t*)calloc(cap, 4);
		r->p->cap = cap;
	} else if (r->p->n_cigar + n_cigar + sizeof(mb_extra_t)/4 > r->p->cap) {
		r->p->cap = r->p->n_cigar + n_cigar + sizeof(mb_extra_t)/4;
		kom_roundup32(r->p->cap);
		r->p = (mb_extra_t*)realloc(r->p, r->p->cap * 4);
	}
}

void mb_append_cigar(mb_hit_t *r, uint32_t n_cigar, const uint32_t *cigar)
{
	mb_extra_t *p;
	if (n_cigar == 0) return;
	mb_enlarge_cigar(r, n_cigar);
	p = r->p;
	if (p->n_cigar > 0 && (p->cigar[p->n_cigar-1]&0xf) == (cigar[0]&0xf)) { // same CIGAR op at the boundary
		p->cigar[p->n_cigar-1] += cigar[0]>>4<<4;
		if (n_cigar > 1) memcpy(p->cigar + p->n_cigar, cigar + 1, (n_cigar - 1) * 4);
		p->n_cigar += n_cigar - 1;
	} else {
		memcpy(p->cigar + p->n_cigar, cigar, n_cigar * 4);
		p->n_cigar += n_cigar;
	}
}

static inline int32_t mb_min_int32(int32_t a, int32_t b)
{
	return a < b? a : b;
}

static inline int32_t max_bw_from_mm(const mb_opt_t *opt, int32_t mm)
{
	int32_t x = mm * (opt->a + opt->b), max2 = 0, max1 = 0;
	if (x >= opt->q  + opt->e)  max1 = (x - opt->q  + opt->e  - 1) / opt->e;
	if (x >= opt->q2 + opt->e2) max2 = (x - opt->q2 + opt->e2 - 1) / opt->e2;
	return max1 > max2? max1 : max2;
}

static void mb_align_pair(void *km, const mb_opt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq,
						  const int8_t *mat, int w, int end_bonus, int zdrop, int ksw_flag, ksw_extz_t *ez)
{
	const int max_bw_adj_len = 100; // don't adjust bandwidth if sequences are too long
	int32_t j, n_mm = -1;
	if ((opt->b_ts != 0 && opt->b != opt->b_ts) || (opt->flag&MB_F_METH))
		ksw_flag |= KSW_EZ_GENERIC_SC;
	if ((ksw_flag & KSW_EZ_EXTZ_ONLY) && tlen >= qlen) { // ungapped extension
		ksw_reset_extz(ez);
		for (j = 0, ez->score = ez->max = 0; j < qlen; ++j) {
			ez->score += mat[tseq[j] * 5 + qseq[j]];
			n_mm += (tseq[j] > 3 || qseq[j] > 3 || mat[tseq[j] * 5 + qseq[j]] < 0);
			if (ez->max < ez->score) ez->max = ez->score, ez->max_q = ez->max_t = j;
		}
		if (n_mm <= 2) {
			ez->mqe = ez->score, ez->mqe_t = qlen - 1;
			if (ez->mqe + end_bonus >= ez->max) {
				ez->reach_end = 1;
				ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, MB_CIGAR_MATCH, qlen);
				return;
			}
		}
	} else if (qlen == tlen && !(ksw_flag & KSW_EZ_EXTZ_ONLY)) { // ungapped alignment
		int32_t max_gapped_score = (qlen - 2) * opt->a - 2 * (opt->q + opt->e);
		ksw_reset_extz(ez);
		for (j = 0, ez->score = 0; j < qlen; ++j) {
			ez->score += mat[tseq[j] * 5 + qseq[j]];
			n_mm += (tseq[j] > 3 || qseq[j] > 3 || mat[tseq[j] * 5 + qseq[j]] < 0);
		}
		if (n_mm <= 3 || ez->score > max_gapped_score) {
			ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, MB_CIGAR_MATCH, qlen);
			return;
		}
	}

	if (n_mm >= 0 && mb_min_int32(qlen, tlen) < max_bw_adj_len) { // n_mm >= 0 => ungapped alignment attempted
		int32_t max_bw;
		max_bw = max_bw_from_mm(opt, n_mm);
		if (w > max_bw + 4) w = max_bw + 4;
	}

	if (opt->max_sw_mat > 0 && (int64_t)tlen * qlen > opt->max_sw_mat) { // too much memory; skip alignment
		ksw_reset_extz(ez);
		ez->zdropped = 1;
	} else if (opt->q == opt->q2 && opt->e == opt->e2) { // affine gap
		ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, w, zdrop * opt->a, end_bonus, ksw_flag, ez);
	} else { // dual affine gap
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, w, zdrop * opt->a, end_bonus, ksw_flag, ez);
		//fprintf(stderr, "D2\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", tlen, qlen, !!(ksw_flag&KSW_EZ_EXTZ_ONLY), ez->max_t, ez->max_q, ez->max, ez->zdropped);
	}
	if (kom_dbg_flag & MB_DBG_ALN_SEQ) {
		int i;
		fprintf(stderr, "===> q=(%d,%d), e=(%d,%d), bw=%d, ksw_flag=0x%x, zdrop=%d, end_bonus=%d <===\n", opt->q, opt->q2, opt->e, opt->e2, w, ksw_flag, opt->zdrop, end_bonus);
		for (i = 0; i < tlen; ++i) fputc("ACGTN"[tseq[i]], stderr);
		fputc('\n', stderr);
		for (i = 0; i < qlen; ++i) fputc("ACGTN"[qseq[i]], stderr);
		fputc('\n', stderr);
		fprintf(stderr, "score=%d, max=%d, cigar=", ez->score, ez->max);
		for (i = 0; i < ez->n_cigar; ++i) fprintf(stderr, "%d%c", ez->cigar[i]>>4, MB_CIGAR_STR[ez->cigar[i]&0xf]);
		fprintf(stderr, "\n");
	}
}

static int *collect_long_gaps(void *km, int as1, int cnt1, mb_anchor_t *a, int min_gap, int *n_)
{
	int i, n, *K;
	*n_ = 0;
	for (i = 1, n = 0; i < cnt1; ++i) { // count the number of gaps longer than min_gap
		int64_t gap = (a[as1 + i].qpos - a[as1 + i - 1].qpos) - (a[as1 + i].tpos - a[as1 + i - 1].tpos);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return 0;
	K = Kmalloc(km, int, n);
	for (i = 1, n = 0; i < cnt1; ++i) { // store the positions of long gaps
		int64_t gap = (a[as1 + i].qpos - a[as1 + i - 1].qpos) - (a[as1 + i].tpos - a[as1 + i - 1].tpos);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	*n_ = n;
	return K;
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mb_anchor_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{ // this function deals with e.g. 1000I20M1000D
	int max_st, max_en, n, i, k, max, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	max = 0, max_st = max_en = -1;
	for (k = 0;; ++k) { // traverse long gaps
		int gap, l, n_ins = 0, n_del = 0, qs, max_diff = 0, max_diff_l = -1;
		int64_t ts;
		if (k == n || k >= max_en) {
			if (max_en > 0)
				for (i = K[max_st]; i < K[max_en]; ++i)
					a[as1 + i].flag |= MB_SEED_IGNORE;
			max = 0, max_st = max_en = -1;
			if (k == n) break;
		}
		i = K[k];
		gap = (a[as1 + i].qpos - a[as1 + i - 1].qpos) - (a[as1 + i].tpos - a[as1 + i - 1].tpos);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = a[as1 + i - 1].qpos;
		ts = a[as1 + i - 1].tpos;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if (a[as1 + j].qpos - a[as1 + j].len - qs > max_ext_len || a[as1 + j].tpos - a[as1 + j].len - ts > max_ext_len) break;
			gap = (a[as1 + j].qpos - a[as1 + j - 1].qpos) - (a[as1 + j].tpos - a[as1 + j - 1].tpos);
			if (gap > 0) n_ins += gap;
			else n_del += -gap;
			diff = n_ins + n_del - abs(n_ins - n_del);
			if (max_diff < diff)
				max_diff = diff, max_diff_l = l;
		}
		if (max_diff > diff_thres && max_diff > max)
			max = max_diff, max_st = k, max_en = max_diff_l;
	}
	kfree(km, K);
}

static void mm_filter_bad_seeds_alt(void *km, int as1, int cnt1, mb_anchor_t *a, int min_gap, int max_ext)
{ // this function deals with e.g. 1000I20M2000I
	int n, k, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	for (k = 0; k < n;) {
		int i = K[k], l;
		int gap1 = (a[as1 + i].qpos - a[as1 + i - 1].qpos) - (a[as1 + i].tpos - a[as1 + i - 1].tpos);
		int64_t te1 = a[as1 + i].tpos;
		int32_t qe1 = a[as1 + i].qpos;
		int32_t left_len = a[as1 + i].len;
		gap1 = gap1 > 0? gap1 : -gap1;
		for (l = k + 1; l < n; ++l) {
			int j = K[l], gap2, m;
			if (a[as1 + j].qpos - qe1 > max_ext || a[as1 + j].tpos - te1 > max_ext) break;
			gap2 = (a[as1 + j].qpos - a[as1 + j - 1].qpos) - (a[as1 + j].tpos - a[as1 + j - 1].tpos);
			int64_t m_t = a[as1 + j - 1].tpos - te1 + left_len;
			int32_t m_q = a[as1 + j - 1].qpos - qe1 + left_len;
			m = m_t < m_q? m_t : m_q;
			gap2 = gap2 > 0? gap2 : -gap2;
			if (m > gap1 + gap2) break;
			te1 = a[as1 + j].tpos;
			qe1 = a[as1 + j].qpos;
			left_len = a[as1 + j].len;
			gap1 = gap2;
		}
		if (l > k + 1) {
			int j, end = K[l - 1];
			for (j = K[k]; j < end; ++j)
				a[as1 + j].flag |= MB_SEED_IGNORE;
			a[as1 + end].flag |= MB_SEED_LONG_JOIN;
		}
		k = l;
	}
	kfree(km, K);
}

static void mm_fix_bad_ends(const mb_hit_t *r, const mb_anchor_t *a, int bw, int min_match, int32_t *as, int32_t *cnt)
{
	int32_t i, l, m;
	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 3) return;
	m = l = a[r->as].len;
	for (i = r->as + 1; i < r->as + r->cnt - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].len;
		if (a[i].flag & MB_SEED_LONG_JOIN) break;
		lr = a[i].tpos - a[i-1].tpos;
		lq = a[i].qpos - a[i-1].qpos;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
	*cnt = r->as + r->cnt - *as;
	m = l = a[r->as + r->cnt - 1].len;
	for (i = r->as + r->cnt - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].len;
		if (a[i+1].flag & MB_SEED_LONG_JOIN) break;
		lr = a[i+1].tpos - a[i].tpos - a[i+1].len + a[i].len;
		lq = a[i+1].qpos - a[i].qpos - a[i+1].len + a[i].len;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
}

static void mb_max_stretch(const mb_hit_t *r, const mb_anchor_t *a, int32_t *as, int32_t *cnt)
{ // find the max ungapped chain
	int32_t i, score, max_score, len, max_i, max_len;

	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 2) return;

	max_score = -1, max_i = -1, max_len = 0;
	score = a[r->as].len, len = 1;
	for (i = r->as + 1; i < r->as + r->cnt; ++i) {
		int32_t lr = a[i].tpos - a[i-1].tpos;
		int32_t lq = a[i].qpos - a[i-1].qpos;
		if (lq == lr) { // ungapped
			score += lq < a[i].len? lq : a[i].len; // in theory, "lq > a[i].len" should always stand
			++len;
		} else { // a gap
			if (score > max_score)
				max_score = score, max_len = len, max_i = i - len;
			score = a[i].len, len = 1;
		}
	}
	if (score > max_score)
		max_score = score, max_len = len, max_i = i - len;
	*as = max_i, *cnt = max_len;
}

static void mb_align1(void *km, const mb_opt_t *opt, const mb_idx_t *mi, int qlen, uint8_t *qseq0[2], l2b_meth_t mt, mb_hit_t *r, mb_hit_t *r2, int n_a, mb_anchor_t *a, ksw_extz_t *ez)
{
	int32_t is_sr, max_back, rev = a[r->as].sid&1, as1, cnt1;
	uint8_t *tseq = 0, *qseq;
	int32_t i, bw, bw_long, dropped = 0, ksw_flag = 0;
	int64_t tid = a[r->as].sid >> 1, l;
	int64_t ts0, te0, ts1, te1, ts, te; // ts0/te0: range of extracted sequence; ts1/te1: range of alignment; ts/te: moving temporary
	int32_t qs0, qe0, qs1, qe1, qs, qe;
	int8_t mat[25];

	is_sr = mb_is_sr_mode(opt, qlen);
	max_back = is_sr? 0 : 10; // for long reads, allow up to 10bp "edges" from chain ends
	r2->cnt = 0;
	if (r->cnt == 0) return;
	if (r->rev) mt = l2b_meth_rev(mt);
	ksw_gen_nt4_mat(mat, opt->a, opt->b, opt->b_ts, opt->b_ambi, (int)mt);
	bw = (int)(opt->bw * 1.5 + 1.);
	if (!is_sr) {
		bw_long = (int)(opt->bw_long * 1.5 + 1.);
		if (bw_long < bw) bw_long = bw;
	} else bw_long = bw; // disable long gap in the short-read mode

	if (is_sr) {
		mb_max_stretch(r, a, &as1, &cnt1);
	} else {
		mm_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2, &as1, &cnt1);
		mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
		mm_filter_bad_seeds_alt(km, as1, cnt1, a, 30, opt->max_gap>>1);
	}
	ts = a[as1].tpos + 1 - a[as1].len + mb_min_int32(a[as1].len>>1, max_back);
	qs = a[as1].qpos + 1 - a[as1].len + mb_min_int32(a[as1].len>>1, max_back);
	te = a[as1+cnt1-1].tpos + 1 - mb_min_int32(a[as1+cnt1-1].len>>1, max_back);
	qe = a[as1+cnt1-1].qpos + 1 - mb_min_int32(a[as1+cnt1-1].len>>1, max_back);
	assert(cnt1 > 0);

	if (kom_dbg_flag & MB_DBG_AN_POS) {
		for (i = 0; i < r->cnt; ++i) {
			int32_t gap = i == 0? 0 : (a[r->as+i].qpos - a[r->as+i-1].qpos) - (a[r->as+i].tpos - a[r->as+i-1].tpos);
			fprintf(stderr, "AF\t%d\t%s\t%ld\t%d\t%d\t%ld\n", r->as, mi->l2b->ctg[tid].name, (long)a[r->as + i].tpos, a[r->as + i].qpos, gap, (long)a[r->as + i].len);
		}
	}

	/* Look for the start and end of regions to perform DP. This sounds easy
	 * but is in fact tricky. Excessively small regions lead to unnecessary
	 * clippings and lose alignable sequences. Excessively large regions
	 * occasionally lead to large overlaps between two chains and may cause
	 * loss of alignments in corner cases. */
	if (is_sr) {
		qs0 = 0, qe0 = qlen;
		l = qs;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		ts0 = ts - l > 0? ts - l : 0;
		l = qlen - qe;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		te0 = te + l < mi->l2b->ctg[tid].len? te + l : mi->l2b->ctg[tid].len;
	} else {
		// compute ts0 and qs0
		ts0 = a[r->as].tpos + 1 - a[r->as].len;
		qs0 = a[r->as].qpos + 1 - a[r->as].len;
		if (ts0 < 0) ts0 = 0;
		assert(qs0 >= 0); // this should never happen, or it is logic error
		ts1 = qs1 = 0;
		if (qs > 0 && ts > 0) {
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs1 = qs1 > qs - l? qs1 : qs - l;
			qs0 = qs0 < qs1? qs0 : qs1; // at least include qs0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < ts? l : ts;
			ts1 = ts1 > ts - l? ts1 : ts - l;
			ts0 = ts0 < ts1? ts0 : ts1;
			ts0 = ts0 < ts? ts0 : ts;
		} else ts0 = ts, qs0 = qs;
		// compute te0 and qe0
		te0 = a[r->as + r->cnt - 1].tpos + 1;
		qe0 = a[r->as + r->cnt - 1].qpos + 1;
		te1 = mi->l2b->ctg[tid].len, qe1 = qlen;
		if (qe < qlen && te < mi->l2b->ctg[tid].len) {
			l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
			qe1 = qe1 < qe + l? qe1 : qe + l;
			qe0 = qe0 > qe1? qe0 : qe1; // at least include qe0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < mi->l2b->ctg[tid].len - te? l : mi->l2b->ctg[tid].len - te;
			te1 = te1 < te + l? te1 : te + l;
			te0 = te0 > te1? te0 : te1;
		} else te0 = te, qe0 = qe;
	}

	assert(te0 > ts0);
	tseq = Kmalloc(km, uint8_t, te0 - ts0);

	if (qs > 0 && ts > 0) { // left extension; probably the condition can be changed to "qs > qs0 && ts > ts0"
		qseq = &qseq0[rev][qs0];
		l2b_getseq(mi->l2b, tid, ts0, ts, tseq);
		mb_seq_rev(qs - qs0, qseq);
		mb_seq_rev(ts - ts0, tseq);
		mb_align_pair(km, opt, qs - qs0, qseq, ts - ts0, tseq, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, ksw_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0) {
			mb_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->reach_end? ez->mqe : ez->max;
		}
		ts1 = ts - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);
		mb_seq_rev(qs - qs0, qseq);
	} else ts1 = ts, qs1 = qs;
	te1 = ts, qe1 = qs;
	assert(qs1 >= 0 && ts1 >= 0);

	{ // adding exact match on the first unfiltered anchor
		te = te1 = a[as1].tpos + 1 - mb_min_int32(a[as1].len>>1, max_back);
		qe = qe1 = a[as1].qpos + 1 - mb_min_int32(a[as1].len>>1, max_back);
		assert(te - ts == qe - qs && te >= ts);
		uint32_t cigar0 = (te - ts) << 4 | MB_CIGAR_MATCH;
		mb_append_cigar(r, 1, &cigar0);
		r->p->dp_score += opt->a * (te - ts);
		ts = te, qs = qe;
	}

	for (i = 1; i < cnt1; ++i) { // gap filling
		const mb_anchor_t *ai = &a[as1 + i];
		if ((ai->flag & MB_SEED_IGNORE) && i != cnt1 - 1) continue;
		te1 = ai->tpos + 1 - mb_min_int32(ai->len>>1, max_back);
		qe1 = ai->qpos + 1 - mb_min_int32(ai->len>>1, max_back);
		if (i == cnt1 - 1 || (a[as1+i].flag&MB_SEED_LONG_JOIN) || (qe1 - qs >= opt->min_ksw_len && te1 - ts >= opt->min_ksw_len)) { // gap filling
			int32_t j, bw1 = bw_long, zdrop_code;
			int64_t d1 = 0; // distance from (qe1,te1) to trim
			// compute ts and te
			if (ai->len > opt->min_len * 2) {
				d1 = te1 - (ai->tpos + 1 - ai->len); // distance to the start of the anchor
				d1 = d1 < qe1 - qs? d1 : qe1 - qs;
				d1 = d1 < te1 - ts? d1 : te1 - ts;
				d1 -= opt->min_len;
				if (d1 < opt->min_len) d1 = 0;
			}
			te = te1 - d1, qe = qe1 - d1;
			// update bandwidth
			if (a[as1+i].flag & MB_SEED_LONG_JOIN)
				bw1 = qe - qs > te - ts? qe - qs : te - ts;
			// perform alignment
			qseq = &qseq0[rev][qs];
			l2b_getseq(mi->l2b, tid, ts, te, tseq);
			mb_align_pair(km, opt, qe - qs, qseq, te - ts, tseq, mat, bw1, -1, opt->zdrop, ksw_flag|KSW_EZ_APPROX_MAX, ez); // first pass: with approximate Z-drop
			// test Z-drop and inversion Z-drop
			if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat, is_sr)) != 0)
				mb_align_pair(km, opt, qe - qs, qseq, te - ts, tseq, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, ksw_flag, ez); // second pass: lift approximate
			if (kom_dbg_flag & MB_DBG_AN_POS) fprintf(stderr, "AD\t%d\t%ld\t%ld\t%d\t%d\t%d\t%d\n", r->as, (long)ts, (long)te, qs, qe, zdrop_code, ez->zdropped);
			// update CIGAR
			if (ez->n_cigar > 0)
				mb_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				int32_t mlen, blen;
				if (!r->p) {
					assert(ez->n_cigar == 0);
					uint32_t cap = sizeof(mb_extra_t)/4;
					kom_roundup32(cap);
					r->p = (mb_extra_t*)calloc(cap, 4);
					r->p->cap = cap;
				}
				for (j = i - 1; j >= 0; --j)
					if (a[as1 + j].tpos <= ts + ez->max_t)
						break;
				dropped = 1;
				if (j < 0) j = 0;
				r->p->dp_score += ez->max;
				te1 = ts + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				mlen = mb_cal_mblen(cnt1 - (j + 1), &a[as1 + j + 1], &blen);
				if (mlen >= opt->min_chain_score) { // TODO: check if this is correct
					mb_split_hit(r, r2, as1 + j + 1 - r->as, qlen, a, mi->l2b);
					if (zdrop_code == 2) r2->split_inv = 1;
				}
				break;
			} else r->p->dp_score += ez->score;
			if (d1 > 0) {
				uint32_t cigar0 = d1 << 4 | MB_CIGAR_MATCH;
				mb_append_cigar(r, 1, &cigar0);
				r->p->dp_score += opt->a * d1;
				te = te1, qe = qe1;
			}
			ts = te, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && te < te0) { // right extension
		qseq = &qseq0[rev][qe];
		l2b_getseq(mi->l2b, tid, te, te0, tseq);
		mb_align_pair(km, opt, qe0 - qe, qseq, te0 - te, tseq, mat, bw, opt->end_bonus, opt->zdrop, ksw_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mb_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->reach_end? ez->mqe : ez->max;
		}
		te1 = te + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);
	}
	assert(qe1 <= qlen);

	r->ts = ts1, r->te = te1;
	if (!rev) r->qs = qs1, r->qe = qe1;
	else r->qs = qlen - qe1, r->qe = qlen - qs1;

	assert(te1 - ts1 <= te0 - ts0);
	if (r->p) {
		l2b_getseq(mi->l2b, tid, ts1, te1, tseq);
		qseq = &qseq0[r->rev][qs1];
		mb_update_extra(km, r, qseq, tseq, mat, opt->q, opt->e, opt->flag, !is_sr);
	}

	kfree(km, tseq);
}

static int mb_align1_inv(void *km, const mb_opt_t *opt, const mb_idx_t *mi, int qlen, uint8_t *qseq0[2], l2b_meth_t mt, const mb_hit_t *r1, const mb_hit_t *r2, mb_hit_t *r_inv, ksw_extz_t *ez)
{ // NB: this doesn't work with the qstrand mode
	int tl, ql, score, ret = 0, q_off, t_off;
	uint8_t *tseq, *qseq;
	int8_t mat[25];
	void *qp;

	memset(r_inv, 0, sizeof(mb_hit_t));
	if (!(r1->split&1) || !(r2->split&2)) return 0;
	if (r1->id != r1->parent && r1->parent != MB_PARENT_TMP_PRI) return 0;
	if (r2->id != r2->parent && r2->parent != MB_PARENT_TMP_PRI) return 0;
	if (r1->tid != r2->tid || r1->rev != r2->rev) return 0;
	ql = r1->rev? r1->qs - r2->qe : r2->qs - r1->qe;
	tl = r2->ts - r1->te;
	if (ql < opt->min_chain_score || ql > opt->max_gap) return 0;
	if (tl < opt->min_chain_score || tl > opt->max_gap) return 0;

	if (!r1->rev) mt = l2b_meth_rev(mt); // TODO: check if this is correct
	ksw_gen_nt4_mat(mat, opt->a, opt->b, opt->b_ts, opt->b_ambi, (int)mt);

	tseq = (uint8_t*)kmalloc(km, tl);
	l2b_getseq(mi->l2b, r1->tid, r1->te, r2->ts, tseq);
	qseq = r1->rev? &qseq0[0][r2->qe] : &qseq0[1][qlen - r2->qs];

	mb_seq_rev(ql, qseq);
	mb_seq_rev(tl, tseq);
	qp = ksw_ll_qinit(km, 2, ql, qseq, 5, mat);
	score = ksw_ll_i16(qp, tl, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, qp);
	mb_seq_rev(ql, qseq);
	mb_seq_rev(tl, tseq);
	if (score < opt->min_dp_max * opt->a) goto end_align1_inv;
	q_off = ql - (q_off + 1), t_off = tl - (t_off + 1);
	mb_align_pair(km, opt, ql - q_off, qseq + q_off, tl - t_off, tseq + t_off, mat, (int)(opt->bw * 1.5), -1, opt->zdrop, KSW_EZ_EXTZ_ONLY, ez);
	if (ez->n_cigar == 0) goto end_align1_inv; // should never be here
	mb_append_cigar(r_inv, ez->n_cigar, ez->cigar);
	r_inv->p->dp_score = ez->max;
	r_inv->id = -1;
	r_inv->parent = MB_PARENT_UNSET;
	r_inv->inv = 1;
	r_inv->rev = !r1->rev;
	r_inv->tid = r1->tid;
	if (r_inv->rev == 0) {
		r_inv->qs = r2->qe + q_off;
		r_inv->qe = r_inv->qs + ez->max_q + 1;
	} else {
		r_inv->qe = r2->qs - q_off;
		r_inv->qs = r_inv->qe - (ez->max_q + 1);
	}
	r_inv->ts = r1->te + t_off;
	r_inv->te = r_inv->ts + ez->max_t + 1;
	mb_update_extra(km, r_inv, &qseq[q_off], &tseq[t_off], mat, opt->q, opt->e, opt->flag, mb_is_sr_mode(opt, qlen));
	ret = 1;
end_align1_inv:
	kfree(km, tseq);
	return ret;
}

static inline mb_hit_t *mb_insert_reg(const mb_hit_t *r, int i, int *n_regs, mb_hit_t *regs)
{
	regs = (mb_hit_t*)realloc(regs, (*n_regs + 1) * sizeof(mb_hit_t));
	if (i + 1 != *n_regs)
		memmove(&regs[i + 2], &regs[i + 1], sizeof(mb_hit_t) * (*n_regs - i - 1));
	regs[i + 1] = *r;
	++*n_regs;
	return regs;
}

static inline void mb_count_gaps(const mb_hit_t *r, int32_t *n_gap_, int32_t *n_gapo_)
{
	uint32_t i;
	int32_t n_gapo = 0, n_gap = 0;
	*n_gap_ = *n_gapo_ = -1;
	if (r->p == 0) return;
	for (i = 0; i < r->p->n_cigar; ++i) {
		int32_t op = r->p->cigar[i] & 0xf, len = r->p->cigar[i] >> 4;
		if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL)
			++n_gapo, n_gap += len;
	}
	*n_gap_ = n_gap, *n_gapo_ = n_gapo;
}

static double mb_event_identity(const mb_hit_t *r)
{
	int32_t n_gap, n_gapo;
	if (r->p == 0) return -1.0f;
	mb_count_gaps(r, &n_gap, &n_gapo);
	return (double)r->mlen / (r->blen + r->p->n_ambi - n_gap + n_gapo);
}

static int32_t mb_recal_max_dp(const mb_hit_t *r, double b2, int32_t match_sc)
{
	uint32_t i;
	int32_t n_gap = 0, n_mis;
	double gap_cost = 0.0;
	if (r->p == 0) return -1;
	for (i = 0; i < r->p->n_cigar; ++i) {
		int32_t op = r->p->cigar[i] & 0xf, len = r->p->cigar[i] >> 4;
		if (op == MB_CIGAR_INS || op == MB_CIGAR_DEL) {
			gap_cost += b2 * mb_log2(1.0 + len);
			n_gap += len;
		}
	}
	n_mis = r->blen + r->p->n_ambi - r->mlen - n_gap;
	return (int32_t)(match_sc * (r->mlen - b2 * n_mis - gap_cost) + .499);
}

void mb_update_dp_max(int qlen, int n_regs, mb_hit_t *regs, double frac, int a, int b)
{
	int32_t max = -1, max2 = -1, i, max_i = -1, max2_i = -1;
	double div, b2;
	if (n_regs < 2) return;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (r->p == 0) continue;
		if (r->p->dp_max > max) max2 = max, max2_i = max_i, max = r->p->dp_max, max_i = i;
		else if (r->p->dp_max > max2) max2 = r->p->dp_max, max2_i = i;
	}
	if (max_i < 0 || max2_i < 0) return;
	if (regs[max_i].qe - regs[max_i].qs < qlen * frac) return;
	if (regs[max2_i].qe - regs[max2_i].qs < (regs[max_i].qe - regs[max_i].qs) * sqrt(frac)) return;
	div = 1. - mb_event_identity(&regs[max_i]);
	if (div < 0.02) div = 0.02;
	b2 = 0.5 / div; // max value: 25
	if (b2 * a < b) b2 = (double)a / b;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (r->p == 0) continue;
		r->p->dp_max = mb_recal_max_dp(r, b2, a);
		if (r->p->dp_max < 0) r->p->dp_max = 0;
	}
}

mb_hit_t *mb_align_skeleton(void *km, const mb_opt_t *opt, const mb_idx_t *mi, int qlen, const uint8_t *qseq, l2b_meth_t mt, int *n_regs_, mb_hit_t *regs, mb_anchor_t *a)
{
	int32_t i, n_regs = *n_regs_, n_a;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	// encode the query sequence
	qseq0[0] = Kmalloc(km, uint8_t, qlen * 2);
	qseq0[1] = qseq0[0] + qlen;
	for (i = 0; i < qlen; ++i)
		qseq0[0][i] = qseq[i], qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;

	// align through seed hits
	n_a = mb_squeeze_a(km, n_regs, regs, a);
	memset(&ez, 0, sizeof(ksw_extz_t));
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t r2; // only used for inversion
		mb_align1(km, opt, mi, qlen, qseq0, mt, &regs[i], &r2, n_a, a, &ez);
		if (r2.cnt > 0) regs = mb_insert_reg(&r2, i, &n_regs, regs);
		if (i > 0 && regs[i].split_inv) {
			if (mb_align1_inv(km, opt, mi, qlen, qseq0, mt, &regs[i-1], &regs[i], &r2, &ez)) {
				regs = mb_insert_reg(&r2, i, &n_regs, regs);
				++i; // skip the inserted INV alignment
			}
		}
	}
	kfree(km, qseq0[0]);
	kfree(km, ez.cigar);
	mb_filter_hits(opt, qlen, &n_regs, regs);
	if (!mb_is_sr_mode(opt, qlen)) {
		mb_update_dp_max(qlen, n_regs, regs, 0.9, opt->a, opt->b);
		mb_filter_hits(opt, qlen, &n_regs, regs);
	}
	mb_hit_sort(km, &n_regs, regs);
	*n_regs_ = n_regs;
	return regs;
}
