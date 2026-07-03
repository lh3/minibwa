#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "ksw2.h"

#if defined(__ARM_NEON)
#include "s2n-lite.h"
#elif defined(__SSE2__)
#include <xmmintrin.h>
#else
#error "Missing SSE2 or NEON intrinsics"
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#ifndef Kgrow
#define Kgrow(km, type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = (type*)krealloc(km, ptr, (__m) * sizeof(type)); \
		} \
	} while (0)
#endif

typedef struct {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
	void *km;
} kswq_t;

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
	kswq_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
	q = (kswq_t*)kmalloc(km, sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	q->km = km;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? -1 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? -1 : ma[query[k]]);
		}
	}
	return q;
}

static inline int ksw_le_u8(__m128i a, __m128i b)
{
#if defined(__ARM_NEON)
	return vmaxvq_u8(vqsubq_u8(a, b)) == 0;
#elif defined(__SSE2__)
	return _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(a, b), _mm_setzero_si128())) == 0xffff;
#endif
}

static inline int ksw_max_u8(__m128i x)
{
#if defined(__ARM_NEON)
	return vmaxvq_u8(x);
#elif defined(__SSE2__)
	x = _mm_max_epu8(x, _mm_srli_si128(x, 8));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 4));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 2));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 1));
	return _mm_extract_epi16(x, 0) & 0x00ff;
#endif
}

ksw_llrst_t ksw_ll_u8_core(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int xtra)
{
	kswq_t *q = (kswq_t*)q_;
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i gapoe, gape, shift, *H0, *H1, *E, *Hmax;
	ksw_llrst_t r = { 0, -1, -1, -1, -1 };

	// initialization
	minsc = (xtra&KSW_LL_SUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_LL_STOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	gapoe = _mm_set1_epi8(_gapo + _gape);
	gape = _mm_set1_epi8(_gape);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	memset(E,    0, slen * sizeof(__m128i));
	memset(H0,   0, slen * sizeof(__m128i));
	memset(Hmax, 0, slen * sizeof(__m128i));
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, h, t, f, max, *S = q->qp + target[i] * slen; // s is the 1st score vector
		f = max = _mm_setzero_si128();
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, gape); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, gapoe); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, gape);
			t = _mm_subs_epu8(h, gapoe); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, gapoe);
				f = _mm_subs_epu8(f, gape);
				if (UNLIKELY(ksw_le_u8(f, h))) goto end_loop_u8;
			}
		}
end_loop_u8:
		imax = ksw_max_u8(max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) Kgrow(q->km, uint64_t, b, n_b, m_b);
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*)Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	kfree(q->km, b);
	return r;
}

static inline int ksw_le_epi16(__m128i a, __m128i b)
{
#if defined(__ARM_NEON)
	return vmaxvq_u16(vcgtq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))) == 0;
#elif defined(__SSE2__)
	return _mm_movemask_epi8(_mm_cmpgt_epi16(a, b)) == 0;
#endif
}

static inline int ksw_max_i16(__m128i x)
{
#if defined(__ARM_NEON)
	return vmaxvq_s16(vreinterpretq_s16_u8(x));
#elif defined(__SSE2__)
	x = _mm_max_epi16(x, _mm_srli_si128(x, 8));
	x = _mm_max_epi16(x, _mm_srli_si128(x, 4));
	x = _mm_max_epi16(x, _mm_srli_si128(x, 2));
	return _mm_extract_epi16(x, 0);
#endif
}

ksw_llrst_t ksw_ll_i16_core(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int xtra)
{
	kswq_t *q = (kswq_t*)q_;
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i gapoe, gape, *H0, *H1, *E, *Hmax;
	ksw_llrst_t r = { 0, -1, -1, -1, -1 };

	// initialization
	minsc = (xtra&KSW_LL_SUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_LL_STOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	gapoe = _mm_set1_epi16(_gapo + _gape);
	gape = _mm_set1_epi16(_gape);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	memset(E,    0, slen * sizeof(__m128i));
	memset(H0,   0, slen * sizeof(__m128i));
	memset(Hmax, 0, slen * sizeof(__m128i));
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f, max, *S = q->qp + target[i] * slen; // s is the 1st score vector
		f = max = _mm_setzero_si128();
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, _mm_load_si128(S++));
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, gape);
			t = _mm_subs_epu16(h, gapoe);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, gape);
			t = _mm_subs_epu16(h, gapoe);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, gapoe);
				f = _mm_subs_epu16(f, gape);
				if (UNLIKELY(ksw_le_epi16(f, h))) goto end_loop_i16;
			}
		}
end_loop_i16:
		imax = ksw_max_i16(max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) Kgrow(q->km, uint64_t, b, n_b, m_b);
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	r.score = gmax; r.te = te;
	{
		int max = -1, tmp, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp;
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	kfree(q->km, b);
	return r;
}

int ksw_ll_i16(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int *qe, int *te)
{
	ksw_llrst_t r;
	r = ksw_ll_i16_core(q_, tlen, target, _gapo, _gape, 0);
	*qe = r.qe, *te = r.te;
	return r.score;
}
