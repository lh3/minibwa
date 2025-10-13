/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "utils.h"
#include "bwt.h"

/**********
 * Macros *
 **********/

// requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#if OCC_INTV_SHIFT == 7
#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(mb_uint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))
#else
#define bwt_bwt(b, k) ((b)->bwt[(k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(mb_uint_t)/4*4) + sizeof(mb_uint_t)/4*4 + (k)%OCC_INTERVAL/16])
#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(mb_uint_t)/4*4)
#error "OCC_INTV_SHIFT MUST BE 7" // although bwt_bwt() and bwt_occ_intv() are correct, other places may be wrong
#endif

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

/********************
 * Basic operations *
 ********************/

static inline void bwt_gen_cnt_table(mb_bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

mb_bwt_t *mb_bwt_init(void)
{
	mb_bwt_t *bwt;
	bwt = mb_calloc(mb_bwt_t, 1);
	bwt_gen_cnt_table(bwt);
	return bwt;
}

void mb_bwt_destroy(mb_bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

/********
 * Rank *
 ********/

static inline int rank_aux1(uint64_t y, uint8_t c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

mb_uint_t mb_bwt_rank11(const mb_bwt_t *bwt, mb_uint_t k0, uint8_t c)
{
	mb_uint_t n, k;
	uint32_t *p, *end;

	if (k0 == 0) return 0;
	if (k0 == bwt->seq_len + 1) return bwt->L2[c+1] - bwt->L2[c];
	k = k0 - 1;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	p = bwt_occ_intv(bwt, k);
	n = ((mb_uint_t*)p)[c];
	p += sizeof(mb_uint_t); // jump to the start of the first BWT cell; sizeof(mb_uint_t) = sizeof(mb_uint_t) * 4 / sizeof(uint32_t)

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += rank_aux1((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += rank_aux1(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

static inline uint32_t rank_aux4(const mb_bwt_t *bwt, uint8_t c)
{
	return bwt->cnt_table[c&0xff] + bwt->cnt_table[c>>8&0xff] + bwt->cnt_table[c>>16&0xff] + bwt->cnt_table[c>>24];
}

void mb_bwt_rank1a(const mb_bwt_t *bwt, mb_uint_t k0, mb_uint_t cnt[4])
{
	mb_uint_t x, k;
	uint32_t *p, tmp, *end;

	memset(cnt, 0, 4 * sizeof(mb_uint_t));
	if (k0 == 0) return;
	k = k0 - 1;
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(mb_uint_t));
	p += sizeof(mb_uint_t); // sizeof(mb_uint_t) = 4*(sizeof(mb_uint_t)/sizeof(uint32_t))
	end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p) x += rank_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += rank_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

/*
// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, mb_uint_t k, mb_uint_t l, uint8_t c, mb_uint_t *ok, mb_uint_t *ol)
{
	mb_uint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (mb_uint_t)(-1) || l == (mb_uint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		mb_uint_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		n = ((mb_uint_t*)(p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(mb_uint_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
			n += rank_aux((uint64_t)p[0]<<32 | p[1], c);
		m = n;
		n += rank_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += rank_aux((uint64_t)p[0]<<32 | p[1], c);
		m += rank_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
		if (c == 0) m -= ~l&31; // corrected for the masked bits
		*ol = m;
	}
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
void bwt_2occ4(const bwt_t *bwt, mb_uint_t k, mb_uint_t l, mb_uint_t cntk[4], mb_uint_t cntl[4])
{
	mb_uint_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (mb_uint_t)(-1) || l == (mb_uint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		mb_uint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(mb_uint_t));
		p += sizeof(mb_uint_t); // sizeof(mb_uint_t) = 4*(sizeof(mb_uint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
		for (x = 0; p < endk; ++p) x += rank_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += rank_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += rank_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += rank_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(mb_uint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

int bwt_match_exact(const bwt_t *bwt, int len, const uint8_t *str, mb_uint_t *sa_begin, mb_uint_t *sa_end)
{
	mb_uint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		uint8_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const uint8_t *str, mb_uint_t *k0, mb_uint_t *l0)
{
	int i;
	mb_uint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		uint8_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}
*/
/*********************
 * Bidirectional BWT *
 *********************/
/*
void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	mb_uint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}
// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			if (c >= 0 && ik.x[2] >= max_intv) bwt_extend(bwt, p, ok, 1);
			if (c < 0 || ik.x[2] < max_intv || ok[c].x[2] < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	return bwt_smem1a(bwt, len, q, x, min_intv, 0, mem, tmpvec);
}

int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem)
{
	int i, c;
	bwtintv_t ik, ok[4];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] < max_intv && i - x >= min_len) {
				*mem = ok[c];
				mem->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	return len;
}
*/
/***************************
 * Suffix array operations *
 ***************************/

static inline mb_uint_t bwt_invPsi(const mb_bwt_t *bwt, mb_uint_t k) // compute inverse CSA
{
	mb_uint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + mb_bwt_rank11(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void mb_gen_sa(mb_bwt_t *bwt, uint32_t intv)
{
	mb_uint_t isa, sa, i, mask; // S(isa) = sa

	mb_assert(intv > 0 && (intv & (intv - 1)) == 0, "SA sample interval is not a power of 2.");
	mb_assert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	for (i = 0; 1U<<i != intv; ++i) {}
	bwt->sa_intv_bit = i;
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) >> bwt->sa_intv_bit;
	bwt->sa = mb_calloc(mb_uint_t, bwt->n_sa);

	// calculate SA value
	isa = 0, sa = bwt->seq_len, mask = intv - 1;
	for (i = 0; i < bwt->seq_len; ++i) {
		if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_intv_bit] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_intv_bit] = sa;
	bwt->sa[0] = (mb_uint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

mb_uint_t mb_bwt_sa(const mb_bwt_t *bwt, mb_uint_t k)
{
	mb_uint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k >> bwt->sa_intv_bit];
}

/*************************
 * Read/write BWT and SA *
 *************************/

void mb_bwt_dump_bwt(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	mb_assert(fp, "failed to write the BWT file.");
	fwrite(&bwt->primary, sizeof(mb_uint_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(mb_uint_t), 4, fp);
	fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	fflush(fp);
	fclose(fp);
}

void mb_bwt_dump_sa(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	mb_assert(fp, "failed to write the suffix array sample file.");
	fwrite(&bwt->primary, sizeof(mb_uint_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(mb_uint_t), 4, fp);
	fwrite(&bwt->sa_intv, sizeof(mb_uint_t), 1, fp);
	fwrite(&bwt->seq_len, sizeof(mb_uint_t), 1, fp);
	fwrite(bwt->sa + 1, sizeof(mb_uint_t), bwt->n_sa - 1, fp);
	fflush(fp);
	fclose(fp);
}

static mb_uint_t fread_huge(FILE *fp, mb_uint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	mb_uint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = fread(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void mb_bwt_restore_sa(const char *fn, mb_bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	mb_uint_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(mb_uint_t), 1, fp);
	mb_assert(primary == bwt->primary, "inconsistent primary.");
	fread(skipped, sizeof(mb_uint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(mb_uint_t), 1, fp);
	fread(&primary, sizeof(mb_uint_t), 1, fp);
	mb_assert(primary == bwt->seq_len, "inconsistent sequence length.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = mb_calloc(mb_uint_t, bwt->n_sa);
	bwt->sa[0] = -1;

	fread_huge(fp, sizeof(mb_uint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}

mb_bwt_t *mb_bwt_restore_bwt(const char *fn)
{
	mb_bwt_t *bwt;
	FILE *fp;

	bwt = mb_bwt_init();
	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(mb_uint_t) * 5) >> 2;
	bwt->bwt = mb_calloc(uint32_t, bwt->bwt_size);
	fseek(fp, 0, SEEK_SET);
	fread(&bwt->primary, sizeof(mb_uint_t), 1, fp);
	fread(bwt->L2 + 1, sizeof(mb_uint_t), 4, fp);
	fread_huge(fp, bwt->bwt_size<<2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	return bwt;
}
