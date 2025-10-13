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
#include "kalloc.h"
#include "bwt.h"

/**********
 * Macros *
 **********/

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

/********************
 * Basic operations *
 ********************/

static void bwt_gen_cnt_table(uint32_t cnt[256])
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		cnt[i] = x;
	}
}

mb_bwt_t *mb_bwt_init(void)
{
	mb_bwt_t *bwt;
	bwt = mb_calloc(mb_bwt_t, 1);
	bwt_gen_cnt_table(bwt->cnt_table);
	return bwt;
}

void mb_bwt_destroy(mb_bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

/******************
 * Encode raw BWT *
 ******************/

#define raw_B00(b, k) ((b)[(k)>>4]>>((~(k)&0xf)<<1)&3)

mb_bwt_t *mb_bwt_init_from_raw(const uint32_t *raw, uint64_t len, uint64_t primary)
{
	uint64_t bwt_len, occ_len, c[4], x[4], i, k;
	mb_bwt_t *bwt;

	bwt = mb_bwt_init();
	bwt->primary = primary;
	bwt->seq_len = len;
	bwt_len = (len + 127) / 128 * 4;
	occ_len = ((len + 127) / 128 + 1) * 4; // +1 for the final counts
	bwt->bwt_size = (bwt_len + occ_len) * sizeof(uint64_t);
	bwt->bwt = mb_calloc(uint64_t, bwt_len + occ_len);

	memset(c, 0, 32);
	for (i = k = 0; i < len; ++i) {
		uint8_t a = raw_B00(raw, i);
		if ((i & 0x7f) == 0) {
			if (i > 0) {
				memcpy(&bwt->bwt[k], x, 32);
				k += 4;
			}
			memcpy(&bwt->bwt[k], c, 32);
			k += 4;
			memset(x, 0, 32);
		}
		++c[a];
		x[(i&0x7f)>>5] |= (uint64_t)a << ((i&0x1f)<<1); // compatible with little endian
	}
	// the last block
	memcpy(&bwt->bwt[k], x, 32);
	k += 4;
	memcpy(&bwt->bwt[k], c, 32);
	k += 4;
	assert(k * sizeof(uint64_t) == bwt->bwt_size);
	for (i = 0, bwt->L2[0] = 0; i < 4; ++i)
		bwt->L2[i+1] = bwt->L2[i] + c[i];
	assert(bwt->L2[4] == len);
	return bwt;
}

/********
 * Rank *
 ********/

#define bwt_block(b, k) ((b)->bwt + ((k)>>7<<3))

static inline int rank_aux1(uint64_t y, uint8_t c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

uint64_t mb_bwt_rank11(const mb_bwt_t *bwt, uint64_t k, uint8_t c)
{
	const uint64_t *p, *end;
	uint64_t n;
	if (k == 0) return 0;
	if (k == bwt->seq_len + 1) return bwt->L2[c+1] - bwt->L2[c];
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_block(bwt, k);
	n = p[c];
	p += 4;
	end = p + ((k&0x7f) >> 5);
	for (; p < end; ++p) n += rank_aux1(*p, c);
	n += rank_aux1(*p << ((~k&0x1f) << 1), c);
	if (c == 0) n -= ~k&0x1f;
	return n;
}

static inline uint32_t rank_aux4(const mb_bwt_t *bwt, uint32_t x)
{
	return bwt->cnt_table[x&0xff] + bwt->cnt_table[x>>8&0xff] + bwt->cnt_table[x>>16&0xff] + bwt->cnt_table[x>>24];
}

void mb_bwt_rank1a(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	const uint32_t *q, *end;
	const uint64_t *p;
	uint32_t x, tmp;
	if (k == 0) {
		memset(cnt, 0, 4 * sizeof(uint64_t));
		return;
	}
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_block(bwt, k); // p points to the block start
	memcpy(cnt, p, 4 * sizeof(uint64_t));
	q = (const uint32_t*)(p + 4); // 8 32-bit integers in each block
	end = q + ((k&0x7f) >> 4);
	for (x = 0; q < end; ++q) x += rank_aux4(bwt, *q); // NB: this assumes little endian
	tmp = *q << ((~k&0xf) << 1);
	x += rank_aux4(bwt, tmp) - (~k&0xf);
	cnt[0] += x&0xff, cnt[1] += x>>8&0xff, cnt[2] += x>>16&0xff, cnt[3] += x>>24;
}

void mb_bwt_rank2a(const mb_bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4])
{
	uint64_t k1 = k - 1, l1 = l - 1;
	k1 -= (k1 >= bwt->primary);
	l1 -= (l1 >= bwt->primary);
	if (k1>>7 != l1>>7 || k == 0 || l == 0) {
		mb_bwt_rank1a(bwt, k, cntk);
		mb_bwt_rank1a(bwt, l, cntl);
	} else {
		const uint64_t *p;
		const uint32_t *q, *endk, *endl;
		uint32_t x, y, tmp;
		k = k1, l = l1;
		p = bwt_block(bwt, k);
		memcpy(cntk, p, 4 * sizeof(uint64_t));
		q = (const uint32_t*)(p + 4);
		// prepare cntk[]
		endk = q + ((k&0x7f) >> 4);
		endl = q + ((l&0x7f) >> 4);
		for (x = 0; q < endk; ++q) x += rank_aux4(bwt, *q);
		y = x;
		tmp = *q << ((~k&0xf) << 1);
		x += rank_aux4(bwt, tmp) - (~k&0xf);
		// calculate cntl[] and finalize cntk[]
		for (; q < endl; ++q) y += rank_aux4(bwt, *q);
		tmp = *q << ((~l&0xf) << 1);
		y += rank_aux4(bwt, tmp) - (~l&0xf);
		memcpy(cntl, cntk, 4 * sizeof(uint64_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

/*********************
 * Bidirectional BWT *
 *********************/

void mb_bwt_extend(const mb_bwt_t *bwt, const mb_sai_t *ik, mb_sai_t ok[4], int is_back)
{
	uint64_t tk[4], tl[4];
	int i;
	mb_bwt_rank2a(bwt, ik->x[!is_back], ik->x[!is_back] + ik->size, tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i]; // +1 for the missing sentinel
		ok[i].size = (tl[i] -= tk[i]);
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->size - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[0].x[is_back] = ok[1].x[is_back] + tl[1];
}

int64_t mb_bwt_smem1(void *km, const mb_bwt_t *f, int64_t min_occ, int64_t min_len, int64_t len, const uint8_t *q, int64_t x, mb_sai_v *mem, int32_t check_long)
{
	int64_t i, j;
	mb_sai_t ik, ok[4], *p;

	assert(len <= INT32_MAX); // this can be relaxed if we define a new struct for mem
	if (len - x < min_len) return len;
	mb_bwt_set_intv(f, q[x + min_len - 1], &ik);
	for (i = x + min_len - 2; i >= x; --i) { // backward extension
		int c = q[i];
		mb_bwt_extend(f, &ik, ok, 1);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	if (i >= x) return i + 1; // no MEM found
	if (check_long) return -1;
	for (j = x + min_len; j < len; ++j) { // forward extension
		int c = 3 - q[j];
		mb_bwt_extend(f, &ik, ok, 0);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	Kgrow(km, mb_sai_t, mem->a, mem->n, mem->m);
	p = &mem->a[mem->n++];
	*p = ik;
	p->info = (uint64_t)x<<32 | j;
	if (j == len) return len;
	mb_bwt_set_intv(f, q[j], &ik);
	for (i = j - 1; i > x; --i) { // backward extension again
		int c = q[i];
		mb_bwt_extend(f, &ik, ok, 1);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	return i + 1;
}

int64_t mb_bwt_smem(void *km, const mb_bwt_t *f, int64_t len, const uint8_t *q, mb_sai_v *mem, int64_t min_occ, int64_t min_len)
{
	int64_t x = 0;
	mem->n = 0;
	do {
		x = mb_bwt_smem1(km, f, min_occ, min_len, len, q, x, mem, 0);
	} while (x < len);
	return mem->n;
}

/*
// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, uint64_t k, uint64_t l, uint8_t c, uint64_t *ok, uint64_t *ol)
{
	uint64_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (uint64_t)(-1) || l == (uint64_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		uint64_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		n = ((uint64_t*)(p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(uint64_t);
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

int bwt_match_exact(const bwt_t *bwt, int len, const uint8_t *str, uint64_t *sa_begin, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
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

int bwt_match_exact_alt(const bwt_t *bwt, int len, const uint8_t *str, uint64_t *k0, uint64_t *l0)
{
	int i;
	uint64_t k, l, ok, ol;
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
/*
static inline uint64_t bwt_invPsi(const mb_bwt_t *bwt, uint64_t k) // compute inverse CSA
{
	uint64_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + mb_bwt_rank11(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void mb_gen_sa(mb_bwt_t *bwt, uint32_t intv)
{
	uint64_t isa, sa, i, mask; // S(isa) = sa

	mb_assert(intv > 0 && (intv & (intv - 1)) == 0, "SA sample interval is not a power of 2.");
	mb_assert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	for (i = 0; 1U<<i != intv; ++i) {}
	bwt->sa_intv_bit = i;
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) >> bwt->sa_intv_bit;
	bwt->sa = mb_calloc(uint64_t, bwt->n_sa);

	// calculate SA value
	isa = 0, sa = bwt->seq_len, mask = intv - 1;
	for (i = 0; i < bwt->seq_len; ++i) {
		if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_intv_bit] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_intv_bit] = sa;
	bwt->sa[0] = (uint64_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

uint64_t mb_bwt_sa(const mb_bwt_t *bwt, uint64_t k)
{
	uint64_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	// without setting bwt->sa[0] = -1, the following line should be
	// changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1)
	return sa + bwt->sa[k >> bwt->sa_intv_bit];
}
*/
/*************************
 * Read/write BWT and SA *
 *************************/

static uint64_t fread_huge(FILE *fp, uint64_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	uint64_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = fread(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

mb_bwt_t *mb_bwt_load_raw(const char *fn)
{
	mb_bwt_t *bwt;
	uint32_t *raw;
	uint64_t L2[5], primary, raw_size;
	FILE *fp;

	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	raw_size = (ftell(fp) - sizeof(uint64_t) * 5) >> 2;
	raw = mb_calloc(uint32_t, raw_size);
	fseek(fp, 0, SEEK_SET);
	fread(&primary, sizeof(uint64_t), 1, fp);
	fread(L2 + 1, sizeof(uint64_t), 4, fp);
	L2[0] = 0;
	fread_huge(fp, raw_size<<2, raw);
	fclose(fp);
	bwt = mb_bwt_init_from_raw(raw, L2[4], primary);
	free(raw);
	return bwt;
}
/*
void mb_bwt_dump_sa(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	mb_assert(fp, "failed to write the suffix array sample file.");
	fwrite(&bwt->primary, sizeof(uint64_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(uint64_t), 4, fp);
	fwrite(&bwt->sa_intv, sizeof(uint64_t), 1, fp);
	fwrite(&bwt->seq_len, sizeof(uint64_t), 1, fp);
	fwrite(bwt->sa + 1, sizeof(uint64_t), bwt->n_sa - 1, fp);
	fflush(fp);
	fclose(fp);
}

void mb_bwt_restore_sa(const char *fn, mb_bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	uint64_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(uint64_t), 1, fp);
	mb_assert(primary == bwt->primary, "inconsistent primary.");
	fread(skipped, sizeof(uint64_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(uint64_t), 1, fp);
	fread(&primary, sizeof(uint64_t), 1, fp);
	mb_assert(primary == bwt->seq_len, "inconsistent sequence length.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = mb_calloc(uint64_t, bwt->n_sa);
	bwt->sa[0] = -1;

	fread_huge(fp, sizeof(uint64_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}
*/
