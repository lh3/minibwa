#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "hl.h"
#include "kalloc.h"
#include "bwt.h"

int mb_verbose = 3;

/**********
 * Macros *
 **********/

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
	bwt = hl_calloc(mb_bwt_t, 1);
	bwt_gen_cnt_table(bwt->cnt_table);
	return bwt;
}

void mb_bwt_destroy(mb_bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->data);
	free(bwt);
}

/******************
 * Encode raw BWT *
 ******************/

#define raw_B00(b, k) ((b)[(k)>>4]>>((~(k)&0xf)<<1)&3)

static uint64_t mb_bwt_data_len(uint64_t len)
{
	uint64_t bwt_len, occ_len;
	bwt_len = (len + 127) / 128 * 4;
	occ_len = ((len + 127) / 128 + 1) * 4; // +1 for the final counts
	return bwt_len + occ_len;
}

mb_bwt_t *mb_bwt_init_from_raw(const uint32_t *raw, uint64_t len, uint64_t primary)
{
	uint64_t c[4], x[4], i, k;
	mb_bwt_t *bwt;

	bwt = mb_bwt_init();
	bwt->primary = primary;
	bwt->seq_len = len;
	bwt->data_len = mb_bwt_data_len(len);
	bwt->data = hl_calloc(uint64_t, bwt->data_len);

	memset(c, 0, 32);
	for (i = k = 0; i < len; ++i) {
		uint8_t a = raw_B00(raw, i);
		if ((i & 0x7f) == 0) {
			if (i > 0) {
				memcpy(&bwt->data[k], x, 32);
				k += 4;
			}
			memcpy(&bwt->data[k], c, 32);
			k += 4;
			memset(x, 0, 32);
		}
		++c[a];
		x[(i&0x7f)>>5] |= (uint64_t)a << ((i&0x1f)<<1); // compatible with little endian
	}
	// the last block
	memcpy(&bwt->data[k], x, 32);
	k += 4;
	memcpy(&bwt->data[k], c, 32);
	k += 4;
	assert(k == bwt->data_len);
	for (i = 0, bwt->L2[0] = 0; i < 4; ++i)
		bwt->L2[i+1] = bwt->L2[i] + c[i];
	assert(bwt->L2[4] == len);
	return bwt;
}

/********
 * Rank *
 ********/

#define bwt_block(b, k) ((b)->data + ((k)>>7<<3))

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

/***************************
 * Suffix array operations *
 ***************************/

// retrieve a character from the $-removed BWT string. Note that mb_bwt_t::data is
// not exactly the BWT string and therefore this macro is called bwt_B0 instead of bwt_B
#define bwt_B0(b, k) (*(bwt_block(b, k) + 4 + ((k)>>5)) >> (((k)&31)<<1) & 3)

static inline uint64_t bwt_invPsi(const mb_bwt_t *bwt, uint64_t k) // compute inverse CSA
{
	uint64_t x = k - (k > bwt->primary);
	int c = bwt_B0(bwt, x);
	x = bwt->L2[c] + mb_bwt_rank11(bwt, k, c);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void mb_gen_sa(mb_bwt_t *bwt, uint32_t intv)
{
	uint64_t isa, sa, i, mask; // S(isa) = sa

	hl_assert(intv > 0 && (intv & (intv - 1)) == 0, "SA sample interval is not a power of 2.");
	hl_assert(bwt->data, "bwt_t::data is not initialized.");

	if (bwt->sa) free(bwt->sa);
	for (i = 0; 1U<<i != intv; ++i) {}
	bwt->sa_intv_bit = i;
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) >> bwt->sa_intv_bit;
	bwt->sa = hl_calloc(uint64_t, bwt->n_sa);

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

/*************************
 * Read/write BWT and SA *
 *************************/

static uint64_t read_huge(FILE *fp, uint64_t size, void *a)
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
	raw = hl_calloc(uint32_t, raw_size);
	fseek(fp, 0, SEEK_SET);
	fread(&primary, sizeof(uint64_t), 1, fp);
	fread(L2 + 1, sizeof(uint64_t), 4, fp);
	L2[0] = 0;
	read_huge(fp, raw_size<<2, raw);
	fclose(fp);
	bwt = mb_bwt_init_from_raw(raw, L2[4], primary);
	free(raw);
	return bwt;
}

int mb_bwt_save(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	uint32_t dummy = 0;

	fp = fopen(fn, "wb");
	if (fp == 0) return -1;
	fwrite(MB_MAGIC, 1, 4, fp);
	fwrite(&dummy, 4, 1, fp);
	fwrite(&bwt->primary, 8, 1, fp);
	fwrite(&bwt->L2[1], 8, 4, fp);
	fwrite(bwt->data, 8, bwt->data_len, fp);
	fclose(fp);
	return 0;
}

mb_bwt_t *mb_bwt_load(const char *fn)
{
	FILE *fp;
	char magic[4];
	uint32_t dummy;
	uint64_t x[5];
	mb_bwt_t *bwt;

	fp = fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, MB_MAGIC, 4) != 0) {
		fclose(fp);
		return 0;
	}
	fread(&dummy, 4, 1, fp);
	fread(x, 8, 5, fp);
	bwt = mb_bwt_init();
	bwt->primary = x[0];
	memcpy(&bwt->L2[1], &x[1], 32);
	bwt->seq_len = bwt->L2[4];
	bwt->data_len = mb_bwt_data_len(bwt->seq_len);
	bwt->data = hl_calloc(uint64_t, bwt->data_len);
	read_huge(fp, bwt->data_len << 3, bwt->data);
	return bwt;
}

/*
void mb_bwt_dump_sa(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	hl_assert(fp, "failed to write the suffix array sample file.");
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
	hl_assert(primary == bwt->primary, "inconsistent primary.");
	fread(skipped, sizeof(uint64_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(uint64_t), 1, fp);
	fread(&primary, sizeof(uint64_t), 1, fp);
	hl_assert(primary == bwt->seq_len, "inconsistent sequence length.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = hl_calloc(uint64_t, bwt->n_sa);
	bwt->sa[0] = -1;

	read_huge(fp, sizeof(uint64_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}
*/
