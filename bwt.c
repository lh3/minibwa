#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "kommon.h"
#include "kalloc.h"
#include "bwt.h"

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
	bwt = kom_calloc(mb_bwt_t, 1);
	bwt->sa_bit = (uint32_t)-1;
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

#define BWT_CNT_SHIFT 56
#define BWT_CNT_MASK ((1ULL<<BWT_CNT_SHIFT) - 1)

static uint64_t mb_bwt_data_len(uint64_t len)
{
	uint64_t bwt_len, occ_len;
	bwt_len = (len + 127) / 128 * 4;
	occ_len = ((len + 127) / 128 + 1) * 4; // +1 for the final counts
	return bwt_len + occ_len;
}

mb_bwt_t *mb_bwt_init_from_raw(int is_byte, const void *raw_, uint64_t len, uint64_t primary)
{
	uint64_t c[4], x[4], i, k, *last_c = 0;
	mb_bwt_t *bwt;
	const uint32_t *raw32 = is_byte? 0 : (const uint32_t*)raw_;
	const uint8_t *raw8 = is_byte? (const uint8_t*)raw_ : 0;

	bwt = mb_bwt_init();
	bwt->primary = primary;
	bwt->seq_len = len;
	bwt->data_len = mb_bwt_data_len(len);
	bwt->data = kom_calloc(uint64_t, bwt->data_len);

	memset(c, 0, 32);
	for (i = k = 0; i < len; ++i) {
		uint8_t j, a = is_byte? raw8[i]&3 : raw_B00(raw32, i);
		if ((i & 0x7f) == 0) { // hard coded: each block encodes 128bp
			if (i > 0) {
				memcpy(&bwt->data[k], x, 32);
				k += 4;
			}
			last_c = &bwt->data[k];
			memcpy(&bwt->data[k], c, 32);
			k += 4;
			memset(x, 0, 32);
		} else if ((i & 0x3f) == 0 && last_c) {
			for (j = 0; j < 4; ++j)
				last_c[j] |= (c[j] - last_c[j]) << BWT_CNT_SHIFT;
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

static inline void mb_bwt_block_prefetch(const mb_bwt_t *bwt, uint64_t k)
{
	if (k > 0) __builtin_prefetch(bwt_block(bwt, k - 1 - (k - 1 >= bwt->primary)));
}

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
	uint64_t n, mask;
	if (k == 0) return 0;
	if (k == bwt->seq_len + 1) return bwt->L2[c+1] - bwt->L2[c];
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	mask = (k&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	p = bwt_block(bwt, k);
	n = (p[c] & BWT_CNT_MASK) + ((p[c] >> BWT_CNT_SHIFT) & mask);
	p += 4; // p points to 2-bit encoded BWT
	end = p + ((k&0x7f) >> 5);
	p += (k&0x7f) >= 64? 2 : 0;
	for (; p < end; ++p) n += rank_aux1(*p, c);
	n += rank_aux1(*p << ((~k&0x1f) << 1), c);
	if (c == 0) n -= ~k&0x1f;
	return n;
}

static inline const uint32_t *seek_block(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	const uint64_t *p;
	uint64_t mask = (k&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	p = bwt_block(bwt, k);
	cnt[0] = (p[0] & BWT_CNT_MASK) + ((p[0] >> BWT_CNT_SHIFT) & mask);
	cnt[1] = (p[1] & BWT_CNT_MASK) + ((p[1] >> BWT_CNT_SHIFT) & mask);
	cnt[2] = (p[2] & BWT_CNT_MASK) + ((p[2] >> BWT_CNT_SHIFT) & mask);
	cnt[3] = (p[3] & BWT_CNT_MASK) + ((p[3] >> BWT_CNT_SHIFT) & mask);
	return (const uint32_t*)(p + 4);
}

static inline uint32_t rank_aux4(const mb_bwt_t *bwt, uint32_t x)
{
	return bwt->cnt_table[x&0xff] + bwt->cnt_table[x>>8&0xff] + bwt->cnt_table[x>>16&0xff] + bwt->cnt_table[x>>24];
}

void mb_bwt_rank1a(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	const uint32_t *q, *end;
	uint32_t x, tmp;
	if (k == 0) {
		memset(cnt, 0, 4 * sizeof(uint64_t));
		return;
	}
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	q = seek_block(bwt, k, cnt);
	end = q + ((k&0x7f) >> 4);
	q += (k&0x7f) >= 64? 4 : 0;
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
	mb_bwt_block_prefetch(bwt, k);
	if (k1>>7 != l1>>7 || k == 0 || l == 0) {
		mb_bwt_block_prefetch(bwt, l);
		mb_bwt_rank1a(bwt, k, cntk);
		mb_bwt_rank1a(bwt, l, cntl);
	} else {
		const uint32_t *q, *endk, *endl;
		uint32_t x, y, tmp;
		k = k1, l = l1;
		q = seek_block(bwt, k, cntk);
		// prepare cntk[]
		endk = q + ((k&0x7f) >> 4);
		endl = q + ((l&0x7f) >> 4);
		q += (k&0x7f) >= 64? 4 : 0;
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
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->size > bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[0].x[is_back] = ok[1].x[is_back] + tl[1];
}

int64_t mb_bwt_smem(const mb_bwt_t *f, uint32_t len, const uint8_t *q, int64_t x, int64_t min_len, int64_t min_occ, int64_t max_occ, mb_sai_t *p)
{
	int64_t i, j, xn, min_occ2;
	mb_sai_t ik, ok[4];

	assert(min_occ <= max_occ);
	assert(len <= INT32_MAX); // this can be relaxed if we define a new struct for mem
	p->size = 0;
	if (len - x < min_len) return len;
	for (i = x, xn = -1; i < x + min_len; ++i) // find the last N in [x,x+min_len)
		if (q[i] > 3) xn = i;
	if (xn >= 0) return xn + 1;
	mb_bwt_set_intv(f, q[x + min_len - 1], &ik);
	for (i = x + min_len - 2; i >= x; --i) { // backward extension
		int c = q[i];
		if (c > 3) break;
		mb_bwt_extend(f, &ik, ok, 1);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	if (i >= x) return i + 1; // no MEM found
	min_occ2 = ik.size < max_occ? ik.size : max_occ; // min_occ2>=min_occ; if min_occ==max_occ, min_occ2=min_occ
	for (j = x + min_len; j < len; ++j) { // forward extension
		int c = 3 - q[j];
		if (q[j] > 3) break;
		mb_bwt_extend(f, &ik, ok, 0);
		if (ok[c].size < min_occ2) break;
		ik = ok[c];
	}
	*p = ik;
	p->info = (uint64_t)x<<32 | j;
	if (j == len) return len;
	mb_bwt_set_intv(f, q[j], &ik);
	for (i = j - 1; i > x; --i) { // backward extension again
		int c = q[i];
		mb_bwt_extend(f, &ik, ok, 1);
		if (ok[c].size < min_occ2) break;
		ik = ok[c];
	}
	return i + 1;
}

/***************************
 * Suffix array operations *
 ***************************/

// retrieve a character from the $-removed BWT string. Note that mb_bwt_t::data is
// not exactly the BWT string and therefore this macro is called bwt_B0 instead of bwt_B
#define bwt_B0(b, k) ((b)->data[((k)>>7<<3) + 4 + (((k)&127)>>5)] >> (((k)&31)<<1) & 3)

static inline uint64_t bwt_invPsi(const mb_bwt_t *bwt, uint64_t k) // compute inverse CSA
{
	uint64_t x = k - (k > bwt->primary);
	int c = bwt_B0(bwt, x);
	x = bwt->L2[c] + 1 + mb_bwt_rank11(bwt, k, c); // +1 to account for the sentinel
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void mb_bwt_gen_sa(mb_bwt_t *bwt, uint32_t sa_bit)
{
	uint64_t isa, sa, i, mask; // S(isa) = sa

	assert(bwt->data);
	if (bwt->sa) free(bwt->sa);
	bwt->sa_bit = sa_bit;
	bwt->n_sa = (bwt->seq_len + (1<<sa_bit)) >> sa_bit;
	bwt->sa = kom_calloc(uint64_t, bwt->n_sa);

	// calculate SA value
	isa = 0, sa = bwt->seq_len, mask = (1ULL<<sa_bit) - 1;
	for (i = 0; i < bwt->seq_len; ++i) {
		if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_bit] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_bit] = sa;
	bwt->sa[0] = (uint64_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

uint64_t mb_bwt_sa(const mb_bwt_t *bwt, uint64_t k)
{
	uint64_t sa = 0, mask = (1ULL<<bwt->sa_bit) - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	// without setting bwt->sa[0] = -1, the following line should be
	// changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1)
	return sa + bwt->sa[k >> bwt->sa_bit];
}

void mb_bwt_sa_batch(void *km, const mb_bwt_t *bwt, int64_t n, uint64_t *x)
{
	uint64_t mask = (1ULL<<bwt->sa_bit) - 1;
	int64_t i, step = 0, r = n;
	kom128_t *z;
	z = Kmalloc(km, kom128_t, n);
	for (i = 0; i < n; ++i)
		z[i].x = x[i], z[i].y = i;
	for (step = 0; r; ++step) {
		int64_t r0;
		for (i = 0; i < r; ++i) {
			if (i + 1 < r) mb_bwt_block_prefetch(bwt, z[i+1].x);
			z[i].x = bwt_invPsi(bwt, z[i].x);
		}
		for (i = 0; i < r; ++i)
			if ((z[i].x & mask) == 0)
				x[z[i].y] = step + 1 + bwt->sa[z[i].x >> bwt->sa_bit];
		r0 = r;
		for (i = 0, r = 0; i < r0; ++i)
			if (z[i].x & mask)
				z[r++] = z[i];
	}
	kfree(km, z);
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
	raw = kom_calloc(uint32_t, raw_size);
	fseek(fp, 0, SEEK_SET);
	fread(&primary, sizeof(uint64_t), 1, fp);
	fread(L2 + 1, sizeof(uint64_t), 4, fp);
	L2[0] = 0;
	read_huge(fp, raw_size<<2, raw);
	fclose(fp);
	bwt = mb_bwt_init_from_raw(0, raw, L2[4], primary);
	free(raw);
	return bwt;
}

int mb_bwt_save(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	if (fp == 0) return -1;
	fwrite(MB_MAGIC, 1, 4, fp);
	fwrite(&bwt->sa_bit, 4, 1, fp);
	fwrite(&bwt->primary, 8, 1, fp);
	fwrite(&bwt->L2[1], 8, 4, fp);
	fwrite(bwt->data, 8, bwt->data_len, fp);
	fwrite(&bwt->n_sa, 8, 1, fp);
	if (bwt->sa_bit != (uint32_t)-1 && bwt->n_sa > 0 && bwt->sa)
		fwrite(bwt->sa, 8, bwt->n_sa, fp);
	fclose(fp);
	return 0;
}

mb_bwt_t *mb_bwt_load(const char *fn)
{
	FILE *fp;
	char magic[4];
	uint64_t x[5];
	mb_bwt_t *bwt;

	fp = fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, MB_MAGIC, 4) != 0) {
		fclose(fp);
		return 0;
	}
	bwt = mb_bwt_init();
	fread(&bwt->sa_bit, 4, 1, fp);
	fread(x, 8, 5, fp);
	bwt->primary = x[0];
	memcpy(&bwt->L2[1], &x[1], 32);
	bwt->seq_len = bwt->L2[4];
	bwt->data_len = mb_bwt_data_len(bwt->seq_len);
	bwt->data = kom_calloc(uint64_t, bwt->data_len);
	read_huge(fp, bwt->data_len << 3, bwt->data);
	fread(&bwt->n_sa, 8, 1, fp);
	if (bwt->sa_bit != (uint32_t)-1 && bwt->n_sa > 0) {
		bwt->sa = kom_malloc(uint64_t, bwt->n_sa);
		fread(bwt->sa, 8, bwt->n_sa, fp);
	}
	fclose(fp);
	return bwt;
}
