#ifndef MB_BWT_H
#define MB_BWT_H

#include <stdint.h>
#include <stddef.h>

#define MB_MAGIC "MBW\1"

typedef struct {
	uint64_t primary; // S^{-1}(0), or the primary index of BWT
	uint64_t L2[5]; // C(), cumulative count
	uint64_t seq_len; // sequence length
	uint64_t data_len;
	uint64_t *data; // BWT
	uint32_t cnt_table[256];
	// suffix array
	uint32_t sa_bit; // sample rate: 1/(1<<sa_bit)
	uint64_t n_sa;
	uint64_t *sa;
} mb_bwt_t;

typedef struct {
	uint64_t x[2];
	uint64_t size, info;
} mb_sai_t;

typedef struct { size_t n, m; mb_sai_t *a; } mb_sai_v;

#ifdef __cplusplus
extern "C" {
#endif

extern int mb_verbose;

mb_bwt_t *mb_bwt_init(void);
void mb_bwt_destroy(mb_bwt_t *bwt);
mb_bwt_t *mb_bwt_load_raw(const char *fn); // from raw bwt_gen.c output
int mb_bwt_save(const char *fn, const mb_bwt_t *bwt);
mb_bwt_t *mb_bwt_load(const char *fn);
mb_bwt_t *mb_bwt_init_from_raw(int is_byte, const void *raw_, uint64_t len, uint64_t primary);

uint64_t mb_bwt_rank11(const mb_bwt_t *bwt, uint64_t k, uint8_t c);
void mb_bwt_rank1a(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4]);
void mb_bwt_rank2a(const mb_bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4]);

int64_t mb_bwt_smem(void *km, const mb_bwt_t *f, int64_t min_len, int64_t min_occ, int64_t len, const uint8_t *q, int64_t x, mb_sai_t *p);

void mb_bwt_gen_sa(mb_bwt_t *bwt, uint32_t sa_bit);
uint64_t mb_bwt_sa(const mb_bwt_t *bwt, uint64_t k);
void mb_bwt_sa2(void *km, const mb_bwt_t *bwt, int64_t n, uint64_t *x);

static inline void mb_bwt_set_intv(const mb_bwt_t *bwt, int c, mb_sai_t *ik)
{
	ik->x[0] = bwt->L2[c] + 1; // +1 for the missing sentinel
	ik->x[1] = bwt->L2[3-c] + 1;
	ik->size = bwt->L2[c+1] - bwt->L2[c];
	ik->info = 0;
}

#ifdef __cplusplus
}
#endif

#endif
