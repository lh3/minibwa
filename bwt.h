#ifndef MB_BWT_H
#define MB_BWT_H

#include <stdint.h>
#include <stddef.h>

typedef uint64_t mb_uint_t;

typedef struct {
	mb_uint_t primary; // S^{-1}(0), or the primary index of BWT
	mb_uint_t L2[5]; // C(), cumulative count
	mb_uint_t seq_len; // sequence length
	mb_uint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT; size of each block: sizeof(mb_uint_t)*4 + (128*2)/8
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int32_t sa_intv, sa_intv_bit;
	mb_uint_t n_sa;
	mb_uint_t *sa;
} mb_bwt_t;

typedef struct {
	mb_uint_t x[3], info;
} mb_intv_t;

typedef struct { size_t n, m; mb_intv_t *a; } mb_intv_v;

#define mb_bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

void mb_bwt_dump_bwt(const char *fn, const mb_bwt_t *bwt);
void mb_bwt_dump_sa(const char *fn, const mb_bwt_t *bwt);

mb_bwt_t *mb_bwt_restore_bwt(const char *fn);
void mb_bwt_restore_sa(const char *fn, mb_bwt_t *bwt);

mb_bwt_t *mb_bwt_init(void);
void mb_bwt_destroy(mb_bwt_t *bwt);

mb_uint_t mb_bwt_rank11(const mb_bwt_t *bwt, mb_uint_t k, uint8_t c);
void mb_bwt_rank1a(const mb_bwt_t *bwt, mb_uint_t k, mb_uint_t cnt[4]);

void bwt_cal_sa(mb_bwt_t *bwt, int intv);
mb_uint_t bwt_sa(const mb_bwt_t *bwt, mb_uint_t k);

#ifdef __cplusplus
}
#endif

#endif
