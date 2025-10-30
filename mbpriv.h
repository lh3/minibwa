#ifndef MBPRIV_H
#define MBPRIV_H

#include "minibwa.h"
#include "l2bit.h"
#include "bwt.h"

#define MB_DBG_NO_KALLOC    (0x1LL)

struct mb_idx_s {
	l2b_t *l2b;
	mb_bwt_t *bwt;
};

#ifdef __cplusplus
extern "C" {
#endif

// defined in bwtgen.c
void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

#ifdef __cplusplus
}
#endif

#endif
