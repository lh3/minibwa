#include <stdlib.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"

struct mb_tbuf_s {
	void *km;
};

mb_tbuf_t *mb_tbuf_init(void)
{
	mb_tbuf_t *b;
	b = kom_calloc(mb_tbuf_t, 1);
	if (!(kom_dbg_flag & MB_DBG_NO_KALLOC))
		b->km = km_init();
	return b;
}

void mb_tbuf_destroy(mb_tbuf_t *b)
{
	if (b->km) km_destroy(b->km);
	free(b);
}

static void mb_seed_intv(void *km, const mb_idx_t *idx, const mb_mopt_t *opt, int qlen, const uint8_t *seq)
{
}

mb_hit_t *mb_map(const mb_idx_t *idx, int64_t qlen, const char *seq0, int32_t *n_hit, mb_tbuf_t *b, const mb_mopt_t *opt, const char *qname)
{
	int64_t i;
	void *km = b->km;
	uint8_t *seq;

	seq = Kcalloc(km, uint8_t, qlen);
	for (i = 0; i < qlen; ++i)
		seq[i] = kom_nt4_table[(uint8_t)seq0[i]];
	free(seq);
}
