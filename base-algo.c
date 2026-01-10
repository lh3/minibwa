#include "bwt.h"
#include "minibwa.h"
#include "kalloc.h"

mb_sai_t *mb_seed_intv(void *km, const mb_seedopt_t *opt, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int64_t *n_a_)
{
	int64_t x = 0;
	int64_t n_a = 0, m_a = 0, i, n_a0;
	mb_sai_t p, *a = 0;

	do {
		x = mb_bwt_smem(bwt, opt->min_len, 1, 1, len, seq, x, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, a, n_a, m_a);
			a[n_a++] = p;
		}
	} while (x < len);

	n_a0 = n_a;
	for (i = 0; i < n_a0; ++i) {
		mb_sai_t *q = &a[i];
		uint32_t st = q->info>>32, en = (uint32_t)q->info;
		if (en - st < opt->min_len * 2 || q->size > opt->max_sub_occ)
			continue;
		x = 0;
		do {
			x = mb_bwt_smem(bwt, opt->min_len, q->size + 1, opt->max_sub_occ * 2, len, seq, x, &p);
			if (p.size > 0) {
				Kgrow(km, mb_sai_t, a, n_a, m_a);
				a[n_a++] = p;
			}
		} while (x < len);
	}
	*n_a_ = n_a;
	return a;
}
