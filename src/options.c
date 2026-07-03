#include <string.h>
#include <math.h>
#include "minibwa.h"

static void mb_opt_reset(mb_opt_t *opt)
{
	memset(opt, 0, sizeof(mb_opt_t));
	opt->min_len = 19;
	opt->max_sr_len = 325;
	// seeding options
	opt->max_sub_occ = 10;
	opt->max_occ = 250;
	// chaining options
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->chain_gap_scale = 0.8f;
	opt->bw_long = 20000;
	opt->pri_ratio = 0.5f;
	// hit processing options
	opt->mask_level = 0.5f;
	opt->mask_len = 0x7fffffff;
	// alignment options
	opt->a = 2,  opt->b = 8;
	opt->q = 12, opt->q2 = 23;
	opt->e = 2,  opt->e2 = 1;
	opt->b_ambi = 1;
	// pairing options
	opt->max_pe_ins = 10000;
	opt->max_rescue = 10;
	opt->pen_unpair = 17;
	opt->pe_avg = 400, opt->pe_std = 100;
	opt->pe_lo = 50, opt->pe_hi = 800;
	// I/O options
	opt->sb_len = 1000000;
	opt->sb_seq = 24;
	opt->n_thread = 1;
	opt->seed = 11;
	opt->xa_max = 5;
	opt->xa_ratio = 0.8f;
	opt->max_sw_mat = 100000000;
	opt->cap_kalloc = 1UL<<28;
	opt->max_mb_size = 1000000000;
#ifndef HAVE_KALLOC
	opt->flag |= MB_F_NO_KALLOC;
#endif
}

void mb_opt_init(mb_opt_t *opt)
{
	mb_opt_reset(opt);
	mb_opt_preset(opt, "adap");
}

int mb_opt_preset(mb_opt_t *opt, const char *preset)
{
	mb_opt_reset(opt);
	if (strcmp(preset, "sr") == 0 || strcmp(preset, "adap") == 0) {
		opt->flag |= MB_F_PE;
		if (strcmp(preset, "adap") == 0) opt->flag |= MB_F_ADAP;
		opt->min_dp_max = 30;
		opt->flag |= MB_F_ADAP;
		opt->bw = 100;
		opt->max_gap = 100;
		opt->zdrop = 80;
		opt->zdrop_inv = 80;
		opt->best_n = 50;
		opt->end_bonus = 10;
		opt->min_chain_score = 25;
		opt->min_ksw_len = 20;
		opt->mb_size = 100000000;
	} else if (strcmp(preset, "lr") == 0) {
		opt->flag |= MB_F_LONG;
		opt->flag &= ~MB_F_PE;
		opt->min_dp_max = 50;
		opt->bw = 500;
		opt->max_gap = 5000;
		opt->zdrop = 400;
		opt->zdrop_inv = 240;
		opt->best_n = 5;
		opt->end_bonus = -1;
		opt->min_chain_score = 40;
		opt->min_ksw_len = 200;
		opt->mb_size = 500000000;
	} else {
		return -1;
	}
	return 0;
}

void mb_opt_adap(const mb_opt_t *opt0, int32_t len, mb_opt_t *opt)
{
	const int32_t min_len = 100, mid_len = 2000;
	double a, b;
	*opt = *opt0;
	if (!(opt0->flag & MB_F_ADAP)) return;
	a = -log(0.5) / (mid_len - min_len);
	b = exp(-a * ((len > min_len? len : min_len) - min_len));
	if (opt0->max_gap < 5000) {
		opt->max_gap = (int32_t)(5000 - (5000 - opt0->max_gap) * b + .499);
		if (opt->max_gap > len) opt->max_gap = len;
	}
	if (opt0->bw < 500)
		opt->bw = (int32_t)(500 - (500 - opt0->bw) * b + .499);
	if (opt->bw_long > len * 5) opt->bw_long = len * 5;
	if (opt->bw_long < opt->bw) opt->bw_long = opt->bw;
	if (opt0->zdrop < 400)
		opt->zdrop = (int32_t)(400 - (400 - opt0->zdrop) * b + .499);
	if (opt0->zdrop_inv < 240)
		opt->zdrop_inv = (int32_t)(240 - (240 - opt0->zdrop_inv) * b + .499);
	if (opt0->best_n > 5)
		opt->best_n = (int32_t)((opt0->best_n - 5) * b + 5 + .499);
	if (opt0->min_dp_max < 50)
		opt->min_dp_max = (int32_t)(50 - (50 - opt0->min_dp_max) * b + .499);
	if (opt0->min_chain_score < 40)
		opt->min_chain_score = (int32_t)(40 - (40 - opt0->min_chain_score) * b + .499);
}
