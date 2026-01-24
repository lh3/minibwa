#include <string.h>
#include "minibwa.h"

void mb_opt_init(mb_opt_t *opt)
{
	memset(opt, 0, sizeof(mb_opt_t));
	// seeding options
	opt->max_sub_occ = 10;
	opt->max_occ = 250;
	// chaining options
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->chain_gap_scale = 0.8f;
	// hit processing options
	opt->mask_level = 0.5f;
	opt->mask_len = 0x7fffffff;
	// alignment options
	opt->a = 2,  opt->b = 7;
	opt->q = 13, opt->q2 = 34;
	opt->e = 3,  opt->e2 = 1;
	opt->b_ambi = 1;
	opt->zdrop = 400;
	opt->zdrop_inv = 200;
	// pairing options
	opt->max_pe_ins = 10000;
	// I/O options
	opt->sb_len = 1000000;
	opt->sb_seq = 24;
	opt->n_thread = 1;
	opt->seed = 11;
	opt->max_sw_mat = 100000000;
	// default to sr
	mb_opt_preset(opt, "sr");
}

int mb_opt_preset(mb_opt_t *opt, const char *preset)
{
	if (strcmp(preset, "sr") == 0) {
		opt->flag |= MB_F_PE;
		opt->min_len = 19;
		opt->bw = opt->bw_long = 150;
		opt->max_gap = 150;
		opt->pri_ratio = 0.5f;
		opt->best_n = 101;
		opt->end_bonus = 10;
		opt->min_chain_score = 25;
		opt->min_ksw_len = 20;
		opt->mb_size = 50000000;
	} else if (strcmp(preset, "lr") == 0) {
		opt->flag |= MB_F_LONG;
		opt->flag &= ~MB_F_PE;
		opt->min_len = 23;
		opt->bw = 500, opt->bw_long = 20000;
		opt->max_gap = 5000;
		opt->pri_ratio = 0.8f;
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
