#include <string.h>
#include "minibwa.h"

void mb_mopt_init(mb_mopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	// seeding options
	opt->min_len = 19;
	opt->max_sub_occ = 10;
	opt->max_occ = 500;
	// general algorithm options
	opt->bw = 500, opt->bw_long = 20000;
	opt->max_gap = 5000;
	// chaining options
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->min_chain_score = 40;
	opt->rmq_inner_dist = 1000;
	opt->rmq_size_cap = 100000;
	opt->chn_pen_skip = 0.05f;
	// hit processing options
	opt->mask_level = 0.5f;
	opt->mask_len = INT32_MAX;
	opt->sub_diff = 0;
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;
	// alignment options
	opt->a = 2,  opt->b = 8;
	opt->q = 12, opt->q2 = 26;
	opt->e = 2,  opt->e2 = 1;
	opt->b_ambi = 1;
	opt->end_bonus = -1;
	opt->min_dp_max = 80;
	opt->zdrop = 400;
	opt->zdrop_inv = 200;
	opt->min_ksw_len = 200;
	// I/O options
	opt->n_thread = 1; // TODO: use 4 by default in future
	opt->seed = 11;
	opt->mb_size = 500000000;
	opt->max_sw_mat = 100000000;
}

int mb_set_preset(const char *preset, mb_mopt_t *opt)
{
	if (preset == 0) {
		mb_mopt_init(opt);
	} else if (strcmp(preset, "lr:hq") == 0) { // to be added
	} else if (strcmp(preset, "asm5") == 0) { // to be added
	} else if (strcmp(preset, "sr") == 0) {
		mb_mopt_init(opt);
		opt->flag |= MB_F_FRAG_MODE;
		opt->max_gap = 100;
		opt->end_bonus = 10;
		opt->bw = opt->bw_long = 100;
		opt->min_chain_score = 25;
		opt->mb_size = 50000000;
	} else {
		return -1;
	}
	return 0;
}
