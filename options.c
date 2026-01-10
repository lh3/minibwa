#include <string.h>
#include "minibwa.h"

void mb_mopt_init(mb_mopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));

	// seed options
	opt->sopt.min_len = 19;
	opt->sopt.max_occ = 500;
	opt->sopt.max_sub_occ = 10;

	opt->n_thread = 4;
	opt->mb_size = 500000000;
}
