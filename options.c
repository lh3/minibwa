#include <string.h>
#include "minibwa.h"

void mb_mopt_init(mb_mopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->min_k = 19;
	opt->n_thread = 4;
	opt->max_occ = 500;
	//opt->max_for_occ = 20;
	opt->mb_size = 500000000;
}
