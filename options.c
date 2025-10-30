#include <string.h>
#include "minibwa.h"

void mb_mopt_init(mb_mopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->min_k = 19;
	opt->n_thread = 1;
	opt->mb_size = 500000000;
}
