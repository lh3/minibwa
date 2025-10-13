#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

void mb_panic(const char *fn, const char *msg)
{
	fprintf(stderr, "[E::%s] %s ABORT!\n", fn, msg);
	abort();
}
