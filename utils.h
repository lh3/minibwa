#ifndef MB_UTILS_H
#define MB_UTILS_H

#include <stdint.h>

#define mb_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define mb_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define mb_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define mb_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define mb_roundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

#define mb_assert(cond, msg) if ((cond) == 0) mb_panic(__func__, (msg))

#ifdef __cplusplus
extern "C" {
#endif

void mb_panic(const char *fn, const char *msg);

static inline uint64_t mb_splitmix64(uint64_t *x)
{
	uint64_t z = ((*x) += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

#ifdef __cplusplus
}
#endif

#endif
