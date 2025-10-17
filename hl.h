#ifndef HL_HEADER
#define HL_HEADER

#include <stdint.h>

#define hl_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define hl_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define hl_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

// make enough room to write ptr[i]
#define hl_grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = hl_realloc(type, (ptr), (__m)); \
		} \
	} while (0)

#define hl_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define hl_roundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

#define hl_assert(cond, msg) if ((cond) == 0) hl_panic(__func__, (msg))

#ifdef __cplusplus
extern "C" {
#endif

char *hl_strdup(const char *src);
int64_t hl_parse_num(const char *str, char **q);
void hl_panic(const char *func, const char *msg);

extern uint8_t hl_nt4_table[256], hl_comp_table[256];
void hl_revcomp(uint64_t len, char *seq);

static inline uint64_t hl_splitmix64(uint64_t *x)
{
	uint64_t z = ((*x) += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

#ifdef __cplusplus
}
#endif

#endif // ~HL_HEADER
