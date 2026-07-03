#ifndef HL_HEADER
#define HL_HEADER

#include <stdint.h>

#define kom_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define kom_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define kom_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

// make enough room to write ptr[i]
#define kom_grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = kom_realloc(type, (ptr), (__m)); \
		} \
	} while (0)

#define kom_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define kom_roundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

#define kom_reverse(type, n, a) do { \
		size_t i; \
		for (i = 0; i < (n)>>1; ++i) { \
			type t = a[(n) - 1 - i]; \
			a[(n) - 1 - i] = a[i], a[i] = t; \
		} \
	} while (0)

#define kom_assert(cond, msg) if ((cond) == 0) kom_panic(__func__, (msg))

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct { uint64_t x, y; } kom128_t;

#ifdef __cplusplus
extern "C" {
#endif

extern int kom_verbose, kom_dbg, kom_dbg_flag;

char *kom_strdup(const char *src);
int64_t kom_parse_num(const char *str, char **q);
void kom_panic(const char *func, const char *msg);

int64_t kom_sprintf_lite(kstring_t *s, const char *fmt, ...);
int64_t km_sprintf_lite(void *km, kstring_t *s, const char *fmt, ...);

double kom_cputime(void);
double kom_realtime(void); // call at the beginning to reset the timer
long kom_peakrss(void);
double kom_percent_cpu(void);

extern uint8_t kom_nt4_table[256], kom_comp_table[256];
void kom_revcomp(uint64_t len, char *seq);

static inline uint64_t kom_splitmix64(uint64_t *x)
{
	uint64_t z = ((*x) += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

static inline double kom_u64todbl(uint64_t x)
{
	union { uint64_t i; double d; } u;
	u.i = 0x3FFULL << 52 | x >> 12;
	return u.d - 1.0;
}

#ifdef __cplusplus
}
#endif

#endif // ~HL_HEADER
