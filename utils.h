#ifndef MB_UTILS_H
#define MB_UTILS_H

#include <stdint.h>

#define mb_assert(cond, msg) if ((cond) == 0) mb_panic(__func__, (msg))

#ifdef __cplusplus
extern "C" {
#endif

void mb_panic(const char *fn, const char *msg);
uint64_t mb_read_huge(FILE *fp, uint64_t size, void *ptr);

#ifdef __cplusplus
}
#endif

#endif
