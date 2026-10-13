/* ksw2_extd2_dispatch.c — runtime ISA dispatch for ksw_extd2_sse.
 *
 * minibwa calls ksw_extd2_sse() throughout align.c. This file makes that symbol a
 * thin dispatcher that, on x86, selects the widest available implementation at
 * RUNTIME (so a single binary runs everywhere and only uses AVX-512/AVX2 where
 * the CPU supports it):
 *
 *   AVX-512BW  -> ksw_extd2_avx512   (ksw2_extd2_wide.c built -mavx512bw)
 *   AVX2       -> ksw_extd2_avx2     (ksw2_extd2_wide.c built -mavx2)
 *   otherwise  -> extd2_ref_impl     (original ksw2_extd2_sse.c, SSE4.2)
 *
 * The wide kernels only widen minibwa's gap-fill path (cigar + approx + left);
 * every other flag combination delegates to extd2_ref_impl, so behaviour is
 * unchanged. On non-x86 (e.g. arm64/NEON) the dispatcher is just extd2_ref_impl.
 *
 * Override for testing the fallback path on a capable CPU:
 *   MINIBWA_EXTD2_ISA=sse|avx2|avx512
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ksw2.h"

#define KSW_EXTD2_ARGS \
    void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, \
    int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2, \
    int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez
#define KSW_EXTD2_CALL \
    km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez

void extd2_ref_impl(KSW_EXTD2_ARGS);

#if defined(__x86_64__) || defined(__i386__)
void ksw_extd2_avx2(KSW_EXTD2_ARGS);
void ksw_extd2_avx512(KSW_EXTD2_ARGS);

typedef void (*ksw_extd2_fn)(KSW_EXTD2_ARGS);
static ksw_extd2_fn g_extd2 = 0;

static ksw_extd2_fn ksw_extd2_resolve(void) {
    const char *force = getenv("MINIBWA_EXTD2_ISA");
    __builtin_cpu_init();
    if (force) {
        if (!strcmp(force, "sse"))    return extd2_ref_impl;
        if (!strcmp(force, "avx2")) {
            if (__builtin_cpu_supports("avx2")) return ksw_extd2_avx2;
            fprintf(stderr, "[W::ksw_extd2] MINIBWA_EXTD2_ISA=avx2 but CPU lacks AVX2; using SSE\n");
            return extd2_ref_impl;
        }
        if (!strcmp(force, "avx512")) {
            if (__builtin_cpu_supports("avx512bw")) return ksw_extd2_avx512;
            fprintf(stderr, "[W::ksw_extd2] MINIBWA_EXTD2_ISA=avx512 but CPU lacks AVX-512BW; falling back\n");
            if (__builtin_cpu_supports("avx2")) return ksw_extd2_avx2;
            return extd2_ref_impl;
        }
        fprintf(stderr, "[W::ksw_extd2] ignoring unrecognized MINIBWA_EXTD2_ISA='%s' "
                        "(expected sse|avx2|avx512)\n", force);
    }
    if (__builtin_cpu_supports("avx512bw")) return ksw_extd2_avx512;
    if (__builtin_cpu_supports("avx2"))     return ksw_extd2_avx2;
    return extd2_ref_impl;
}

/* resolve once before main() (no per-call branch, no thread race) */
__attribute__((constructor)) static void ksw_extd2_init(void) { g_extd2 = ksw_extd2_resolve(); }

void ksw_extd2_sse(KSW_EXTD2_ARGS) {
    if (__builtin_expect(g_extd2 == 0, 0)) g_extd2 = ksw_extd2_resolve(); /* belt-and-suspenders */
    g_extd2(KSW_EXTD2_CALL);
}
#else
void ksw_extd2_sse(KSW_EXTD2_ARGS) { extd2_ref_impl(KSW_EXTD2_CALL); }
#endif
