/* ksw2_wide_simd.h — shared SIMD width/ISA abstraction for the wide ksw2 kernels.
 *
 * Two back-ends, selected at compile time:
 *   AVX-512BW : W = 64 bytes  (no VBMI required)
 *   AVX2      : W = 32 bytes
 *
 * Provides byte-wise ops plus the one width-specific primitive the anti-diagonal
 * recurrence needs: a FULL-REGISTER "shift left by one byte" (the SSE
 * _mm_slli_si128 is per-128-bit-lane, so on AVX2/AVX-512 it must cross lanes).
 */
#ifndef KSW2_WIDE_SIMD_H
#define KSW2_WIDE_SIMD_H

#include <stdint.h>

#if defined(__AVX512BW__)
#  include <immintrin.h>
#  define KSW_VBYTES 64
typedef __m512i ksw_vec;
static inline __m512i ksw_vslli1(__m512i a) {        /* result[p]=a[p-1], [0]=0 */
    __m512i lo = _mm512_bslli_epi128(a, 1);
    __m512i hb = _mm512_bsrli_epi128(a, 15);
    const __m512i idx = _mm512_set_epi64(5, 4, 3, 2, 1, 0, 0, 0);
    __m512i hb_sh = _mm512_maskz_permutexvar_epi64(0xFC, idx, hb);
    return _mm512_or_si512(lo, hb_sh);
}
static inline __m512i ksw_setlane0(int b) {
    return _mm512_zextsi128_si512(_mm_cvtsi32_si128((unsigned char)b));
}
#  define KSW_LDU(p)        _mm512_loadu_si512((const void *)(p))
#  define KSW_STU(p, x)     _mm512_storeu_si512((void *)(p), (x))
#  define KSW_SET1(x)       _mm512_set1_epi8(x)
#  define KSW_ZERO()        _mm512_setzero_si512()
#  define KSW_ADD8(a, b)    _mm512_add_epi8((a), (b))
#  define KSW_SUB8(a, b)    _mm512_sub_epi8((a), (b))
#  define KSW_MAXS8(a, b)   _mm512_max_epi8((a), (b))
#  define KSW_MAXU8(a, b)   _mm512_max_epu8((a), (b))
#  define KSW_MINU8(a, b)   _mm512_min_epu8((a), (b))
#  define KSW_MINS8(a, b)   _mm512_min_epi8((a), (b))
#  define KSW_OR(a, b)      _mm512_or_si512((a), (b))

#elif defined(__AVX2__)
#  include <immintrin.h>
#  define KSW_VBYTES 32
typedef __m256i ksw_vec;
static inline __m256i ksw_vslli1(__m256i a) {
    __m256i perm = _mm256_permute2x128_si256(a, a, 0x08);
    return _mm256_alignr_epi8(a, perm, 15);
}
static inline __m256i ksw_setlane0(int b) {
    return _mm256_zextsi128_si256(_mm_cvtsi32_si128((unsigned char)b));
}
#  define KSW_LDU(p)        _mm256_loadu_si256((const __m256i *)(p))
#  define KSW_STU(p, x)     _mm256_storeu_si256((__m256i *)(p), (x))
#  define KSW_SET1(x)       _mm256_set1_epi8(x)
#  define KSW_ZERO()        _mm256_setzero_si256()
#  define KSW_ADD8(a, b)    _mm256_add_epi8((a), (b))
#  define KSW_SUB8(a, b)    _mm256_sub_epi8((a), (b))
#  define KSW_MAXS8(a, b)   _mm256_max_epi8((a), (b))
#  define KSW_MAXU8(a, b)   _mm256_max_epu8((a), (b))
#  define KSW_MINU8(a, b)   _mm256_min_epu8((a), (b))
#  define KSW_MINS8(a, b)   _mm256_min_epi8((a), (b))
#  define KSW_OR(a, b)      _mm256_or_si256((a), (b))

#else
#  error "wide ksw2 kernels require -mavx2 or -mavx512bw"
#endif

#endif /* KSW2_WIDE_SIMD_H */
