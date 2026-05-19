#include <stdlib.h>
#include "ksw2.h"

static int ksw_x86_simd_flag = -1;

#define SIMD_SSE      0x1
#define SIMD_SSE2     0x2
#define SIMD_SSE3     0x4
#define SIMD_SSSE3    0x8
#define SIMD_SSE4_1   0x10
#define SIMD_SSE4_2   0x20
#define SIMD_AVX      0x40
#define SIMD_AVX2     0x80
#define SIMD_AVX512F  0x100
#define SIMD_AVX512BW 0x200

#if !defined(_MSC_VER) && !defined(__aarch64__)
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id)
{
#if defined(__x86_64__)
	__asm__ volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	__asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif
#include <stdio.h>
int ksw_x86_simd(void)
{
#ifdef __aarch64__
	return -1;
#else
	int flag = 0, cpuid[4], max_id;
	__cpuidex(cpuid, 0, 0);
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE4_1;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE4_2;
	if ((cpuid[2]>>27&1) && (cpuid[2]>>28&1)) { // OSXSAVE && AVX CPU bits
		uint32_t xcr0;
#if defined(_MSC_VER)
		xcr0 = (uint32_t)_xgetbv(0);
#else
		__asm__ volatile("xgetbv" : "=a"(xcr0) : "c"(0) : "edx");
#endif
		if ((xcr0 & 0x6) == 0x6) { // XMM and YMM state enabled by OS
			flag |= SIMD_AVX;
			if (max_id >= 7) {
				__cpuidex(cpuid, 7, 0);
				if (cpuid[1]>>5 &1) flag |= SIMD_AVX2;
				if ((xcr0 & 0xE6) == 0xE6) { // AVX-512 opmask and ZMM state enabled
					if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
					if (cpuid[1]>>30&1) flag |= SIMD_AVX512BW;
				}
			}
		}
	}
	return flag;
#endif
}

void ksw_set_simd(void)
{
	ksw_x86_simd_flag = ksw_x86_simd();
}

void ksw_set_avx2(void)
{
	ksw_x86_simd_flag = ksw_x86_simd();
	if (ksw_x86_simd_flag >= 0)
		ksw_x86_simd_flag |= SIMD_SSE|SIMD_SSE2|SIMD_SSE3|SIMD_SSSE3|SIMD_SSE4_1|SIMD_SSE4_2|SIMD_AVX|SIMD_AVX2;
}

void ksw_extd2_simd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
					int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
	extern void ksw_extd2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
	if (ksw_x86_simd_flag >= 0 && (ksw_x86_simd_flag & SIMD_AVX2))
		ksw_extd2_avx2(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
	else
		ksw_extd2_sse(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
}
