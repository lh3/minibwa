#ifndef S2N_LITE_H
#define S2N_LITE_H

#include <arm_neon.h>

typedef uint8x16_t __m128i;

static inline __m128i _mm_load_si128(const __m128i *ptr) { return vld1q_u8((const uint8_t*)ptr); }
static inline __m128i _mm_loadu_si128(const __m128i *ptr) { return vld1q_u8((const uint8_t*)ptr); }
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { vst1q_u8((uint8_t*)ptr, a); }
static inline void _mm_storeu_si128(__m128i *ptr, __m128i a) { vst1q_u8((uint8_t*)ptr, a); }
static inline __m128i _mm_setzero_si128(void) { return vdupq_n_u8(0); }
static inline __m128i _mm_or_si128(__m128i a, __m128i b) { return vorrq_u8(a, b); }
static inline __m128i _mm_and_si128(__m128i a, __m128i b) { return vandq_u8(a, b); }
static inline __m128i _mm_andnot_si128(__m128i a, __m128i b) { return vbicq_u8(b, a); }

#define _mm_slli_si128(a, imm8) vextq_u8(_mm_setzero_si128(), (a), 16 - (imm8))
#define _mm_srli_si128(a, imm8) vextq_u8((a), _mm_setzero_si128(), (imm8))

static inline __m128i _mm_blendv_epi8(__m128i a, __m128i b, __m128i mask) { return vbslq_u8(vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(mask), 7)), b, a); }

static inline __m128i _mm_set1_epi8(int a) { return vdupq_n_u8(a); }
static inline __m128i _mm_add_epi8(__m128i a, __m128i b) { return vaddq_u8(a, b); }
static inline __m128i _mm_adds_epu8(__m128i a, __m128i b) { return vqaddq_u8(a, b); }
static inline __m128i _mm_sub_epi8(__m128i a, __m128i b) { return vsubq_u8(a, b); }
static inline __m128i _mm_subs_epu8(__m128i a, __m128i b) { return vqsubq_u8(a, b); }
static inline __m128i _mm_cmpeq_epi8(__m128i a, __m128i b) { return vceqq_u8(a, b); }
static inline __m128i _mm_cmpgt_epi8(__m128i a, __m128i b) { return vcgtq_s8(vreinterpretq_s8_u8(a), vreinterpretq_s8_u8(b)); }
static inline __m128i _mm_max_epi8(__m128i a, __m128i b) { return vreinterpretq_u8_s8(vmaxq_s8(vreinterpretq_s8_u8(a), vreinterpretq_s8_u8(b))); }
static inline __m128i _mm_min_epi8(__m128i a, __m128i b) { return vreinterpretq_u8_s8(vminq_s8(vreinterpretq_s8_u8(a), vreinterpretq_s8_u8(b))); }
static inline __m128i _mm_max_epu8(__m128i a, __m128i b) { return vmaxq_u8(a, b); }
static inline __m128i _mm_min_epu8(__m128i a, __m128i b) { return vminq_u8(a, b); }

static inline __m128i _mm_set1_epi16(int a) { return vreinterpretq_u8_s16(vdupq_n_s16(a)); }
static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_u16(vcgtq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_max_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vmaxq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_adds_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vqaddq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_subs_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vqsubq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_subs_epu16(__m128i a, __m128i b) { return vreinterpretq_u8_u16(vqsubq_u16(vreinterpretq_u16_u8(a), vreinterpretq_u16_u8(b))); }

#define _mm_extract_epi16(a, imm8) vgetq_lane_s16(vreinterpretq_s16_u8(a), (imm8))
#define _mm_insert_epi16(a, b, imm8) vreinterpretq_u8_s16(vsetq_lane_s16((b), vreinterpretq_s16_u8(a), (imm8)))

static inline __m128i _mm_set1_epi32(int a) { return vreinterpretq_u8_s32(vdupq_n_s32(a)); }
static inline __m128i _mm_cvtsi32_si128(int a) { return vreinterpretq_u8_s32(vsetq_lane_s32(a, vdupq_n_s32(0), 0)); }
static inline __m128i _mm_setr_epi32(int a, int b, int c, int d) {
	int32_t x[4] = {a, b, c, d};
	return vld1q_u8((const uint8_t*)x);
}
static inline __m128i _mm_cmpgt_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_u32(vcgtq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_max_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vmaxq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_add_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vaddq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_sub_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vsubq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }

#define _mm_insert_epi32(a, b, imm8) vreinterpretq_u8_s32(vsetq_lane_s32((b), vreinterpretq_s32_u8(a), (imm8)))

#define _mm_prefetch(p, i) __builtin_prefetch(p, 0, (i))

#endif
