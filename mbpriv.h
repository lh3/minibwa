#ifndef MBPRIV_H
#define MBPRIV_H

#include "minibwa.h"
#include "l2bit.h"
#include "bwt.h"
#include "kommon.h"

#define MB_DBG_ALN_SEQ     (0x1LL)
#define MB_DBG_ANCHOR      (0x2LL)
#define MB_DBG_SEED        (0x4LL)
#define MB_DBG_QNAME       (0x8LL)

#define MB_SEED_LONG_JOIN  0x1
#define MB_SEED_IGNORE     0x2
#define MB_SEED_SELF       0x4 // only when a sequence is aligned against itself

struct mb_idx_s {
	l2b_t *l2b;
	mb_bwt_t *bwt;
};

typedef struct {
	int32_t sid; // tid<<1|rev
	int32_t len; // length of the anchor
	int32_t qpos; // the query coordinate of the last base in the anchor; the start base is qpos+1-len
	uint32_t flag:31, flt:1;
	int64_t tpos; // target/contig coordinate
} mb_anchor_t;

typedef struct { int64_t n, m; mb_anchor_t *a; } mb_anchor_v;

typedef struct { uint64_t x, y; } mb128_t;
void radix_sort_mb128x(mb128_t *beg, mb128_t *end);

#ifdef __cplusplus
extern "C" {
#endif

// defined in bwtgen.c
void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

// defined in seed.c
void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v);
void mb_seed_intv_batch(void *km, const mb_bwt_t *bwt, int32_t n_seq, int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v);
void mb_anchor(void *km, const mb_idx_t *idx, mb_sai_v *u, int32_t qlen, int32_t max_occ, mb_anchor_v *v);

// defined in lchain.c
mb_anchor_t *mb_lchain_dp(void *km, int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_sc, float chn_pen_gap, float chn_pen_skip,
						  int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u);
mb_anchor_t *mb_lchain_rmq(void *km, int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_sc, float chn_pen_gap, float chn_pen_skip,
						   int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u);

// defined in map-algo.c
void *mb_tbuf_km(mb_tbuf_t *b);
int32_t mb_cal_mblen(int32_t n, const mb_anchor_t *a, int32_t *blen_);
mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const l2b_t *l2b, int n_u, uint64_t *u, mb_anchor_t *a);
void mb_set_parent(void *km, float mask_level, int mask_len, int n, mb_hit_t *r, int sub_diff, int hard_mask_level);
void mb_hit_sort(void *km, int *n_regs, mb_hit_t *r);
void mb_sync_hits(void *km, int n_regs, mb_hit_t *regs);
void mb_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mb_hit_t *r);
void mb_filter_hits(const mb_opt_t *opt, int qlen, int *n_regs, mb_hit_t *regs);
int mb_squeeze_a(void *km, int n_regs, mb_hit_t *regs, mb_anchor_t *a);
void mb_split_hit(mb_hit_t *r, mb_hit_t *r2, int n, int qlen, mb_anchor_t *a, const l2b_t *l2b);

mb_hit_t *mb_map_sai(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const uint8_t *seq, mb_sai_v *u, int32_t *n_hit_, mb_tbuf_t *b, const char *qname);

// defined in format.c
void mb_fmt_paf_basic(kstring_t *s, const l2b_t *l2b, int64_t qlen, const mb_hit_t *p, const char *qname);

// defined in align.c
mb_hit_t *mb_align_skeleton(void *km, const mb_opt_t *opt, const mb_idx_t *mi, int qlen, const uint8_t *seq, int *n_regs_, mb_hit_t *regs, mb_anchor_t *a);

// Fast log2 approximation (from minimap2)
static inline float mb_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

static inline uint64_t mb_hash64(uint64_t x)
{
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

static inline uint32_t mb_hash_str(const char *s)
{
	uint32_t h = 2166136261U;
	const unsigned char *t = (const unsigned char*)s;
	for (; *t; ++t)
		h ^= *t, h *= 16777619;
	return h;
}

#ifdef __cplusplus
}
#endif

#endif
