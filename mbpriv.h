#ifndef MBPRIV_H
#define MBPRIV_H

#include "minibwa.h"
#include "l2bit.h"
#include "bwt.h"
#include "kommon.h"
#include "bseq.h"

#define MB_DBG_ALN_SEQ     (0x1LL)
#define MB_DBG_ANCHOR      (0x2LL)
#define MB_DBG_SEED        (0x4LL)
#define MB_DBG_QNAME       (0x8LL)
#define MB_DBG_ALN_PE      (0x10LL)
#define MB_DBG_AN_POS      (0x20LL)

#define MB_SEED_LONG_JOIN  0x1
#define MB_SEED_IGNORE     0x2

struct mb_idx_s {
	int32_t is_meth;
	l2b_t *l2b;
	mb_bwt_t *bwt;
};

typedef struct {
	int32_t sid; // tid<<1|rev
	int32_t len; // length of the anchor
	int32_t qpos; // the query coordinate of the last base in the anchor; the start base is qpos+1-len; flipped by rev
	uint32_t flag:31, flt:1;
	int64_t tpos; // target/contig coordinate
} mb_anchor_t;

typedef struct { int64_t n, m; mb_anchor_t *a; } mb_anchor_v;

typedef struct {
	int lo, hi;      // lower and upper bounds within which a read pair is considered to be properly paired
	int failed;      // non-zero if the orientation is not supported by sufficient data
	double avg, std; // mean and stddev of the insert size distribution
} mb_pestat_t;

typedef struct { uint64_t x, y; } mb128_t;
void radix_sort_mb128x(mb128_t *beg, mb128_t *end);

#ifdef __cplusplus
extern "C" {
#endif

// defined in options.c
void mb_opt_adap(const mb_opt_t *opt0, int32_t len, mb_opt_t *opt);

// defined in bwtgen.c
void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

// defined in seed.c
void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v);
void mb_seed_intv_batch(void *km, const mb_bwt_t *bwt, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v);
double mb_anchor(void *km, const mb_idx_t *idx, mb_sai_v *u, int32_t min_len, int32_t qlen, const uint8_t *qseq, l2b_meth_t mt, int32_t max_occ, mb_anchor_v *v);
void mb_anchor_sort(const l2b_t *l2b, int64_t n_a, mb_anchor_t *a);

// defined in lchain.c
mb_anchor_t *mb_lchain_dp(void *km, const l2b_t *l2b, int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_sc, float chn_pen_gap,
						  int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u);

// defined in map-algo.c
void *mb_tbuf_km(mb_tbuf_t *b);
int32_t mb_cal_mblen(int32_t n, const mb_anchor_t *a, int32_t *blen_);
mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const l2b_t *l2b, int n_u, uint64_t *u, mb_anchor_t *a);
void mb_sync_high_cov(int32_t n, mb_hit_t *h);
void mb_set_parent(void *km, float mask_level, int mask_len, int n, mb_hit_t *r, int sub_diff, int hard_mask_level);
void mb_set_sam_pri(int32_t n, mb_hit_t *r, int32_t is_primary5);
void mb_hit_sort(void *km, int *n_regs, mb_hit_t *r);
void mb_sync_hits(void *km, int n_regs, mb_hit_t *regs);
void mb_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mb_hit_t *r);
void mb_filter_hits(const mb_opt_t *opt, int qlen, int *n_regs, mb_hit_t *regs);
int mb_squeeze_a(void *km, int n_regs, mb_hit_t *regs, mb_anchor_t *a);
void mb_split_hit(mb_hit_t *r, mb_hit_t *r2, int n, int qlen, mb_anchor_t *a, const l2b_t *l2b);
void mb_set_mapq(void *km, int32_t qlen, int n_regs, mb_hit_t *regs, int min_chain_sc, int match_sc, int is_sr, int max_sr_len);

mb_hit_t *mb_map_sai(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const char *seq, l2b_meth_t mt, mb_sai_v *u, int32_t *n_hit_, mb_tbuf_t *b, const char *qname);

void radix_sort_mb64(uint64_t *st, uint64_t *en);
void radix_sort_mb128x(mb128_t *st, mb128_t *en);

// in cs.c
void mb_write_cs_ds(void *km, kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mb_hit_t *r, int is_ds);
void mb_write_MD(void *km, kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mb_hit_t *r);

// defined in format.c
void mb_fmt_paf(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx);
int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[]);
void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen);
char *mb_escape(char *s);

// defined in align.c
mb_hit_t *mb_align_skeleton(void *km, const mb_opt_t *opt, const mb_idx_t *mi, int qlen, const uint8_t *seq, l2b_meth_t mt, int *n_regs_, mb_hit_t *regs, mb_anchor_t *a);
void mb_append_cigar(mb_hit_t *r, uint32_t n_cigar, const uint32_t *cigar);
void mb_update_extra(void *km, mb_hit_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e, uint64_t opt_flag, int log_gap);

// defined in pe.c
void mb_pestat(void *km, const mb_opt_t *opt, int32_t n_seg, const int32_t *seg_off, const int32_t *seg_cnt, const int32_t *n_hit, mb_hit_t *const *hit, mb_pestat_t pes[4]);
void mb_pair(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], int32_t qlen[2], char *const qseq[2]);

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

static inline void mb_seq_rev(uint32_t len, uint8_t *seq)
{
	kom_reverse(uint8_t, len, seq);
}

static inline int32_t mb_is_sr_mode(const mb_opt_t *opt, int32_t qlen)
{
	return (opt->flag & MB_F_LONG) || ((opt->flag & MB_F_ADAP) && qlen > opt->max_sr_len)? 0 : 1;
}

#ifdef __cplusplus
}
#endif

#endif
