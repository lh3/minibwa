#ifndef MINIBWA_H
#define MINIBWA_H

#include <stdint.h>

#define MB_VERSION "0.0"

#define MB_F_SAM              (0x1LL)    // output in the SAM format
#define MB_F_WRITE_UNMAP      (0x2LL)    // output unmapped query sequences
#define MB_F_COPY_COMMENT     (0x4LL)    // copy FASTX comments to output (SAM only)
#define MB_F_FRAG_MODE        (0x8LL)    // fragment/paired-end mode

typedef struct {
	uint64_t flag;
	// seeding options
	int32_t min_len; // min seed length
	int32_t max_sub_occ; // look for shorter seed if smem occ below this value
	int32_t max_occ; // max interval occurrence
	// general algorithm options
	int32_t bw, bw_long; // bandwidth
	int32_t max_gap; // break a chain if there are no seeds in a max_gap window
	// chaining options
	int32_t max_chain_skip, max_chain_iter;
	int32_t min_chain_score; // min chaining score
	int32_t rmq_inner_dist; // RMQ inner distance
	int32_t rmq_size_cap; // RMQ size cap
	float chn_pen_gap; // gap penalty coefficient
	float chn_pen_skip; // skip penalty coefficient
	// hit processing options
	float mask_level;
	int32_t mask_len;
	int32_t sub_diff;
	float pri_ratio;
	int32_t best_n;
	float alt_diff_frac;
	// input/output options
	int32_t n_thread; // number of worker threads, excluding I/O threads
	int64_t mb_size;  // mini-batch size
} mb_mopt_t;

struct mb_idx_s;
typedef struct mb_idx_s mb_idx_t;

typedef struct {
	uint32_t cap;
	int32_t n_cigar;
	uint64_t cigar[];
} mb_extra_t;

#define MB_PARENT_UNSET   (-1)
#define MB_PARENT_TMP_PRI (-2)

typedef struct {
	int32_t id;             // ID for internal uses
	int32_t cnt;            // number of anchors
	int32_t score, score0;  // chaining score; score0 is the original chaining score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t qs, qe;         // query start and end
	int64_t tid;            // target ID (the original tid, NOT stranded)
	int64_t ts, te;         // target start and end
	int32_t parent, n_sub, subsc;
	int32_t mlen, blen;
	uint32_t hash;
	uint32_t rev:1, sam_pri:1, inv:1, split:2, seg_split:1, strand_retained:1, dummy:25;
	int32_t seg_id;
	mb_extra_t *p;
} mb_hit_t;

struct mb_tbuf_s;
typedef struct mb_tbuf_s mb_tbuf_t;

#ifdef __cplusplus
extern "C" {
#endif

mb_idx_t *mb_idx_load(const char *prefix);
void mb_idx_destroy(mb_idx_t *idx);

void mb_mopt_init(mb_mopt_t *opt);
int mb_set_preset(const char *preset, mb_mopt_t *opt);

mb_tbuf_t *mb_tbuf_init(void);
void mb_tbuf_destroy(mb_tbuf_t *b);

#ifdef __cplusplus
}
#endif

#endif
