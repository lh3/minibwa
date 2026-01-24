#ifndef MINIBWA_H
#define MINIBWA_H

#include <stdint.h>

#define MB_VERSION "0.0-r169-dirty"

#define MB_F_SAM              (0x1LL)    // output in the SAM format
#define MB_F_WRITE_UNMAP      (0x2LL)    // output unmapped query sequences
#define MB_F_COPY_COMMENT     (0x4LL)    // copy FASTX comments to output (SAM only)
#define MB_F_PE               (0x8LL)    // paired-end mode
#define MB_F_LONG             (0x10LL)   // long-sequence mode
#define MB_F_EQX              (0x20LL)   // = in CIGAR
#define MB_F_NO_KALLOC        (0x40LL)   // disable kalloc
#define MB_F_NO_ALN           (0x80LL)   // skip base alignment

#define MB_CIGAR_MATCH      0
#define MB_CIGAR_INS        1
#define MB_CIGAR_DEL        2
#define MB_CIGAR_N_SKIP     3
#define MB_CIGAR_SOFTCLIP   4
#define MB_CIGAR_HARDCLIP   5
#define MB_CIGAR_PADDING    6
#define MB_CIGAR_EQ_MATCH   7
#define MB_CIGAR_X_MISMATCH 8

#define MB_CIGAR_STR  "MIDNSHP=XB"

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
	int32_t max_chain_skip;
	int32_t max_chain_iter;
	int32_t min_chain_score; // min chaining score
	float chain_gap_scale;
	// hit processing options
	float mask_level;
	int32_t mask_len;
	float pri_ratio;
	int32_t best_n;
	// alignment options
	int32_t a, b;     // match, mismatch
	int32_t b_ts;     // transition mismatch
	int32_t b_ambi;   // ambiguous mismatch
	int32_t q, q2;    // gap open, long gap open
	int32_t e, e2;    // gap extension, long gap extension
	int32_t end_bonus;
	int32_t min_dp_max;
	int32_t zdrop;
	int32_t zdrop_inv;
	int32_t min_ksw_len;
	// pairing options
	int32_t max_pe_ins;
	// input/output options
	int32_t sb_len;   // number of bases for batch smem
	int32_t sb_seq;   // number of sequences for batch smem
	int32_t n_thread; // number of worker threads, excluding I/O threads
	int32_t out_n;    // max number of secondary alignments to output
	int32_t seed;
	int64_t mb_size;  // mini-batch size
	int64_t max_sw_mat;
} mb_opt_t;

struct mb_idx_s;
typedef struct mb_idx_s mb_idx_t;

typedef struct {
	uint32_t cap;                      // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2; // DP score; score of the max-scoring segment; score of the best alternate mappings
	int32_t dp_max0;                   // DP score before mb_update_dp_max() adjustment
	uint32_t n_ambi;                   // number of ambiguous bases;
	int32_t n_cigar;                   // number of cigar operations in cigar[]
	uint32_t cigar[];
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
	int32_t mapq;
	uint32_t hash;
	uint32_t rev:1, sam_pri:1, flt:1, inv:1, split:2, split_inv:1, dummy:25;
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

void mb_opt_init(mb_opt_t *opt);
int mb_opt_preset(mb_opt_t *opt, const char *preset);

mb_tbuf_t *mb_tbuf_init(int no_kalloc);
void mb_tbuf_destroy(mb_tbuf_t *b);

mb_hit_t *mb_map(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const char *seq0, int32_t *n_hit_, mb_tbuf_t *b, const char *qname);

#ifdef __cplusplus
}
#endif

#endif
