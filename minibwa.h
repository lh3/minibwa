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
	int32_t min_k;    // min seed/k-mer length
	int32_t n_thread; // number of worker threads, excluding I/O threads
	int32_t max_occ;  // max SA interval size during seeding
	int64_t mb_size;  // mini-batch size
} mb_mopt_t;

struct mb_idx_s;
typedef struct mb_idx_s mb_idx_t;

typedef struct {
	uint32_t cap;
	int32_t n_cigar;
	uint64_t cigar[];
} mb_extra_t;

typedef struct {
	int64_t tid;
	int64_t ts, te;
	int32_t qs, qe;
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

mb_tbuf_t *mb_tbuf_init(void);
void mb_tbuf_destroy(mb_tbuf_t *b);

#ifdef __cplusplus
}
#endif

#endif
