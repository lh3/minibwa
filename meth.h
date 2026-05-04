#ifndef MB_METH_H
#define MB_METH_H

#include <stdint.h>

#include "l2bit.h"
#include "minibwa.h"
#include "bseq.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Tag prefixes (5 chars each, NUL-terminated) used by --meth output. */
#define MB_METH_YS_PREFIX "YS:Z:"
#define MB_METH_YC_PREFIX "YC:Z:"
#define MB_METH_YD_PREFIX "YD:Z:"

/* Chrom-consolidation map: f<name>/r<name> → <name>; built once from l2b. */
typedef struct mb_meth_cmap_s {
	int      n_internal;    /* == l2b->n_ctg */
	int      n_output;      /* consolidated names */
	int     *out_tid;       /* [n_internal] -> output index */
	char    *direction;     /* [n_internal] 'f', 'r', or 0 */
	char   **output_names;  /* [n_output] */
	int64_t *output_lens;   /* [n_output] */
} mb_meth_cmap_t;

mb_meth_cmap_t *mb_meth_cmap_build(const l2b_t *l2b);
void            mb_meth_cmap_free(mb_meth_cmap_t *m);

/* Bisulfite read ingest: project R1 C->T / R2 G->A in place; stash pre-projection
 * SEQ in YS:Z and conversion direction in YC:Z, appended to t->comment. */
void mb_meth_ingest(int n_seq, mb_bseq1_t *seq, int frag_mode, int split_files);

/* Chimera QC: flag matching-strand reads (--set-as-failed) and short-longest-M chimeras
 * (unless --do-not-penalize-chimeras); propagate 0x200 across both segments. */
void mb_meth_apply_qc(const mb_opt_t *opt, const mb_meth_cmap_t *cmap,
                      const mb_bseq1_t *seq, int32_t n_seg,
                      const int32_t *n_hit, mb_hit_t *const *hit);

#ifdef __cplusplus
}
#endif

#endif
