#ifndef MB_BAM_H
#define MB_BAM_H

#include <zlib.h>
#include "bseq.h"

#ifdef __cplusplus
extern "C" {
#endif

// reading only: minibwa outputs SAM text, so there is no BAM writer here
struct mb_bam_file_s;
typedef struct mb_bam_file_s mb_bam_file_t;

mb_bam_file_t *mb_bam_open(gzFile fp); // expects that the  "BAM\1" magic is already consumed
void mb_bam_close(mb_bam_file_t *b);   // doesn't close the underlying gzFile
const char *mb_bam_hdr_text(const mb_bam_file_t *b);
int mb_bam_is_pe(const mb_bam_file_t *b); // 1 if the first record is flagged paired
int mb_bam_eof(const mb_bam_file_t *b);

// read one record, skipping secondary and supplementary ones and undoing the
// reverse complement on strand 0x10; return the length, -1 at EOF, -2 on error
int mb_bam_read1(mb_bam_file_t *b, mb_bseq1_t *s, int with_qual, int with_aux);

#ifdef __cplusplus
}
#endif

#endif
