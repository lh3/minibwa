#ifndef MB_BSEQ_H
#define MB_BSEQ_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

struct mb_bseq_file_s;
typedef struct mb_bseq_file_s mb_bseq_file_t;

typedef struct {
	uint64_t l_seq, id; // FIXME: bseq doesn't fully support 64-bit integers yet
	char *name, *seq, *qual, *comment;
} mb_bseq1_t;

mb_bseq_file_t *mb_bseq_open(const char *fn);
void mb_bseq_close(mb_bseq_file_t *fp);
mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_);
mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_);
int mb_bseq_eof(mb_bseq_file_t *fp);

extern unsigned char mb_nt4_table[256];
extern unsigned char mb_comp_table[256];

static inline int mb_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int mb_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = mb_qname_len(s1);
	l2 = mb_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

static inline void mb_revcomp_bseq(mb_bseq1_t *s)
{
	int i, t, l = s->l_seq;
	for (i = 0; i < l>>1; ++i) {
		t = s->seq[l - i - 1];
		s->seq[l - i - 1] = mb_comp_table[(uint8_t)s->seq[i]];
		s->seq[i] = mb_comp_table[t];
	}
	if (l&1) s->seq[l>>1] = mb_comp_table[(uint8_t)s->seq[l>>1]];
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

#ifdef __cplusplus
}
#endif

#endif
