#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minibwa.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static void process_batch(const mb_opt_t *opt, const mb_idx_t *idx, mb_tbuf_t *tbuf, int32_t n_seq, int32_t *qlen, char **seq, char **name)
{
	mb_hit_t **hit;
	int32_t *n_hit, k, i, j;
	n_hit = (int32_t*)calloc(n_seq, sizeof(int32_t));
	hit = mb_map_batch(opt, idx, n_seq, qlen, (const char**)seq, 0, n_hit, tbuf, (const char**)name);
	for (k = 0; k < n_seq; ++k) {
		for (j = 0; j < n_hit[k]; ++j) {
			mb_hit_t *h = &hit[k][j];
			const char *ctg_name = mb_idx_ctg_name(idx, h->tid);
			int64_t ctg_len = mb_idx_ctg_len(idx, h->tid);
			printf("%s\t%d\t%d\t%d\t%c\t", name[k], qlen[k], h->qs, h->qe, "+-"[h->rev]);
			printf("%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\tcg:Z:", ctg_name, (long)ctg_len, (long)h->ts, (long)h->te, h->mlen, h->blen, h->mapq);
			for (i = 0; i < h->p->n_cigar; ++i)
				printf("%d%c", h->p->cigar[i]>>4, MB_CIGAR_STR[h->p->cigar[i]&0xf]);
			putchar('\n');
			free(h->p);
		}
		free(hit[k]); free(seq[k]); free(name[k]);
	}
	free(hit);
	free(n_hit);
}

int main(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	gzFile f;
	kseq_t *ks;
	int32_t n_seq, *qlen;
	char **seq, **name;

	mb_opt_init(&opt);
	if (argc < 3) {
		fprintf(stderr, "Usage: mbmap-batch <idxPrefix> <query.fa>\n");
		return 1;
	}

	f = gzopen(argv[2], "r");
	assert(f);
	ks = kseq_init(f);
	idx = mb_idx_load(argv[1]);
	assert(idx);

	qlen = (int32_t*)malloc(opt.sb_seq * sizeof(int32_t));
	seq = (char**)malloc(opt.sb_seq * sizeof(char*));
	name = (char**)malloc(opt.sb_seq * sizeof(char*));

	n_seq = 0;
	while (kseq_read(ks) >= 0) {
		qlen[n_seq] = ks->seq.l;
		seq[n_seq] = strdup(ks->seq.s);
		name[n_seq] = strdup(ks->name.s);
		if (++n_seq >= opt.sb_seq) {
			process_batch(&opt, idx, 0, n_seq, qlen, seq, name);
			n_seq = 0;
		}
	}
	if (n_seq > 0)
		process_batch(&opt, idx, 0, n_seq, qlen, seq, name);

	free(qlen); free(seq); free(name);
	mb_idx_destroy(idx);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
