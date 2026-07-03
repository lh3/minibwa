#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minibwa.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_opt_init(&opt);

	if (argc < 3) {
		fprintf(stderr, "Usage: mbmap-sgl <idxPrefix> <query.fa>\n");
		return 1;
	}

	// open query file for reading; you may use your favorite FASTA/Q parser
	gzFile f = gzopen(argv[2], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);

	// open index reader
	mb_idx_t *idx = mb_idx_load(argv[1], 0);
	assert(idx);
	while (kseq_read(ks) >= 0) {
		mb_hit_t *hit;
		int32_t i, j, n_hit;
		hit = mb_map(&opt, idx, ks->seq.l, ks->seq.s, 0, &n_hit, 0, ks->name.s);
		for (j = 0; j < n_hit; ++j) {
			mb_hit_t *h = &hit[j];
			printf("%s\t%lu\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, h->qs, h->qe, "+-"[h->rev]);
			printf("%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\tcg:Z:", mb_idx_ctg_name(idx, h->tid), (long)mb_idx_ctg_len(idx, h->tid), (long)h->ts, (long)h->te, h->mlen, h->blen, h->mapq);
			for (i = 0; i < h->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
				printf("%d%c", h->p->cigar[i]>>4, MB_CIGAR_STR[h->p->cigar[i]&0xf]);
			putchar('\n');
			free(h->p);
		}
		free(hit);
	}
	mb_idx_destroy(idx);
	kseq_destroy(ks); // close the query file
	gzclose(f);
	return 0;
}
