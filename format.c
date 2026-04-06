#include "mbpriv.h"
#include "kommon.h"

void mb_fmt_paf(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx)
{
	kom_sprintf_lite(s, "%s", t->name);
	if (n_seg > 1 && seg_idx >= 0)
		kom_sprintf_lite(s, "/%d", seg_idx + 1);
	kom_sprintf_lite(s, "\t%ld", (long)t->l_seq);
	if (p == 0) { // for unmapped reads
		kom_sprintf_lite(s, "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n");
		return;
	}
	kom_sprintf_lite(s, "\t%d\t%d\t%c\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\ttp:A:%c\ts1:i:%d\tcm:i:%d",
		p->qs, p->qe, p->rev? '-' : '+', l2b->ctg[p->tid].name, (long)l2b->ctg[p->tid].len, (long)p->ts, (long)p->te,
		p->mlen, p->blen, p->mapq, p->parent == p->id? 'P' : 'S', p->score, p->cnt);
	if (p->parent == p->id) kom_sprintf_lite(s, "\ts2:i:%d", p->subsc);
	if (p->p) {
		int32_t nm = p->blen - p->mlen + p->p->n_ambi;
		kom_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tms:i:%d\tm2:i:%d", nm, p->p->dp_score, p->p->dp_max, p->p->dp_max2);
		if (p->p->n_cigar > 0) {
			int32_t i;
			kom_sprintf_lite(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i)
				kom_sprintf_lite(s, "%d%c", p->p->cigar[i]>>4, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
		}
	}
	if (opt_flag & (MB_F_WRITE_DS|MB_F_WRITE_CS|MB_F_WRITE_MD))
		kom_sprintf_lite(s, "\t%s", (char*)&p->p->cigar[p->p->n_cigar]);
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
}
