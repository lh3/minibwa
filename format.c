#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mbpriv.h"
#include "kommon.h"

static char mb_rg_id[256];

/**************
 * PAF output *
 **************/

static inline void write_tags(kstring_t *s, const mb_hit_t *p)
{
	int32_t nm = p->blen - p->mlen + p->p->n_ambi;
	kom_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tms:i:%d\tmd:i:%d", nm, p->p->dp_score, p->p->dp_max0, p->p->dp_max - p->p->dp_max2);
}

/* Bismark XR/XG/XM. Under directional --meth: XR follows mate (R1=CT,
 * R2=GA, SE=CT); XG = "GA" iff (is_R2 XOR r->rev), reproducing Bismark's
 * (XR,XG) -> {OT,OB,CTOT,CTOB} encoding.
 *
 * XM is a per-base methylation-call string of length qlen, written in
 * SAM SEQ orientation. SAM SEQ is always in top-strand-ref orientation
 * (sam_write_sq revcomps when r->rev), so a CIGAR-ordered walk over SEQ
 * positions produces a string that already aligns with SEQ in display
 * order; no end-of-walk reversal is needed.
 *
 * Alphabet: . non-cytosine; z/Z CpG un/methylated; x/X CHG; h/H CHH;
 * u/U cytosine in unknown context (off the reference window or N base). */
static void write_meth_tags(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t,
							int n_seg, int seg_idx, const mb_hit_t *r)
{
	int is_r2 = (n_seg == 2 && seg_idx == 1);
	int xg_is_ga = is_r2 ^ r->rev;
	const char *xr = is_r2 ? "GA" : "CT";
	const char *xg = xg_is_ga ? "GA" : "CT";
	int qlen = t->l_seq;
	int64_t tlen = l2b->ctg[r->tid].len;
	int64_t fetch_st = r->ts >= 2 ? r->ts - 2 : 0;
	int64_t fetch_en = r->te + 2 <= tlen ? r->te + 2 : tlen;
	uint8_t *ref = (uint8_t*)kom_malloc(uint8_t, fetch_en - fetch_st);
	char *xm = (char*)kom_malloc(char, qlen + 1);
	int seq_pos = r->qs;
	int64_t ref_pos = r->ts;
	uint32_t k;
	int j;

	l2b_getseq(l2b, r->tid, fetch_st, fetch_en, ref);
	for (j = 0; j < qlen; ++j) xm[j] = '.';
	xm[qlen] = 0;

	for (k = 0; k < r->p->n_cigar; ++k) {
		int op = r->p->cigar[k] & 0xf;
		int len = r->p->cigar[k] >> 4;
		if (op == 0 || op == 7 || op == 8) { /* M/=/X */
			for (j = 0; j < len; ++j) {
				int read_b, off, c1, c2;
				char meth = '.', unmeth = '.', call = '.';
				if (ref_pos < fetch_st || ref_pos >= fetch_en) goto next;
				off = (int)(ref_pos - fetch_st);
				/* Top-strand-ref-frame read base. SEQ at seq_pos is
				 * t->seq[seq_pos] for r->rev=0 and comp(t->seq[qlen-1-seq_pos])
				 * for r->rev=1; either way it aligns to top_strand_ref[ref_pos]. */
				if (r->rev) {
					int b = kom_nt4_table[(uint8_t)t->seq[qlen - 1 - seq_pos]];
					read_b = b < 4 ? 3 - b : 4;
				} else {
					read_b = kom_nt4_table[(uint8_t)t->seq[seq_pos]];
				}
				if (!xg_is_ga) {
					/* XG=CT: top-strand C; context on top strand. */
					if (ref[off] != 1) goto next;
					c1 = (off + 1 < fetch_en - fetch_st) ? ref[off + 1] : 4;
					c2 = (off + 2 < fetch_en - fetch_st) ? ref[off + 2] : 4;
					if (c1 == 2)        { meth = 'Z'; unmeth = 'z'; }
					else if (c1 == 4)   { meth = 'U'; unmeth = 'u'; }
					else if (c2 == 2)   { meth = 'X'; unmeth = 'x'; }
					else if (c2 == 4)   { meth = 'U'; unmeth = 'u'; }
					else                { meth = 'H'; unmeth = 'h'; }
					if      (read_b == 1) call = meth;
					else if (read_b == 3) call = unmeth;
				} else {
					/* XG=GA: bottom-strand C, top-strand G. Bottom-strand
					 * 5'->3' context base is comp(top[ref_pos-1]); CpG iff
					 * top[ref_pos-1]=C, etc. */
					if (ref[off] != 2) goto next;
					c1 = (ref_pos - 1 >= fetch_st) ? ref[off - 1] : 4;
					c2 = (ref_pos - 2 >= fetch_st) ? ref[off - 2] : 4;
					if (c1 == 1)        { meth = 'Z'; unmeth = 'z'; }
					else if (c1 == 4)   { meth = 'U'; unmeth = 'u'; }
					else if (c2 == 1)   { meth = 'X'; unmeth = 'x'; }
					else if (c2 == 4)   { meth = 'U'; unmeth = 'u'; }
					else                { meth = 'H'; unmeth = 'h'; }
					/* meth bottom C preserved -> comp = G (read_b=2);
					 * unmeth bottom C->T -> comp = A (read_b=0). */
					if      (read_b == 2) call = meth;
					else if (read_b == 0) call = unmeth;
				}
next:
				xm[seq_pos] = call;
				++seq_pos;
				++ref_pos;
			}
		} else if (op == 1) { /* I */
			seq_pos += len;
		} else if (op == 2 || op == 3) { /* D, N */
			ref_pos += len;
		}
	}

	kom_sprintf_lite(s, "\tXR:Z:%s\tXG:Z:%s\tXM:Z:%s", xr, xg, xm);
	free(ref);
	free(xm);
}

void mb_fmt_paf(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx)
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
	if (p->parent == p->id) kom_sprintf_lite(s, "\ts2:i:%d", p->subsc >= 0? p->subsc : 0);
	if (p->p) {
		write_tags(s, p);
		if (p->p->n_cigar > 0) {
			int32_t i;
			kom_sprintf_lite(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i)
				kom_sprintf_lite(s, "%d%c", p->p->cigar[i]>>4, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
		}
		if (p->p->cs) kom_sprintf_lite(s, "\t%s", (char*)&p->p->cigar[p->p->n_cigar]);
	}
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
}

/**************
 * SAM header *
 **************/

static char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

static int sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mb_rg_id, 0, 256);
	if (s == 0) return 0;
	if (strstr(s, "@RG") != s) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = kom_strdup(s);
	mb_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mb_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	kom_sprintf_lite(str, "%s\n", rg_line);
	free(rg_line);
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	if (rg) ret = sam_write_rg_line(str, rg);
	kom_sprintf_lite(str, "@PG\tID:minibwa\tPN:minibwa");
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 1; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	kom_sprintf_lite(str, "\n");
	return ret;
}

/**************
 * SAM output *
 **************/

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kom_roundup64(s->m);
		s->s = kom_realloc(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? kom_comp_table[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

static inline const mb_hit_t *get_sam_pri(int n_hit, const mb_hit_t *hit)
{
	int i;
	for (i = 0; i < n_hit; ++i)
		if (hit[i].sam_pri)
			return &hit[i];
	assert(n_hit == 0);
	return NULL;
}

static void write_sam_cigar(kstring_t *s, int sam_flag, int in_tag, int qlen, const mb_hit_t *r, int64_t opt_flag)
{
	if (r->p == 0) {
		kom_sprintf_lite(s, "*");
	} else {
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag) {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 5 : 4;
			kom_sprintf_lite(s, "\tCG:B:I");
			if (clip_len[0]) kom_sprintf_lite(s, ",%u", clip_len[0]<<4|clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				kom_sprintf_lite(s, ",%u", r->p->cigar[k]);
			if (clip_len[1]) kom_sprintf_lite(s, ",%u", clip_len[1]<<4|clip_char);
		} else {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 'H' : 'S';
			assert(clip_len[0] < qlen && clip_len[1] < qlen);
			if (clip_len[0]) kom_sprintf_lite(s, "%d%c", clip_len[0], clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				kom_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, MB_CIGAR_STR[r->p->cigar[k]&0xf]);
			if (clip_len[1]) kom_sprintf_lite(s, "%d%c", clip_len[1], clip_char);
		}
	}
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx, int32_t mate_qlen)
{
	int flag, n_h = n_hit[seg_idx];
	int this_tid = -1, this_pos = -1;
	const mb_hit_t *h = hit[seg_idx], *r_prev = NULL, *r_next;
	const mb_hit_t *r = n_h > 0 && hit_idx < n_h && hit_idx >= 0? &h[hit_idx] : NULL;

	assert(n_seg == 1 || n_seg == 2);

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1) {
		int next_sid = (seg_idx + 1) % n_seg;
		r_prev = r_next = get_sam_pri(n_hit[next_sid], hit[next_sid]);
	} else r_prev = r_next = NULL;

	// write QNAME and FLAG
	kom_sprintf_lite(s, "%s", t->name);
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0) {
		flag |= 0x4;
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
	}
	if (n_seg > 1) {
		if (r && r->proper_pair) flag |= 0x2;
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	}
	kom_sprintf_lite(s, "\t%d", flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0) {
		if (r_prev) {
			this_tid = r_prev->tid, this_pos = r_prev->ts;
			kom_sprintf_lite(s, "\t%s\t%d\t0\t*", l2b->ctg[this_tid].name, this_pos+1);
		} else kom_sprintf_lite(s, "\t*\t0\t0\t*");
	} else {
		this_tid = r->tid, this_pos = r->ts;
		kom_sprintf_lite(s, "\t%s\t%d\t%d\t", l2b->ctg[r->tid].name, r->ts+1, r->mapq);
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt_flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_tid >= 0 && r_next) {
			if (this_tid == r_next->tid) {
				if (r) {
					int this_pos5 = r->rev? r->te - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->te - 1 : r_next->ts;
					tlen = next_pos5 - this_pos5;
				}
				kom_sprintf_lite(s, "\t=\t");
			} else kom_sprintf_lite(s, "\t%s\t", l2b->ctg[r_next->tid].name);
			kom_sprintf_lite(s, "%d\t", r_next->ts + 1);
		} else if (r_next) { // && this_tid < 0
			kom_sprintf_lite(s, "\t%s\t%d\t", l2b->ctg[r_next->tid].name, r_next->ts + 1);
		} else if (this_tid >= 0) { // && r_next == NULL
			kom_sprintf_lite(s, "\t=\t%d\t", this_pos + 1); // next segment will take r's coordinate
		} else kom_sprintf_lite(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		kom_sprintf_lite(s, "%d\t", tlen);
	} else kom_sprintf_lite(s, "\t*\t0\t0\t");

	// write SEQ and QUAL
	if (r == 0) {
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		kom_sprintf_lite(s, "\t");
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else kom_sprintf_lite(s, "*");
	} else {
		if ((flag & 0x900) == 0 || (opt_flag & MB_F_SUPP_SOFT)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		} else if ((flag & 0x100) && !(opt_flag & MB_F_2ND_SEQ)){
			kom_sprintf_lite(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		}
	}

	// write tags
	if (mb_rg_id[0]) kom_sprintf_lite(s, "\tRG:Z:%s", mb_rg_id);
	if (n_seg > 2) kom_sprintf_lite(s, "\tFI:i:%d", seg_idx);
	if (r) {
		write_tags(s, r);
		if (opt_flag & MB_F_METH) write_meth_tags(s, l2b, t, n_seg, seg_idx, r);
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			kom_sprintf_lite(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt_flag);
			kom_sprintf_lite(s, "\tMQ:i:%d", r_next->mapq);
		}
		if (r->p->cs) kom_sprintf_lite(s, "\t%s", (char*)&r->p->cigar[r->p->n_cigar]);
		if (r->parent == r->id && r->p && n_h > 1 && h && r >= h && r - h < n_h) { // supplementary aln may exist
			int i, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_h; ++i)
				if (i != r - h && h[i].parent == h[i].id && h[i].p)
					++n_sa;
			if (n_sa > 0) {
				kom_sprintf_lite(s, "\tSA:Z:");
				for (i = 0; i < n_h; ++i) {
					const mb_hit_t *q = &h[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0) continue;
					if (q->qe - q->qs < q->te - q->ts) l_M = q->qe - q->qs, l_D = (q->te - q->ts) - l_M;
					else l_M = q->te - q->ts, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					kom_sprintf_lite(s, "%s,%d,%c,", l2b->ctg[q->tid].name, q->ts+1, "+-"[q->rev]);
					if (clip5) kom_sprintf_lite(s, "%dS", clip5);
					if (l_M) kom_sprintf_lite(s, "%dM", l_M);
					if (l_I) kom_sprintf_lite(s, "%dI", l_I);
					if (l_D) kom_sprintf_lite(s, "%dD", l_D);
					if (clip3) kom_sprintf_lite(s, "%dS", clip3);
					kom_sprintf_lite(s, ",%d,%d;", q->mapq, q->blen - q->mlen + q->p->n_ambi);
				}
			}
		}
	}

	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx, int32_t mate_qlen)
{
	if (!(opt_flag & MB_F_PAF))
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt_flag, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt_flag, n_seg, seg_idx);
}
