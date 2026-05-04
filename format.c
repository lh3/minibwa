#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const mb_meth_cmap_t *cmap, const char *rg, const char *ver, int argc, char *argv[])
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (cmap) { // --meth: one @SQ per consolidated contig
		for (i = 0; i < cmap->n_output; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n",
				cmap->output_names[i], (long)cmap->output_lens[i]);
	} else if (idx) {
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	}
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

// Append `comment` to `s` (preceded by '\t'). If `filter`, drop segments that don't conform to
// SAM aux-tag shape `XX:T:...` (T one of A/i/f/Z/H/B). --meth uses filter=1 because YS:Z/YC:Z
// are stuffed alongside the user's FASTQ comment; non-meth -y is verbatim.
static void sam_write_comment(kstring_t *s, const char *comment, int filter)
{
	const char *p = comment, *seg_start;
	if (!filter) { kom_sprintf_lite(s, "\t%s", comment); return; }
	while (*p) {
		int n, ok;
		seg_start = p;
		while (*p && *p != '\t') ++p;
		n = (int)(p - seg_start);
		ok = n >= 5
			&& ((seg_start[0] >= 'A' && seg_start[0] <= 'Z') || (seg_start[0] >= 'a' && seg_start[0] <= 'z'))
			&& ((seg_start[1] >= 'A' && seg_start[1] <= 'Z') || (seg_start[1] >= 'a' && seg_start[1] <= 'z') || (seg_start[1] >= '0' && seg_start[1] <= '9'))
			&& seg_start[2] == ':'
			&& (seg_start[3] == 'A' || seg_start[3] == 'i' || seg_start[3] == 'f' || seg_start[3] == 'Z' || seg_start[3] == 'H' || seg_start[3] == 'B')
			&& seg_start[4] == ':';
		if (ok) {
			str_enlarge(s, n + 1);
			s->s[s->l++] = '\t';
			memcpy(&s->s[s->l], seg_start, n);
			s->l += n;
			s->s[s->l] = 0;
		}
		if (*p == '\t') ++p;
	}
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

// like sam_write_sq, but uppercases ACGT* (consumers expect uppercase under --meth)
static void sam_write_sq_upper(kstring_t *s, char *seq, int l, int rev, int comp)
{
	int i;
	str_enlarge(s, l);
	if (rev) {
		for (i = 0; i < l; ++i) {
			int c = (unsigned char)seq[l - 1 - i];
			char out = (c < 128 && comp) ? (char)kom_comp_table[c] : (char)c;
			if (out >= 'a' && out <= 'z') out = (char)(out - 32);
			s->s[s->l + i] = out;
		}
	} else {
		for (i = 0; i < l; ++i) {
			char out = seq[i];
			if (out >= 'a' && out <= 'z') out = (char)(out - 32);
			s->s[s->l + i] = out;
		}
	}
	s->l += l;
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

// Pre-c2t bases stashed as YS:Z prefix of t->comment by mb_meth_ingest. NULL on missing/malformed.
static inline const char *meth_orig_seq(const mb_bseq1_t *t)
{
	if (t == 0 || t->comment == 0 || t->l_seq <= 0) return 0;
	if (memcmp(t->comment, MB_METH_YS_PREFIX, 5) != 0) return 0;
	if (strlen(t->comment) < (size_t)5 + (size_t)t->l_seq) return 0;
	return t->comment + 5;
}

static inline const char *sam_ctg_name(const l2b_t *l2b, const mb_meth_cmap_t *cmap, int64_t tid)
{
	if (cmap && tid >= 0 && tid < cmap->n_internal) {
		int out = cmap->out_tid[tid];
		if (out >= 0 && out < cmap->n_output) return cmap->output_names[out];
	}
	return l2b->ctg[tid].name;
}

// True if internal tids a, b collapse to the same output contig.
static inline int sam_same_output_tid(const mb_meth_cmap_t *cmap, int64_t a, int64_t b)
{
	if (a == b) return 1;
	if (cmap && a >= 0 && b >= 0 && a < cmap->n_internal && b < cmap->n_internal)
		return cmap->out_tid[a] == cmap->out_tid[b];
	return 0;
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_meth_cmap_t *cmap, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx, int32_t mate_qlen)
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
		if (r->qc_fail) flag |= 0x200; // set by mb_meth_apply_qc (chimera / set-as-failed)
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
			kom_sprintf_lite(s, "\t%s\t%d\t0\t*", sam_ctg_name(l2b, cmap, this_tid), this_pos+1);
		} else kom_sprintf_lite(s, "\t*\t0\t0\t*");
	} else {
		this_tid = r->tid, this_pos = r->ts;
		kom_sprintf_lite(s, "\t%s\t%d\t%d\t", sam_ctg_name(l2b, cmap, r->tid), r->ts+1, r->mapq);
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt_flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_tid >= 0 && r_next) {
			// f<chr> + r<chr> mates collapse onto one output contig; require both CIGARs for TLEN
			if (sam_same_output_tid(cmap, this_tid, r_next->tid)) {
				if (r && r->p && r_next->p) {
					int this_pos5 = r->rev? r->te - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->te - 1 : r_next->ts;
					tlen = next_pos5 - this_pos5;
				}
				kom_sprintf_lite(s, "\t=\t");
			} else kom_sprintf_lite(s, "\t%s\t", sam_ctg_name(l2b, cmap, r_next->tid));
			kom_sprintf_lite(s, "%d\t", r_next->ts + 1);
		} else if (r_next) { // && this_tid < 0
			kom_sprintf_lite(s, "\t%s\t%d\t", sam_ctg_name(l2b, cmap, r_next->tid), r_next->ts + 1);
		} else if (this_tid >= 0) { // && r_next == NULL
			kom_sprintf_lite(s, "\t=\t%d\t", this_pos + 1); // next segment will take r's coordinate
		} else kom_sprintf_lite(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		kom_sprintf_lite(s, "%d\t", tlen);
	} else kom_sprintf_lite(s, "\t*\t0\t0\t");

	// --meth: restore pre-c2t bases from YS:Z so callers see real C/T content (else c2t-projected)
	const char *seq_src = t->seq;
	int meth_seq = 0;
	if (opt_flag & MB_F_METH) {
		const char *o = meth_orig_seq(t);
		if (o) seq_src = o, meth_seq = 1;
	}
	void (*write_sq)(kstring_t *, char *, int, int, int) =
		meth_seq? sam_write_sq_upper : sam_write_sq;

	// write SEQ and QUAL
	if (r == 0) {
		write_sq(s, (char *)seq_src, t->l_seq, 0, 0);
		kom_sprintf_lite(s, "\t");
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else kom_sprintf_lite(s, "*");
	} else {
		if ((flag & 0x900) == 0 || (opt_flag & MB_F_SUPP_SOFT)) {
			write_sq(s, (char *)seq_src, t->l_seq, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		} else if ((flag & 0x100) && !(opt_flag & MB_F_2ND_SEQ)){
			kom_sprintf_lite(s, "*\t*");
		} else {
			write_sq(s, (char *)seq_src + r->qs, r->qe - r->qs, r->rev, r->rev);
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
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			kom_sprintf_lite(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt_flag);
			kom_sprintf_lite(s, "\tMQ:i:%d", r_next->mapq);
		}
		if (r->p->cs) kom_sprintf_lite(s, "\t%s", (char*)&r->p->cigar[r->p->n_cigar]);
		// YD:Z direction tag from cmap (or contig-name prefix as fallback when cmap is absent)
		if ((opt_flag & MB_F_METH) && r->tid >= 0) {
			char dir = 0;
			if (cmap && r->tid < cmap->n_internal) dir = cmap->direction[r->tid];
			if (dir == 0) {
				const char *cname = l2b->ctg[r->tid].name;
				if (cname && (cname[0] == 'f' || cname[0] == 'r')) dir = cname[0];
			}
			if (dir != 0) kom_sprintf_lite(s, "\tYD:Z:%c", dir);
		}
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
					kom_sprintf_lite(s, "%s,%d,%c,", sam_ctg_name(l2b, cmap, q->tid), q->ts+1, "+-"[q->rev]);
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
		sam_write_comment(s, t->comment, !!(opt_flag & MB_F_METH));
	kom_sprintf_lite(s, "\n");
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_meth_cmap_t *cmap, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx, int32_t mate_qlen)
{
	if (opt_flag & MB_F_SAM)
		mb_fmt_sam(km, s, l2b, cmap, t, n_seg, n_hit, hit, hit_idx, opt_flag, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt_flag, n_seg, seg_idx);
}
