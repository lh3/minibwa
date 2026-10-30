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

/*************
 * Utilities *
 *************/

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

/**************
 * SAM header *
 **************/

char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
			else if (*p == '\0') break;
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

// get the value of a two-letter tag on a header line; return its length
static int hdr_get_tag(const char *line, int len, const char *tag, char *out, int max)
{
	int i, l = 0;
	for (i = 3; i + 3 < len; ++i) {
		if (line[i] != '\t' || line[i+1] != tag[0] || line[i+2] != tag[1] || line[i+3] != ':') continue;
		for (i += 4; i < len && line[i] != '\t' && l < max - 1; ++i)
			out[l++] = line[i];
		break;
	}
	out[l] = 0;
	return l;
}

// 1 if str already holds line[0,len)
static int hdr_has_line(const kstring_t *str, const char *line, int len)
{
	size_t i;
	for (i = 0; i + len < str->l; ++i)
		if ((i == 0 || str->s[i-1] == '\n') && str->s[i+len] == '\n' && strncmp(&str->s[i], line, len) == 0)
			return 1;
	return 0;
}

// 1 if str already holds a line of the given type (e.g. "@RG") with this ID
static int hdr_has_id(const kstring_t *str, const char *type, const char *id)
{
	size_t i;
	int l_id = strlen(id);
	for (i = 0; i + 4 < str->l; ++i) {
		char buf[256];
		size_t j;
		if (!(i == 0 || str->s[i-1] == '\n')) continue;
		if (strncmp(&str->s[i], type, 3) != 0 || str->s[i+3] != '\t') continue;
		for (j = i; j < str->l && str->s[j] != '\n'; ++j) {}
		if (hdr_get_tag(&str->s[i], j - i, "ID", buf, sizeof(buf)) == l_id && strcmp(buf, id) == 0)
			return 1;
	}
	return 0;
}

// carry @RG (deduplicated by ID), @PG and @CO of an input header over to rg[]
// and other[]; @HD and @SQ are dropped as minibwa writes its own. last_pg takes
// the ID of the last @PG, for PP chaining.
static void sam_copy_hdr(kstring_t *rg, kstring_t *other, const char *text, char *last_pg, int max_pg)
{
	const char *p, *q;
	for (p = text; *p; p = *q? q + 1 : q) {
		int len;
		for (q = p; *q && *q != '\n'; ++q) {}
		len = q - p;
		while (len > 0 && p[len-1] == '\r') --len; // tolerate CRLF
		if (len < 3 || p[0] != '@') continue;
		if (strncmp(p, "@HD", 3) == 0 || strncmp(p, "@SQ", 3) == 0) continue;
		if (strncmp(p, "@RG", 3) == 0) {
			char id[256];
			if (hdr_get_tag(p, len, "ID", id, sizeof(id)) == 0) {
				if (kom_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m an @RG line in the input header has no ID; dropped.\033[0m\n");
				continue;
			}
			if (hdr_has_id(rg, "@RG", id)) continue;
			str_copy(rg, p, p + len);
			kom_sprintf_lite(rg, "\n");
		} else {
			if (hdr_has_line(other, p, len)) continue;
			str_copy(other, p, p + len);
			kom_sprintf_lite(other, "\n");
			if (strncmp(p, "@PG", 3) == 0) {
				char id[256];
				if (hdr_get_tag(p, len, "ID", id, sizeof(id)) > 0) // keep the chain if an @PG has no ID
					snprintf(last_pg, max_pg, "%s", id);
			}
		}
	}
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[], int32_t n_hdr, const char *const *hdr_text)
{
	int i, ret = 0, n_dup = 0;
	kstring_t in_rg = {0,0,0}, in_other = {0,0,0};
	char last_pg[256], pg_id[256];

	last_pg[0] = 0;
	strcpy(pg_id, "minibwa");
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	for (i = 0; i < n_hdr; ++i)
		if (hdr_text && hdr_text[i])
			sam_copy_hdr(&in_rg, &in_other, hdr_text[i], last_pg, sizeof(last_pg));
	if (rg) in_rg.l = 0; // -R applies to every record, so the input @RG lines describe nothing
	if (in_rg.l > 0) kom_sprintf_lite(str, "%s", in_rg.s);
	if (rg) ret = sam_write_rg_line(str, rg);
	if (in_other.l > 0) kom_sprintf_lite(str, "%s", in_other.s);
	while (hdr_has_id(&in_other, "@PG", pg_id)) // @PG IDs must be unique, e.g. on re-aligning minibwa output
		snprintf(pg_id, sizeof(pg_id), "minibwa.%d", ++n_dup);
	free(in_rg.s); free(in_other.s);
	kom_sprintf_lite(str, "@PG\tID:%s\tPN:minibwa", pg_id);
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (last_pg[0]) kom_sprintf_lite(str, "\tPP:%s", last_pg);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 0; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	kom_sprintf_lite(str, "\n");
	return ret;
}

/**************
 * SAM output *
 **************/

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

// append the tags carried over from a BAM input, skipping any minibwa has
// already written: SAM forbids a repeated TAG on one line. bam_aux_drop[] in
// bam.c removes the expected ones; this keeps the output valid regardless.
#define MB_MAX_OWN_TAG 64 // minibwa writes ~14; this only has to be an upper bound

static inline uint16_t sam_tag_key(const char *tag) // a two-letter tag name packed into an int
{
	return (uint16_t)((uint8_t)tag[0]<<8 | (uint8_t)tag[1]);
}

static inline int sam_tag_seen(const uint16_t *key, int32_t n_key, uint16_t t)
{
	int32_t i;
	for (i = 0; i < n_key; ++i)
		if (key[i] == t) return 1;
	return 0;
}

static void sam_write_aux(kstring_t *s, size_t tag_st, const char *aux)
{
	static int warned = 0;
	uint16_t key[MB_MAX_OWN_TAG];
	int32_t n_key = 0;
	size_t i;
	const char *p, *q;

	for (i = tag_st; i + 3 < s->l; ++i) // collect once; rescanning s per tag would be quadratic
		if (s->s[i] == '\t' && s->s[i+3] == ':' && n_key < MB_MAX_OWN_TAG)
			key[n_key++] = sam_tag_key(&s->s[i+1]);
	for (p = aux; *p == '\t'; p = q) {
		for (q = p + 1; *q && *q != '\t'; ++q) {}
		if (q - p < 4 || p[3] != ':') continue; // not a \tXX:T:VALUE field
		if (!sam_tag_seen(key, n_key, sam_tag_key(p + 1))) str_copy(s, p, q);
		else if (!warned && kom_verbose >= 2) {
			warned = 1;
			fprintf(stderr, "[WARNING]\033[1;31m input tag %c%c replaced by minibwa's own. Reported once.\033[0m\n", p[1], p[2]);
		}
	}
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

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	int flag, n_h = n_hit[seg_idx];
	size_t tag_st;
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
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	} else if (t->flag & 0x1) flag |= 0x1 | 0x8; // a BAM segment whose mate is absent from this run
	if (t->flag & 0xc0) flag |= t->flag & 0xc0; // the input says which segment this is
	else if (n_seg > 1) {
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
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
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt->flag);
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
		} else if (this_tid >= 0) { // && r_next == NULL. In this case, mate is unmapped and placed at this segment's primary alignment
			const mb_hit_t *r_pri;
			r_pri = get_sam_pri(n_h, h);
			if (r->tid == r_pri->tid)
				kom_sprintf_lite(s, "\t=\t%d\t", r_pri->ts + 1);
			else
				kom_sprintf_lite(s, "\t%s\t%d\t", l2b->ctg[r_pri->tid].name, r_pri->ts + 1);
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
		if ((flag & 0x900) == 0 || (opt->flag & MB_F_SUPP_SOFT)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		} else if ((flag & 0x100) && !(opt->flag & MB_F_2ND_SEQ)){
			kom_sprintf_lite(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		}
	}

	// write tags
	tag_st = s->l;
	if (mb_rg_id[0]) kom_sprintf_lite(s, "\tRG:Z:%s", mb_rg_id); // -R replaces a record's own RG
	if (n_seg > 2) kom_sprintf_lite(s, "\tFI:i:%d", seg_idx);
	if (r) {
		write_tags(s, r);
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			kom_sprintf_lite(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt->flag);
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
			if (opt->xa_max > 0) {
				int i, n_xa = 0;
				for (i = 0; i < n_h; ++i)
					if (i != r - h && h[i].parent == r - h && h[i].p->dp_max >= (double)opt->out_s * r->p->dp_max)
						++n_xa;
				if (n_xa > 0) kom_sprintf_lite(s, "\tn2:i:%d", n_xa);
				if (n_xa > 0 && n_xa <= opt->xa_max) {
					kom_sprintf_lite(s, "\tXA:Z:");
					for (i = 0; i < n_h; ++i) {
						const mb_hit_t *q = &h[i];
						if (i != r - h && q->parent == r - h && q->p->dp_max >= (double)opt->out_s * r->p->dp_max) {
							kom_sprintf_lite(s, "%s,%c%d,", l2b->ctg[q->tid].name, "+-"[q->rev], q->ts+1);
							write_sam_cigar(s, 0, 0, t->l_seq, q, opt->flag);
							kom_sprintf_lite(s, ",%d;", q->blen - q->mlen + q->p->n_ambi);
						}
					}
				}
			}
		}
	}

	if (t->aux) sam_write_aux(s, tag_st, t->aux); // tags carried over from a BAM input

	if ((opt->flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	if (!(opt->flag & MB_F_PAF))
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt->flag, n_seg, seg_idx);
}
