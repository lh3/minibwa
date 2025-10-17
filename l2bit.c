#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "kommon.h"
#include "l2bit.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

static void l2b_format_seq(uint64_t len, char *seq, uint64_t *rng)
{
	uint64_t i;
	for (i = 0; i < len; ++i) {
		int b = (uint8_t)seq[i];
		uint8_t c = kom_nt4_table[b];
		if (c > 4) c = 4;
		if (c == 4) c |= kom_splitmix64(rng) & 3;
		if (b < 'A' || b > 'Z') c |= 1<<3;
		seq[i] = c;
	}
}

static void l2b_add_seq(l2b_t *l2b, uint64_t len, const char *seq, const char *name, const char *comm, uint64_t *rng)
{
	uint64_t i, ambi_len, mask_len, off, m_pac_old;
	l2b_ctg_t *ctg;

	off = l2b->tot_len;
	kom_grow(l2b_ctg_t, l2b->ctg, l2b->n_ctg, l2b->m_ctg);
	ctg = &l2b->ctg[l2b->n_ctg++];
	ctg->name = kom_strdup(name);
	ctg->comm = comm? kom_strdup(comm) : 0;
	ctg->len = len;
	ctg->off = l2b->tot_len;
	l2b->tot_len += len;

	m_pac_old = l2b->m_pac;
	l2b->n_pac = (l2b->tot_len + 31) / 32;
	kom_grow(uint64_t, l2b->pac, l2b->n_pac, l2b->m_pac);
	if (m_pac_old < l2b->m_pac) // zero out newly allocated part
		memset(&l2b->pac[m_pac_old], 0, (l2b->m_pac - m_pac_old) * 8);

	for (i = 0, ambi_len = mask_len = 0; i < len; ++i) {
		uint64_t c = (uint8_t)seq[i], x = off + i;
		if (c & 1<<3) { // soft-masked base
			++mask_len;
		} else if (mask_len > 0) {
			kom_grow(l2b_intv_t, l2b->mask, l2b->n_mask, l2b->m_mask);
			l2b->mask[l2b->n_mask].st = x - mask_len;
			l2b->mask[l2b->n_mask].en = x;
			l2b->n_mask++;
			mask_len = 0;
		}
		if (c & 1<<2) { // ambiguous base
			++ambi_len;
		} else if (ambi_len > 0) {
			kom_grow(l2b_intv_t, l2b->ambi, l2b->n_ambi, l2b->m_ambi);
			l2b->ambi[l2b->n_ambi].st = x - ambi_len;
			l2b->ambi[l2b->n_ambi].en = x;
			l2b->n_ambi++;
			ambi_len = 0;
		}
		l2b->pac[x>>5] |= (c&3) << (x&0x1f)*2;
	}
}

static void l2b_collate_str(l2b_t *l2b)
{
	uint64_t i, tot_name = 0, tot_comm = 0;
	char *p_name, *p_comm;
	if (l2b->cat_name || l2b->cat_comm) return;
	for (i = 0; i < l2b->n_ctg; ++i) {
		l2b_ctg_t *ctg = &l2b->ctg[i];
		tot_name += strlen(ctg->name) + 1;
		tot_comm += ctg->comm? strlen(ctg->comm) + 1 : 1;
	}
	p_name = l2b->cat_name = kom_calloc(char, tot_name);
	p_comm = l2b->cat_comm = kom_calloc(char, tot_comm);
	for (i = 0; i < l2b->n_ctg; ++i) {
		l2b_ctg_t *ctg = &l2b->ctg[i];
		uint64_t len;
		len = strlen(ctg->name);
		memcpy(p_name, ctg->name, len + 1);
		free(ctg->name);
		ctg->name = p_name;
		p_name += len + 1;
		if (ctg->comm) {
			len = strlen(ctg->comm);
			memcpy(p_comm, ctg->comm, len + 1);
			free(ctg->comm);
			ctg->comm = p_comm;
		} else len = 1;
		p_comm += len + 1;
	}
}

l2b_t *l2b_import(const char *fn, uint64_t seed)
{
	gzFile fp;
	kseq_t *ks;
	l2b_t *l2b;
	uint64_t rng = seed;

	fp = fn == 0 || strcmp(fn, "-") == 0? gzdopen(0, "r") : gzopen(fn, "r");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	l2b = kom_calloc(l2b_t, 1);
	while (kseq_read(ks) >= 0) {
		l2b_format_seq(ks->seq.l, ks->seq.s, &rng);
		l2b_add_seq(l2b, ks->seq.l, ks->seq.s, ks->name.s, ks->comment.l? ks->comment.s : 0, &rng);
	}
	kseq_destroy(ks);
	gzclose(fp);
	l2b_collate_str(l2b);
	return l2b;
}

void l2b_destroy(l2b_t *l2b)
{
	free(l2b->cat_name); free(l2b->cat_comm);
	free(l2b->pac); free(l2b->ambi); free(l2b->ctg); free(l2b);
}

int l2b_save(const char *fn, const l2b_t *l2b)
{
	FILE *fp;
	uint64_t i, len_name = 0, len_comm = 0;
	uint32_t dummy = 0;
	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;
	for (i = 0; i < l2b->n_ctg; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[i];
		len_name += strlen(ctg->name) + 1;
		len_comm += ctg->comm? strlen(ctg->comm) + 1 : 1;
	}
	fwrite(L2B_MAGIC, 1, 4, fp);
	fwrite(&dummy, 4, 1, fp);
	fwrite(&l2b->n_ctg, 8, 1, fp);
	fwrite(&l2b->tot_len, 8, 1, fp);
	fwrite(&l2b->n_ambi, 8, 1, fp);
	fwrite(&l2b->n_mask, 8, 1, fp);
	fwrite(&len_name, 8, 1, fp);
	fwrite(&len_comm, 8, 1, fp);
	fwrite(&l2b->n_pac, 8, 1, fp);
	for (i = 0; i < l2b->n_ctg; ++i)
		fwrite(&l2b->ctg[i].len, 8, 1, fp);
	fwrite(l2b->ambi, 16, l2b->n_ambi, fp);
	fwrite(l2b->mask, 16, l2b->n_mask, fp);
	fwrite(l2b->pac, 8, l2b->n_pac, fp);
	fwrite(l2b->cat_name, 1, len_name, fp); // put strings at the end to make sure uint64_t are all aligned
	fwrite(l2b->cat_comm, 1, len_comm, fp);
	fclose(fp);
	return 0;
}

l2b_t *l2b_load(const char *fn)
{
	FILE *fp;
	char magic[4], *p_name, *p_comm;
	uint32_t dummy;
	uint64_t off, i, len_name, len_comm;
	l2b_t *l2b;
	fp = fn == 0 || strcmp(fn, "-") == 0? stdin : fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, L2B_MAGIC, 4) != 0) return 0;
	l2b = kom_calloc(l2b_t, 1);
	fread(&dummy, 4, 1, fp);
	fread(&l2b->n_ctg, 8, 1, fp);
	fread(&l2b->tot_len, 8, 1, fp);
	fread(&l2b->n_ambi, 8, 1, fp);
	fread(&l2b->n_mask, 8, 1, fp);
	fread(&len_name, 8, 1, fp);
	fread(&len_comm, 8, 1, fp);
	fread(&l2b->n_pac, 8, 1, fp);
	l2b->ctg = kom_calloc(l2b_ctg_t, l2b->n_ctg);
	for (i = 0, off = 0; i < l2b->n_ctg; ++i) { // read contig lengths
		fread(&l2b->ctg[i].len, 8, 1, fp);
		l2b->ctg[i].off = off;
		off += l2b->ctg[i].len;
	}
	if (off != l2b->tot_len) goto load_failure;
	l2b->ambi = kom_malloc(l2b_intv_t, l2b->n_ambi);
	l2b->mask = kom_malloc(l2b_intv_t, l2b->n_mask);
	l2b->pac = kom_malloc(uint64_t, l2b->n_pac);
	l2b->cat_name = kom_malloc(char, len_name);
	l2b->cat_comm = kom_malloc(char, len_comm);
	fread(l2b->ambi, 16, l2b->n_ambi, fp);
	fread(l2b->mask, 16, l2b->n_mask, fp);
	fread(l2b->pac, 8, l2b->n_pac, fp);
	fread(l2b->cat_name, 1, len_name, fp);
	fread(l2b->cat_comm, 1, len_comm, fp);
	p_name = l2b->cat_name, p_comm = l2b->cat_comm;
	for (i = 0; i < l2b->n_ctg; ++i) { // synchronize contig names and comments
		l2b_ctg_t *ctg = &l2b->ctg[i];
		ctg->name = p_name;
		p_name += strlen(p_name) + 1;
		ctg->comm = *p_comm? p_comm : 0;
		p_comm += *p_comm? strlen(p_comm) + 1 : 1;
	}
	if (p_name - l2b->cat_name != len_name || p_comm - l2b->cat_comm != len_comm) goto load_failure;
	return l2b;
load_failure:
	l2b_destroy(l2b);
	return 0;
}

int l2b_save_pac(const char *fn, const l2b_t *l2b, int both_strand)
{
	FILE *fp;
	uint64_t n_pac, x;
	int64_t i;
	uint8_t *pac, ct;

	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;

	// fill pac[]
	n_pac = ((both_strand? l2b->tot_len * 2 : l2b->tot_len) + 3) / 4;
	pac = kom_calloc(uint8_t, n_pac);
	for (i = 0, x = 0; i < l2b->tot_len; ++i, ++x)
		pac[x>>2] |= l2b_get0(l2b, i) << (~x&3) * 2;
	if (both_strand)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x)
			pac[x>>2] |= (3 - l2b_get0(l2b, i)) << (~x&3) * 2;

	// write pac
	fwrite(pac, 1, (x>>2) + ((x&3) == 0? 0 : 1), fp);
	// the following codes make the pac file size always (x/4+1+1)
	if (x % 4 == 0) {
		ct = 0;
		fwrite(&ct, 1, 1, fp);
	}
	ct = x % 4;
	fwrite(&ct, 1, 1, fp);

	fclose(fp);
	free(pac);
	return 0;
}
