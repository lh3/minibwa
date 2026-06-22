#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "kommon.h"
#include "l2bit.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

static int64_t l2b_pos2cid(const l2b_t *l2b, int64_t s, int64_t len, int64_t *cst)
{
	int64_t lo = 0, hi = l2b->n_ctg, mid;
	while (lo < hi) {
		const l2b_ctg_t *ctg;
		mid = (lo + hi) / 2;
		ctg = &l2b->ctg[mid];
		if (ctg->off <= s && s < ctg->off + ctg->len) {
			*cst = s - ctg->off;
			return s + len <= ctg->off + ctg->len? mid : -1;
		} else if (s < ctg->off) hi = mid;
		else lo = mid + 1;
	}
	return -1; // s is in no contig
}

int64_t l2b_intv2cid(const l2b_t *l2b, uint64_t st, uint64_t en, int64_t *cst, int *rev)
{
	int64_t s;
	assert(st < en);
	if (en > l2b->tot_len * 2) return -3;
	if (st < l2b->tot_len && l2b->tot_len < en) return -2;
	*rev = (st >= l2b->tot_len);
	s = st < l2b->tot_len? st : l2b->tot_len * 2 - en;
	return l2b_pos2cid(l2b, s, en - st, cst);
}

int64_t l2b_intv2cid_meth(const l2b_t *l2b, uint64_t st, uint64_t en, l2b_meth_t *mt, int64_t *cst, int *rev)
{
	int64_t s, len = en - st;
	int32_t copy;
	uint64_t tot_len = l2b->tot_len;

	assert(st < en);
	if (en > tot_len * 4) return -3;
	copy = st / tot_len; // 0: c2t_f, 1: g2a_f, 2: g2a_r, 3: c2t_r
	*mt = (copy == 0 || copy == 3)? L2B_METH_C2T : L2B_METH_G2A;
	*rev = (copy >= 2);
	s = st - tot_len * copy;
	if (copy >= 2) s = tot_len - len - s; // flip for reverse copies
	return s < 0? -2 : l2b_pos2cid(l2b, s, len, cst);
}

int64_t l2b_getseq(const l2b_t *l2b, int64_t tid, int64_t st, int64_t en, uint8_t *seq)
{
	const l2b_ctg_t *ctg;
	int64_t i, aid;
	int32_t n_ambi;
	if (tid < 0 || tid >= l2b->n_ctg) return -1;
	ctg = &l2b->ctg[tid];
	if (st < 0) st = 0;
	if (en > ctg->len) en = ctg->len;
	st += ctg->off;
	en += ctg->off;
	for (i = st; i < en; ++i)
		seq[i - st] = l2b_get0(l2b, i);
	// retrieve ambiguous bases
	aid = l2b_getambi(l2b, tid, st - ctg->off, en - ctg->off, &n_ambi);
	if (aid >= 0) {
		for (i = 0; i < n_ambi; ++i) {
			const l2b_intv_t *iv = &l2b->ambi[aid + i];
			int64_t s = iv->st, e = iv->en;
			if (s < st) s = st;
			if (e > en) e = en;
			if (s < e) memset(&seq[s - st], 4, e - s);
		}
	}
	return en - st;
}

void l2b_meth_convert(l2b_meth_t mt, int64_t len, uint8_t *seq)
{
	int64_t i;
	if (mt == L2B_METH_C2T) {
		for (i = 0; i < len; ++i)
			if (seq[i] == 1) seq[i] = 3;
	} else if (mt == L2B_METH_G2A) {
		for (i = 0; i < len; ++i)
			if (seq[i] == 2) seq[i] = 0;
	}
}

int64_t l2b_getambi(const l2b_t *l2b, int64_t tid, int64_t st, int64_t en, int32_t *n_ambi)
{
	int64_t g_beg, g_end, lo, hi, mid, i_st, i_en;
	*n_ambi = 0;
	if (tid < 0 || tid >= l2b->n_ctg) return -1;
	if (st < 0) st = 0;
	if (en > l2b->ctg[tid].len) en = l2b->ctg[tid].len;
	if (st >= en) return -1;
	g_beg = l2b->ctg[tid].off + st;
	g_end = l2b->ctg[tid].off + en;

	lo = 0, hi = l2b->n_ambi;
	while (lo < hi) {
		mid = (lo + hi) / 2;
		if (l2b->ambi[mid].en > g_beg) hi = mid;
		else lo = mid + 1;
	}
	i_st = lo;

	lo = i_st, hi = l2b->n_ambi;
	while (lo < hi) {
		mid = (lo + hi) / 2;
		if (l2b->ambi[mid].st >= g_end) hi = mid;
		else lo = mid + 1;
	}
	i_en = lo;

	*n_ambi = i_en - i_st;
	if (*n_ambi == 0) return -1;
	return i_st;
}

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
		} else len = 0;
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
	free(l2b->pac); free(l2b->ambi); free(l2b->mask); free(l2b->ctg); free(l2b);
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
	if (strncmp(magic, L2B_MAGIC, 4) != 0) {
		if (fp != stdin) fclose(fp);
		return 0;
	}
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
	if (fp != stdin) fclose(fp);
	return l2b;
load_failure:
	if (fp != stdin) fclose(fp);
	l2b_destroy(l2b);
	return 0;
}

/******************************
 * Save .paf files for bwtgen *
 ******************************/

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

static inline uint8_t l2b_c2t(uint8_t b) { return b == 1? 3 : b; } // C(1) -> T(3)
static inline uint8_t l2b_g2a(uint8_t b) { return b == 2? 0 : b; } // G(2) -> A(0)

int l2b_save_pac_meth(const char *fn, const l2b_t *l2b, int both_strand)
{
	FILE *fp;
	uint64_t n_pac, x, len;
	int64_t i;
	uint8_t *pac, ct, b;

	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;

	len = l2b->tot_len * 2; // c2t + g2a
	if (both_strand) len *= 2; // forward + reverse
	n_pac = (len + 3) / 4;
	pac = kom_calloc(uint8_t, n_pac);

	// c2t forward
	for (i = 0, x = 0; i < l2b->tot_len; ++i, ++x) {
		b = l2b_c2t(l2b_get0(l2b, i));
		pac[x>>2] |= b << ((~x&3) * 2);
	}
	// g2a forward
	for (i = 0; i < l2b->tot_len; ++i, ++x) {
		b = l2b_g2a(l2b_get0(l2b, i));
		pac[x>>2] |= b << ((~x&3) * 2);
	}
	if (both_strand) {
		// g2a reverse (reverse complement of g2a converted sequence)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x) {
			b = 3 - l2b_g2a(l2b_get0(l2b, i));
			pac[x>>2] |= b << ((~x&3) * 2);
		}
		// c2t reverse (reverse complement of c2t converted sequence)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x) {
			b = 3 - l2b_c2t(l2b_get0(l2b, i));
			pac[x>>2] |= b << ((~x&3) * 2);
		}
	}

	fwrite(pac, 1, (x>>2) + ((x&3) == 0? 0 : 1), fp);
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
