#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#include <zlib.h>
#include "kommon.h"
#include "mbpriv.h"
#include "bseq.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef struct {
	int32_t n_fp, n_threads;
	int64_t n_base, n_seq;
	int64_t mb_size; // mini-batch size
	const mb_opt_t *opt;
	mb_bseq_file_t **fp;
	const mb_idx_t *idx;
	FILE *fp_out;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int32_t n_seq, n_frag, n_sb, n_pe;
	int32_t *n_hit, *seg_off, *seg_cnt;
	int32_t *sb_off, *sb_cnt;
	mb_pestat_t pes[4];
	mb_bseq1_t *seq;
	mb_hit_t **hit;
	mb_tbuf_t **tbuf;
} step_t;

static void worker_for_se_batch(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	const mb_opt_t *opt = s->p->opt;
	const mb_idx_t *idx = s->p->idx;
	mb_tbuf_t *b = s->tbuf[tid];
	void *km;
	int64_t tot;
	int32_t n, j, k, l, p, *len;
	uint8_t **seq, *buf;
	mb_sai_v *sai;

	km = mb_tbuf_km(b);
	for (k = 0, n = 0, tot = 0; k < s->sb_cnt[i]; ++k) {
		n += s->seg_cnt[s->sb_off[i] + k];
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j)
			tot += s->seq[off + j].l_seq;
	}
	buf = Kmalloc(km, uint8_t, tot);
	len = Kmalloc(km, int32_t, n);
	seq = Kmalloc(km, uint8_t*, n);
	sai = Kcalloc(km, mb_sai_v, n);
	for (k = p = 0, tot = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j) {
			const mb_bseq1_t *t = &s->seq[off + j];
			l2b_meth_t mt = !idx->is_meth? L2B_METH_NONE : (j&1) == 0? L2B_METH_C2T : L2B_METH_G2A;
			len[p] = t->l_seq;
			seq[p] = &buf[tot], tot += t->l_seq;
			for (l = 0; l < t->l_seq; ++l)
				seq[p][l] = kom_nt4_table[(uint8_t)t->seq[l]];
			if (mt != L2B_METH_NONE) l2b_meth_convert(mt, len[p], seq[p]);
			++p;
		}
	}
	assert(p == n);
	mb_seed_intv_batch(km, idx->bwt, n, len, seq, opt->min_len, opt->max_sub_occ, sai);
	kfree(km, seq);
	kfree(km, len);
	kfree(km, buf);

	for (k = p = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j) {
			const mb_bseq1_t *t = &s->seq[off + j];
			mb_opt_t opt_adap;
			l2b_meth_t mt = !idx->is_meth? L2B_METH_NONE : (j&1) == 0? L2B_METH_C2T : L2B_METH_G2A;
			mb_opt_adap(opt, t->l_seq, &opt_adap);
			s->hit[off+j] = mb_map_sai(&opt_adap, idx, t->l_seq, t->seq, mt, &sai[p], &s->n_hit[off+j], b, t->name);
			++p;
		}
	}
	kfree(km, sai);
	mb_tbuf_reset(b, opt->cap_kalloc);
}

static void worker_for_pe(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	mb_tbuf_t *b = s->tbuf[tid];
	int32_t off, r, len[2];
	char *seq[2];
	if (s->seg_cnt[i] != 2) return;
	off = s->seg_off[i];
	if (kom_dbg_flag & MB_DBG_QNAME) fprintf(stderr, "QP\t%s\t%d\n", s->seq[off].name, tid);
	for (r = 0; r < 2; ++r) {
		seq[r] = s->seq[off + r].seq;
		len[r] = s->seq[off + r].l_seq;
	}
	mb_pair(mb_tbuf_km(b), s->p->opt, s->p->idx->l2b, &s->n_hit[off], &s->hit[off], s->pes, len, seq);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	const int min_read_cnt = 40000; // should be smaller than opt->mb_size / (read_len * 2) for typical paired-end reads
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
	const mb_opt_t *opt = p->opt;
    if (step == 0) { // step 0: read sequences
		int with_qual = !(opt->flag & MB_F_PAF);
		int with_comment = !!(opt->flag & MB_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(opt->flag & MB_F_PE));
        step_t *s;
        s = kom_calloc(step_t, 1);
		if (p->n_fp > 1) s->seq = mb_bseq_read_frag(p->n_fp, p->fp, p->mb_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mb_bseq_read(p->fp[0], p->mb_size, with_qual, with_comment, frag_mode, min_read_cnt, opt->max_mb_size, &s->n_seq);
		if (s->seq) {
			int32_t sb_len, sb_off;
			s->p = p;
			for (i = 0; i < 4; ++i) s->pes[i].failed = 1;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].id = p->n_seq++;
			s->tbuf = kom_calloc(mb_tbuf_t*, opt->n_thread);
			for (i = 0; i < opt->n_thread; ++i)
				s->tbuf[i] = mb_tbuf_init(opt->flag&MB_F_NO_KALLOC);
			s->n_hit = kom_calloc(int32_t, 5 * s->n_seq); // maybe over allocation as the following 4 members may not need this much
			s->seg_off = s->n_hit   + s->n_seq;
			s->seg_cnt = s->seg_off + s->n_seq;
			s->sb_off  = s->seg_cnt + s->n_seq;
			s->sb_cnt  = s->sb_off  + s->n_seq;
			s->hit = kom_calloc(mb_hit_t*, s->n_seq);
			// set seg_cnt[] and seg_off[]
			for (i = 1, j = 0; i <= s->n_seq; ++i) {
				if (i == s->n_seq || !frag_mode || !mb_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					assert(i - j <= 2);
					s->seg_cnt[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					if (i - j == 2) s->n_pe++;
					j = i;
				}
			}
			if (s->n_pe > 0) { // trim trailing "/[12]" for paired-end reads
				for (i = 0; i < s->n_frag; ++i) {
					int32_t j0, j1, l0, l1;
					if (s->seg_cnt[i] != 2) continue;
					j0 = s->seg_off[i] + 0, l0 = strlen(s->seq[j0].name);
					j1 = s->seg_off[i] + 1, l1 = strlen(s->seq[j1].name);
					if (l0 < 3 || l0 != l1) continue;
					if (s->seq[j0].name[l0 - 1] != s->seq[j1].name[l1 - 1] && s->seq[j0].name[l0 - 2] == '/')
						s->seq[j0].name[l0 - 2] = s->seq[j1].name[l1 - 2] = 0; // truncate
				}
			}
			// set sb_cnt[] and sb_off[]
			for (i = 0, sb_len = sb_off = 0; i < s->n_frag; ++i) {
				if (sb_len >= opt->sb_len || i - sb_off >= opt->sb_seq) {
					s->sb_off[s->n_sb] = sb_off;
					s->sb_cnt[s->n_sb++] = i - sb_off;
					sb_len = 0;
					sb_off = i;
				}
				for (j = 0; j < s->seg_cnt[i]; ++j)
					sb_len += s->seq[s->seg_off[i] + j].l_seq;
			}
			s->sb_off[s->n_sb] = sb_off;
			s->sb_cnt[s->n_sb++] = i - sb_off;
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		step_t *s = (step_t*)in;
		kt_for(opt->n_thread, worker_for_se_batch, in, s->n_sb);
		if ((opt->flag & MB_F_PE) && s->n_frag < s->n_seq && !(opt->flag & MB_F_NO_PAIRING)) { // PE mode
			if ((opt->flag & MB_F_PE_PREDEF) || s->n_pe < 20) { // use predefined PE stats
				s->pes[1].failed = 0;
				s->pes[1].avg = opt->pe_avg, s->pes[1].std = opt->pe_std;
				s->pes[1].lo = opt->pe_lo, s->pes[1].hi = opt->pe_hi;
			} else { // estimate PE stats from data
				void *km;
				km = opt->flag & MB_F_NO_KALLOC? 0 : km_init();
				mb_pestat(km, opt, s->n_frag, s->seg_off, s->seg_cnt, s->n_hit, s->hit, s->pes);
				if (km) km_destroy(km);
			}
			kt_for(opt->n_thread, worker_for_pe, in, s->n_frag);
		}
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mb_idx_t *idx = p->idx;
		kstring_t out = {0,0,0};
		int64_t tot_len = 0;

		for (i = 0; i < opt->n_thread; ++i)
			mb_tbuf_destroy(s->tbuf[i]);
		free(s->tbuf);
		if (!(opt->flag & MB_F_NO_KALLOC)) km = km_init();

		for (k = 0; k < s->n_frag; ++k) {
			int32_t seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->seg_cnt[k];
			out.l = 0;
			for (i = seg_st; i < seg_en; ++i) {
				mb_bseq1_t *t = &s->seq[i];
				int32_t mate_qlen = 0; // mate's l_seq for MC:Z; 0 suppresses MC/MQ
				if (seg_en - seg_st > 1) {
					int32_t mate_idx = i != seg_en - 1? i + 1 : seg_st; // wrap around if i is the last segment
					mate_qlen = s->seq[mate_idx].l_seq;
				}
				tot_len += t->l_seq;
				if (s->n_hit[i] > 0) { // the query has at least one hit
					int32_t n_sec = 0;
					for (j = 0; j < s->n_hit[i]; ++j) {
						const mb_hit_t *h = &s->hit[i][j];
						if (h->parent == h->id || n_sec < opt->out_n)
							mb_format(km, &out, idx->l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], j, opt, i - seg_st, mate_qlen);
						n_sec += (h->parent != h->id);
					}
				} else if (!(opt->flag & MB_F_NO_UNMAP)) {
					mb_format(km, &out, idx->l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], -1, opt, i - seg_st, mate_qlen);
				}
			}
			fwrite(out.s, 1, out.l, s->p->fp_out);
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_hit[i]; ++j) free(s->hit[i][j].p);
				free(s->hit[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}

		free(out.s);
		free(s->hit); free(s->n_hit); free(s->seq);
		km_destroy(km);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %ld bp in %ld sequences\n", __func__, kom_realtime(), kom_percent_cpu(), (long)tot_len, (long)s->n_seq);
		free(s);
	}
    return 0;
}

static mb_bseq_file_t **mb_open_bseqs(int n, const char **fn)
{
	mb_bseq_file_t **fp;
	int32_t i, j;
	fp = kom_calloc(mb_bseq_file_t*, n);
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mb_bseq_open(fn[i])) == 0) {
			if (kom_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mb_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int32_t mb_map_file(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n, const char **fn, const char *hdr, const char *fn_out)
{
	int32_t i, pl_thread;
	pipeline_t pl;
	if (n < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp_out = fn_out == 0 || strcmp(fn_out, "-") == 0? stdout : fopen(fn_out, "wb");
	if (pl.fp_out == 0) return -1;
	if (hdr) fputs(hdr, pl.fp_out);
	pl.n_fp = n;
	pl.fp = mb_open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.idx = idx;
	pl.mb_size = opt->mb_size;
	pl_thread = opt->n_thread <= 2? opt->n_thread : 3;

	kt_pipeline(pl_thread, worker_pipeline, &pl, 3);

	if (pl.fp_out != stdout) fclose(pl.fp_out);
	for (i = 0; i < n; ++i)
		mb_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

/********************
 * Arbitrary header *
 ********************/

static void mb_insert_hdr(kstring_t *out, const char *str)
{
	out->l = 0;
	if (str[0] == '@' && strlen(str) > 4) {
		char *dup;
		dup = mb_escape(kom_strdup(str));
		kom_sprintf_lite(out, "%s\n", dup);
		free(dup);
	} else if (str[0] != '@') {
		int dret;
		kstring_t tmp = {0,0,0};
		gzFile fp;
		kstream_t *ks;
		fp = gzopen(str, "r");
		if (fp == 0) return;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &tmp, &dret) >= 0)
			if (tmp.s[0] == '@' && tmp.l > 4 && tmp.s[3] == '\t')
				kom_sprintf_lite(out, "%s\n", tmp.s);
		free(tmp.s);
		ks_destroy(ks);
		gzclose(fp);
	}
}

/*******
 * CLI *
 *******/

static ko_longopt_t long_options[] = {
	{ "kalloc",       ko_no_argument,       301 },
	{ "outn",         ko_required_argument, 302 },
	{ "pe-predef",    ko_optional_argument, 303 },
	{ "rescue",       ko_required_argument, 304 },
	{ "eqx",          ko_no_argument,       305 },
	{ "pe",           ko_required_argument, 306 },
	{ "long",         ko_optional_argument, 307 },
	{ "adap",         ko_required_argument, 308 },
	{ "chain-only",   ko_no_argument,       309 },
	{ "meth",         ko_no_argument,       310 },
	{ "hic",          ko_no_argument,       311 },
	{ "xa",           ko_required_argument, 312 },
	{ "dbg-aln-seq",  ko_no_argument,       601 },
	{ "dbg-anchor",   ko_no_argument,       602 },
	{ "dbg-seed",     ko_no_argument,       603 },
	{ "dbg-qname",    ko_no_argument,       604 },
	{ "dbg-aln-pe",   ko_no_argument,       605 },
	{ "dbg-an-pos",   ko_no_argument,       606 }, // anchor position
	{ "version",      ko_no_argument,       901 },
	{ "help",         ko_no_argument,       902 },
	{ 0, 0, 0 }
};

static int usage_map(FILE *fp, const mb_opt_t *opt)
{
	fprintf(fp, "Usage: minibwa map [options] <in.idx> <in.fastq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Common:\n");
	fprintf(fp, "    -f               output PAF (SAM by default)\n");
	fprintf(fp, "    -t INT           number of worker threads [%d]\n", opt->n_thread);
	fprintf(fp, "    -l NUM           treat reads <NUM as short reads in the default adaptive mode [%d]\n", opt->max_sr_len);
	fprintf(fp, "    -R STR           SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
	fprintf(fp, "    -b STR           output a base alignment tag: cs, ds or MD []\n");
	fprintf(fp, "    --hic            map Hi-C reads; equivalent to option -5P\n");
	fprintf(fp, "    --meth           map *directional* bisulfite sequencing reads\n");
	fprintf(fp, "  Mapping:\n");
	fprintf(fp, "    -k INT           min seed length [%d]\n", opt->min_len);
	fprintf(fp, "    -c NUM           max seed occurrences [%d]\n", opt->max_occ);
	fprintf(fp, "    -g NUM           max gap size, controlling extension and chain breaking [%d]\n", opt->max_gap);
	fprintf(fp, "    -w NUM           bandwidth [%d]\n", opt->bw);
	fprintf(fp, "    -W NUM           long bandwidth (for long reads or the adaptive mode) [%d]\n", opt->bw_long);
    fprintf(fp, "    -m INT           min chaining score [%d]\n", opt->min_chain_score);
	fprintf(fp, "    -p FLOAT         min secondary-to-primary score ratio [%g]\n", opt->pri_ratio);
	fprintf(fp, "    -N INT           retain at most INT secondary alignments [%d]\n", opt->best_n);
	fprintf(fp, "    --chain-only     perform chaining only without base alignment\n");
	fprintf(fp, "    -x STR           preset (sr, lr or adap for mixed short/long reads) [adap]\n");
	fprintf(fp, "  Alignment:\n");
	fprintf(fp, "    -A INT           matching score [%d]\n", opt->a);
	fprintf(fp, "    -B INT           mismatching openalty [%d]\n", opt->b);
	fprintf(fp, "    -O INT1[,INT2]   gap open penalty [%d,%d]\n", opt->q, opt->q2);
	fprintf(fp, "    -E INT1[,INT2]   gap extension penalty [%d,%d]\n", opt->e, opt->e2);
	fprintf(fp, "    -s INT           suppress alignment with DP score lower than INT*{-A} [%d]\n", opt->min_dp_max);
	fprintf(fp, "  Paired-end:\n");
	fprintf(fp, "    -P               skip pairing and mate resuce\n");
	fprintf(fp, "    --rescue=INT     mate rescue for up to INT candidates; 0 to skip rescue [%d]\n", opt->max_rescue);
	fprintf(fp, "  Input/Output:\n");
	fprintf(fp, "    -o FILE          output file name [stdout]\n");
	fprintf(fp, "    -u               don't output unmapped reads\n");
	fprintf(fp, "    --outn=NUM       output up to INT secondary alignments [0]\n");
	fprintf(fp, "    --xa=NUM         if <=NUM hits with score >%g%% of the best hit, output them to XA [%d]\n", opt->xa_ratio*100.0, opt->xa_max);
	fprintf(fp, "    -y               copy FASTA/Q comments to output\n");
	fprintf(fp, "    -Y               use soft clipping for supplementary alignments\n");
	fprintf(fp, "    -H STR           if STR starts with @, insert to header; or insert lines in file STR []\n");
	fprintf(fp, "    -5               take the alignment with the smallest query position as primary\n");
	fprintf(fp, "    -K NUM1[,NUM2]   process NUM1-NUM2 bp of query sequences in a batch [100m,1g]\n");
	fprintf(fp, "    --version        print version number\n");
	fprintf(fp, "    --help           print this help message\n");
	return fp == stdout? 0 : 1;
}

static inline void yes_or_no(mb_opt_t *opt, uint64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int main_map(int argc, char *argv[])
{
	const char *opt_str = "x:o:k:c:m:p:A:B:U:b:O:E:t:K:N:PyYR:H:aul:w:W:g:5s:f";
	int32_t c;
	mb_idx_t *idx;
	mb_opt_t mo;
	char *fn_out = 0, *rg_line = 0, *s;
	ketopt_t o = KETOPT_INIT;
	kstring_t hdr_ins = {0,0,0}, hdr = {0,0,0};

	mb_opt_init(&mo);
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mb_opt_preset(&mo, o.arg) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'k') mo.min_len = atoi(o.arg);
		else if (c == 'c') mo.max_occ = kom_parse_num(o.arg, 0);
		else if (c == 'p') mo.pri_ratio = atof(o.arg);
		else if (c == 'm') mo.min_chain_score = atoi(o.arg);
		else if (c == 'N') mo.best_n = atoi(o.arg);
		else if (c == 'A') mo.a = atoi(o.arg);
		else if (c == 'B') mo.b = atoi(o.arg);
		else if (c == 'U') mo.pen_unpair = atoi(o.arg);
		else if (c == 'l') mo.max_sr_len = kom_parse_num(o.arg, 0);
		else if (c == 'g') mo.max_gap = kom_parse_num(o.arg, 0);
		else if (c == 'w') mo.bw = kom_parse_num(o.arg, 0);
		else if (c == 'W') mo.bw_long = kom_parse_num(o.arg, 0);
		else if (c == 'a') mo.flag &= ~MB_F_PAF;
		else if (c == 'f') mo.flag |= MB_F_PAF;
		else if (c == 'u') mo.flag |= MB_F_NO_UNMAP;
		else if (c == 'y') mo.flag |= MB_F_COPY_COMMENT;
		else if (c == 'Y') mo.flag |= MB_F_SUPP_SOFT;
		else if (c == '5') mo.flag |= MB_F_PRIMARY5;
		else if (c == 'P') mo.flag |= MB_F_NO_PAIRING;
		else if (c == 's') mo.min_dp_max = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == 't') mo.n_thread = atoi(o.arg);
		else if (c == 'R') rg_line = o.arg;
		else if (c == 'H') mb_insert_hdr(&hdr_ins, o.arg);
		else if (c == 301) { // --kalloc
			yes_or_no(&mo, MB_F_NO_KALLOC, o.longidx, o.arg, 0);
		} else if (c == 302) { // --outn
			mo.out_n = kom_parse_num(o.arg, 0);
		} else if (c == 303) { // --pe-predef
			mo.flag |= MB_F_PE_PREDEF;
		} else if (c == 304) { // --rescue
			mo.max_rescue = atoi(o.arg);
		} else if (c == 305) { // --eqx
			mo.flag |= MB_F_EQX;
		} else if (c == 306) { // --pe
			yes_or_no(&mo, MB_F_PE, o.longidx, o.arg, 1);
		} else if (c == 307) { // --long
			if (o.arg == 0) mo.flag |= MB_F_LONG;
			else yes_or_no(&mo, MB_F_LONG, o.longidx, o.arg, 1);
		} else if (c == 308) { // --adap
			yes_or_no(&mo, MB_F_ADAP, o.longidx, o.arg, 1);
		} else if (c == 309) { // --chain-only
			mo.flag |= MB_F_NO_ALN;
		} else if (c == 310) { // --meth
			mo.flag |= MB_F_METH;
		} else if (c == 311) { // --hic
			mo.flag |= MB_F_PRIMARY5 | MB_F_NO_PAIRING;
		} else if (c == 312) { // --xa
			mo.xa_max = kom_parse_num(o.arg, 0);
		} else if (c == 601) { // --dbg-aln-seq
			kom_dbg_flag |= MB_DBG_ALN_SEQ;
		} else if (c == 602) { // --dbg-anchor
			kom_dbg_flag |= MB_DBG_ANCHOR;
		} else if (c == 603) { // --dbg-seed
			kom_dbg_flag |= MB_DBG_SEED;
		} else if (c == 604) { // --dbg-qname
			kom_dbg_flag |= MB_DBG_QNAME;
		} else if (c == 605) { // --dbg-aln-pe
			kom_dbg_flag |= MB_DBG_ALN_PE;
		} else if (c == 606) { // --dbg-an-pos
			kom_dbg_flag |= MB_DBG_AN_POS;
		} else if (c == 'K') {
			mo.mb_size = mo.max_mb_size = kom_parse_num(o.arg, &s);
			if (*s == ',') mo.max_mb_size = kom_parse_num(s + 1, &s);
			if (mo.max_mb_size < mo.mb_size)
				mo.max_mb_size = mo.mb_size;
		} else if (c == 'O') {
			mo.q = mo.q2 = strtol(o.arg, &s, 10);
			if (*s == ',') mo.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			mo.e = mo.e2 = strtol(o.arg, &s, 10);
			if (*s == ',') mo.e2 = strtol(s + 1, &s, 10);
		} else if (c == 'b') {
			mo.flag &= ~(MB_F_WRITE_CS|MB_F_WRITE_DS|MB_F_WRITE_MD);
			if (strcmp(o.arg, "cs") == 0) mo.flag |= MB_F_WRITE_CS;
			else if (strcmp(o.arg, "ds") == 0) mo.flag |= MB_F_WRITE_DS;
			else if (strcmp(o.arg, "MD") == 0 || strcmp(o.arg, "md") == 0) mo.flag |= MB_F_WRITE_MD;
			else if (kom_verbose >= 2) {
				mo.flag |= MB_F_WRITE_CS;
				fprintf(stderr, "[WARNING]\033[1;31m -b only takes 'cs', 'ds' or 'MD'. Invalid values are assumed to be 'cs'.\033[0m\n");
			}
		} else if (c == 901) { // --version
			puts(MB_VERSION);
			exit(0);
		} else if (c == 902) { // --help
			return usage_map(stdout, &mo);
		}
	}
	if (argc - o.ind < 2)
		return usage_map(stderr, &mo);

	idx = mb_idx_load(argv[o.ind], !!(mo.flag & MB_F_METH));
	kom_assert(idx, "failed to load the index.");
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] index loaded\n", __func__, kom_realtime(), kom_percent_cpu());

	if (!(mo.flag & MB_F_PAF)) {
		int ret;
		ret = mb_fmt_sam_hdr(&hdr, idx->l2b, rg_line, MB_VERSION, argc, argv);
		if (ret < 0) return 1; // TODO: free idx and out.s
		if (hdr_ins.l > 0) kom_sprintf_lite(&hdr, "%s", hdr_ins.s);
	}
	if (hdr_ins.s) free(hdr_ins.s);

	mb_map_file(&mo, idx, argc - (o.ind + 1), (const char**)&argv[o.ind+1], hdr.s, fn_out);
	free(hdr.s);
	mb_idx_destroy(idx);
	return 0;
}

static int usage_mem(FILE *fp, const mb_opt_t *opt)
{
	fprintf(fp, "Usage: minibwa mem [options] <in.idx> <in.fastq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Algorithm:\n");
	fprintf(fp, "    -t INT         number of threads [%d]\n", opt->n_thread);
	fprintf(fp, "    -k INT         minimum seed length [%d]\n", opt->min_len);
	fprintf(fp, "    -w NUM         bandwidth [%d]\n", opt->bw);
	fprintf(fp, "    *r FLOAT       inside seed ratio []\n");
	fprintf(fp, "    *y INT         seed occurrence for the 3rd round seeding []\n");
	fprintf(fp, "    -c NUM         max seed occurrences [%d]\n", opt->max_occ);
	fprintf(fp, "    -D FLOAT       mask level [%g]\n", opt->mask_level);
	fprintf(fp, "    -W INT         min chaining score [%d]\n", opt->min_chain_score);
	fprintf(fp, "    -m INT         mate rescue for up to INT candidates; 0 to skip rescue [%d]\n", opt->max_rescue);
	fprintf(fp, "    -S             skip mate rescue\n");
	fprintf(fp, "    -P             skip pairing\n");
	fprintf(fp, "  Scoring:\n");
	fprintf(fp, "    *A INT         matching score []\n");
	fprintf(fp, "    *B INT         mismatching penalty []\n");
	fprintf(fp, "    *O INT         gap open penalty []\n");
	fprintf(fp, "    *E INT         gap extension penalty []\n");
	fprintf(fp, "    *L INT         end bonus []\n");
	fprintf(fp, "    *U INT         penalty for an unpaired read pair []\n");
	fprintf(fp, "    *d INT         off-diagonal X-dropoff []\n");
	fprintf(fp, "    *x STR         preset []\n");
	fprintf(fp, "  Input/output:\n");
	fprintf(fp, "    *p             smart pairing\n");
	fprintf(fp, "    -R STR         SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
	fprintf(fp, "    -H STR         if STR starts with @, insert to header; or insert lines in file STR []\n");
	fprintf(fp, "    -o FILE        output file name [stdout]\n");
	fprintf(fp, "    *j             treat ALT contigs as part of the primary assembly\n");
	fprintf(fp, "    -5             take the alignment with the smallest query position as primary\n");
	fprintf(fp, "    *q             don't modify mapQ of supplementary alignments\n");
	fprintf(fp, "    *K NUM         batch size []\n");
	fprintf(fp, "    -v INT         verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", kom_verbose);
	fprintf(fp, "    -T INT         suppress alignment with DP score lower than INT*{-A} [%d]\n", opt->min_dp_max);
	fprintf(fp, "    -h INT         max secondary alignments at the XA tag [%d]\n", opt->xa_max);
	fprintf(fp, "    -z FLOAT       ignore in XA is score < FLOAT*bestScore [%g]\n", opt->xa_ratio);
	fprintf(fp, "    -a             output all alignments\n");
	fprintf(fp, "    -C             copy FASTA/Q comments to output\n");
	fprintf(fp, "    *V             output the reference FASTA header in the XR tag\n");
	fprintf(fp, "    -Y             use soft clipping for supplementary alignments\n");
	fprintf(fp, "    *M             mark shorter split hits as secondary\n");
	fprintf(fp, "    -I FLOAT[,FLOAT[,INT[,INT]]]\n");
	fprintf(fp, "                   mean, stddev, max and min of isize distribution [inferred]\n");
	fprintf(fp, "    *u             output XB instead of XA\n");
	fprintf(fp, "Notes:\n");
	fprintf(fp, "  - \"minibwa mem\" aims to match the \"bwa mem\" command-line interface\n");
	fprintf(fp, "  - '*' options are ignored as they are missing or incompatible with minibwa\n");
	fprintf(fp, "  - minibwa and bwa-mem may output different alignments\n");
	return fp == stdout? 0 : 1;
}

int main_mem(int argc, char *argv[])
{
	int32_t c, ret;
	ketopt_t o = KETOPT_INIT;
	mb_opt_t mo;
	char *fn_out = 0, *rg_line = 0;
	kstring_t hdr_ins = {0,0,0}, hdr = {0,0,0};
	mb_idx_t *idx;

	mb_opt_init(&mo);
	mo.flag |= MB_F_WRITE_MD; // bwa-mem always writes MD
	while ((c = ketopt(&o, argc, argv, 1, "t:k:w:d:r:y:c:D:W:m:SPA:B:O:E:L:U:x:pR:H:o:j5qK:v:T:h:z:aCVYMuI:", 0)) >= 0) {
		// algorithm
		if (c == 't') mo.n_thread = atoi(o.arg);
		else if (c == 'k') mo.min_len = atoi(o.arg);
		else if (c == 'w') mo.bw = kom_parse_num(o.arg, 0);
		else if (c == 'r' || c == 'y') {}
		else if (c == 'W') mo.min_chain_score = atoi(o.arg);
		else if (c == 'c') mo.max_occ = kom_parse_num(o.arg, 0);
		else if (c == 'D') mo.mask_level = atof(o.arg);
		else if (c == 'm') mo.max_rescue = atoi(o.arg);
		else if (c == 'S') mo.max_rescue = 0;
		else if (c == 'P') mo.flag |= MB_F_NO_PAIRING;
		// scoring
		else if (c == 'A' || c == 'B' || c == 'O' || c == 'E' || c == 'L' || c == 'U' || c == 'd' || c == 'x') {}
		// input/output
		else if (c == 'p' || c == 'j' || c == 'q' || c == 'K' || c == 'V' || c == 'M' || c == 'u') {}
		else if (c == 'a') mo.out_n = 1000000;
		else if (c == 'h') mo.xa_max = atoi(o.arg);
		else if (c == 'z') mo.xa_ratio = atof(o.arg);
		else if (c == 'R') rg_line = o.arg;
		else if (c == 'H') mb_insert_hdr(&hdr_ins, o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == '5') mo.flag |= MB_F_PRIMARY5;
		else if (c == 'v') kom_verbose = atoi(o.arg);
		else if (c == 'T') mo.min_dp_max = atoi(o.arg);
		else if (c == 'C') mo.flag |= MB_F_COPY_COMMENT;
		else if (c == 'Y') mo.flag |= MB_F_SUPP_SOFT;
		else if (c == 'I') {
			char *q;
			mo.pe_avg = strtod(o.arg, &q);
			mo.pe_std = (int32_t)(.499 + (*q == ','? strtod(q + 1, &q) : mo.pe_avg * 0.1));
			mo.pe_hi = *q == ','? strtol(q + 1, &q, 10) : mo.pe_avg + mo.pe_std * 4;
			mo.pe_lo = *q == ','? strtol(q + 1, &q, 10) : mo.pe_avg - mo.pe_std * 4;
			if (mo.pe_lo < 1) mo.pe_lo = 1;
		}
	}
	if (argc - o.ind < 2) return usage_mem(stderr, &mo);

	idx = mb_idx_load(argv[o.ind], !!(mo.flag & MB_F_METH));
	kom_assert(idx, "failed to load the index.");
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] index loaded\n", __func__, kom_realtime(), kom_percent_cpu());

	ret = mb_fmt_sam_hdr(&hdr, idx->l2b, rg_line, MB_VERSION, argc, argv);
	if (ret < 0) return 1; // TODO: free idx and out.s
	if (hdr_ins.l > 0) kom_sprintf_lite(&hdr, "%s", hdr_ins.s);
	if (hdr_ins.s) free(hdr_ins.s);
	mb_map_file(&mo, idx, argc - (o.ind + 1), (const char**)&argv[o.ind+1], hdr.s, fn_out);
	free(hdr.s);
	mb_idx_destroy(idx);
	return 0;
}
