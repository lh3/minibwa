#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#include "kommon.h"
#include "mbpriv.h"
#include "bseq.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"

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
	int32_t n, j, k, l, p, *len;
	uint8_t **seq;
	mb_sai_v *sai;

	km = mb_tbuf_km(b);
	for (k = 0, n = 0; k < s->sb_cnt[i]; ++k)
		n += s->seg_cnt[s->sb_off[i] + k];
	len = Kmalloc(km, int32_t, n);
	seq = Kmalloc(km, uint8_t*, n);
	sai = Kcalloc(km, mb_sai_v, n);
	for (k = p = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j) {
			const mb_bseq1_t *t = &s->seq[off + j];
			len[p] = t->l_seq;
			seq[p] = Kcalloc(km, uint8_t, t->l_seq);
			for (l = 0; l < t->l_seq; ++l)
				seq[p][l] = kom_nt4_table[(uint8_t)t->seq[l]];
			++p;
		}
	}
	assert(p == n);
	mb_seed_intv_batch(km, idx->bwt, n, len, seq, opt->min_len, opt->max_sub_occ, sai);
	for (k = p = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j) {
			const mb_bseq1_t *t = &s->seq[off + j];
			mb_opt_t opt_adap;
			mb_opt_adap(opt, len[p], &opt_adap);
			s->hit[off+j] = mb_map_sai(&opt_adap, idx, len[p], seq[p], &sai[p], &s->n_hit[off+j], b, t->name);
			kfree(km, seq[p]);
			++p;
		}
	}
	kfree(km, sai);
	kfree(km, seq);
	kfree(km, len);
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
		int with_qual = !!(opt->flag & MB_F_SAM);
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
				tot_len += t->l_seq;
				if (s->n_hit[i] > 0) { // the query has at least one hit
					int32_t n_sec = 0;
					for (j = 0; j < s->n_hit[i]; ++j) {
						const mb_hit_t *h = &s->hit[i][j];
						if (h->parent == h->id || n_sec < opt->out_n)
							mb_format(km, &out, idx->l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], j, opt->flag, i - seg_st);
						n_sec += (h->parent != h->id);
					}
				} else if (!(opt->flag & MB_F_NO_UNMAP)) { // TODO: output unmapped reads
					mb_format(km, &out, idx->l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], -1, opt->flag, i - seg_st);
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

int32_t mb_map_file(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n, const char **fn, const char *fn_out)
{
	int32_t i, pl_thread;
	pipeline_t pl;
	if (n < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp_out = fn_out == 0 || strcmp(fn_out, "-") == 0? stdout : fopen(fn_out, "wb");
	if (pl.fp_out == 0) return -1;
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

static int usage(FILE *fp, const mb_opt_t *opt)
{
	fprintf(fp, "Usage: minibwa map [options] <in.idx> <in.fastq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Common:\n");
	fprintf(fp, "    -a               output SAM (PAF by default)\n");
	fprintf(fp, "    -t INT           number of worker threads [%d]\n", opt->n_thread);
	fprintf(fp, "    -l NUM           treat reads <NUM as short reads in the default adaptive mode [%d]\n", opt->max_sr_len);
	fprintf(fp, "    -R STR           SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
	fprintf(fp, "    -b STR           output a base alignment tag: cs, ds or MD []\n");
	fprintf(fp, "    -5               take the alignment with the smallest query position as primary\n");
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
	fprintf(fp, "    --outn=INT       output up to INT secondary alignments [0]\n");
	fprintf(fp, "    -y               copy FASTA/Q comments to output\n");
	fprintf(fp, "    -Y               use soft clipping for supplementary alignments\n");
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
	const char *opt_str = "x:o:k:c:m:p:A:B:b:O:E:t:K:N:PyYR:aul:w:W:g:5s:";
	int32_t c;
	mb_idx_t *idx;
	mb_opt_t mo;
	char *fn_out = 0, *rg_line = 0, *s;
	ketopt_t o = KETOPT_INIT;

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
		else if (c == 'l') mo.max_sr_len = kom_parse_num(o.arg, 0);
		else if (c == 'g') mo.max_gap = kom_parse_num(o.arg, 0);
		else if (c == 'w') mo.bw = kom_parse_num(o.arg, 0);
		else if (c == 'W') mo.bw_long = kom_parse_num(o.arg, 0);
		else if (c == 'a') mo.flag |= MB_F_SAM;
		else if (c == 'u') mo.flag |= MB_F_NO_UNMAP;
		else if (c == 'y') mo.flag |= MB_F_COPY_COMMENT;
		else if (c == 'Y') mo.flag |= MB_F_SUPP_SOFT;
		else if (c == '5') mo.flag |= MB_F_PRIMARY5;
		else if (c == 'P') mo.flag |= MB_F_NO_PAIRING;
		else if (c == 's') mo.min_dp_max = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == 't') mo.n_thread = atoi(o.arg);
		else if (c == 'R') rg_line = o.arg;
		else if (c == 301) { // --kalloc
			yes_or_no(&mo, MB_F_NO_KALLOC, o.longidx, o.arg, 0);
		} else if (c == 302) { // --outn
			mo.out_n = atoi(o.arg);
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
			return usage(stdout, &mo);
		}
	}
	if (argc - o.ind < 2)
		return usage(stderr, &mo);

	idx = mb_idx_load(argv[o.ind]);
	kom_assert(idx, "failed to load the index.");
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] index loaded\n", __func__, kom_realtime(), kom_percent_cpu());

	if (mo.flag & MB_F_SAM) {
		int ret;
		kstring_t out = {0,0,0};
		ret = mb_fmt_sam_hdr(&out, idx->l2b, rg_line, MB_VERSION, argc, argv);
		if (ret < 0) return 1; // TODO: free idx and out.s
		fwrite(out.s, 1, out.l, stdout);
		free(out.s);
	}

	mb_map_file(&mo, idx, argc - (o.ind + 1), (const char**)&argv[o.ind+1], fn_out);
	mb_idx_destroy(idx);
	return 0;
}
