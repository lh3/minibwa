#include <stdio.h>
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
	const mb_mopt_t *opt;
	mb_bseq_file_t **fp;
	const mb_idx_t *idx;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int32_t n_seq, n_frag;
	int32_t *n_hit, *seg_off, *n_seg;
	mb_bseq1_t *seq;
	mb_hit_t **hit;
	mb_tbuf_t **tbuf;
} step_t;

static void worker_for(void *_data, long i, int tid)
{
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = !!(p->opt->flag & MB_F_SAM);
		int with_comment = (!!(p->opt->flag & MB_F_SAM) && !!(p->opt->flag & MB_F_COPY_COMMENT));
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MB_F_FRAG_MODE));
        step_t *s;
        s = kom_calloc(step_t, 1);
		if (p->n_fp > 1) s->seq = mb_bseq_read_frag(p->n_fp, p->fp, p->mb_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mb_bseq_read(p->fp[0], p->mb_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].id = p->n_seq++;
			s->tbuf = kom_calloc(mb_tbuf_t*, p->opt->n_thread);
			for (i = 0; i < p->opt->n_thread; ++i)
				s->tbuf[i] = mb_tbuf_init();
			s->n_hit = kom_calloc(int32_t, 3 * s->n_seq);
			s->n_seg = s->n_hit + s->n_seq;
			s->seg_off = s->n_seg + s->n_seq;
			s->hit = kom_calloc(mb_hit_t*, s->n_seq);
			for (i = 1, j = 0; i <= s->n_seq; ++i) { // set n_seg[] and seg_off[]
				if (i == s->n_seq || !frag_mode || !mb_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->opt->n_thread, worker_for, in, ((step_t*)in)->n_frag);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mb_idx_t *idx = p->idx;
		for (i = 0; i < p->opt->n_thread; ++i)
			mb_tbuf_destroy(s->tbuf[i]);
		free(s->tbuf);
		if (!(kom_dbg_flag & MB_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int32_t seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mb_bseq1_t *t = &s->seq[i];
				if (s->n_hit[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_hit[i]; ++j) {
						const mb_hit_t *h = &s->hit[i][j];
					}
				} else if (p->opt->flag & MB_F_WRITE_UNMAP) {
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_hit[i]; ++j) free(s->hit[i][j].p);
				free(s->hit[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->hit); free(s->n_hit); free(s->seq);
		km_destroy(km);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %ld sequences\n", __func__, kom_realtime(), kom_percent_cpu(), (long)s->n_seq);
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

int32_t mb_map_file(const mb_idx_t *idx, int32_t n, const char **fn, const mb_mopt_t *opt)
{
	int32_t i;
	pipeline_t pl;
	if (n < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n;
	pl.fp = mb_open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.idx = idx;
	pl.mb_size = opt->mb_size;

	kt_pipeline(3, worker_pipeline, &pl, 3);

	free(pl.str.s);
	for (i = 0; i < n; ++i)
		mb_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

static ko_longopt_t long_options[] = {
	{ "version",      ko_no_argument,     901 },
	{ 0, 0, 0 }
};

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minibwa map [options] <in.idx> <in.fastq>\n");
	return fp == stdout? 0 : 1;
}

int main_map(int argc, char *argv[])
{
	int32_t c;
	mb_idx_t *idx;
	mb_mopt_t mo;
	ketopt_t o = KETOPT_INIT;

	kom_realtime(); // reset the timer
	mb_mopt_init(&mo);
	while ((c = ketopt(&o, argc, argv, 1, "k:", long_options)) >= 0) {
		if (c == 'k') mo.min_k = atoi(o.arg);
		else if (c == 901) { // --version
			puts(MB_VERSION);
			exit(0);
		}
	}
	if (argc - o.ind < 2)
		return usage(stderr);

	idx = mb_idx_load(argv[o.ind]);
	kom_assert(idx, "failed to load the index.");
	mb_map_file(idx, argc - o.ind - 1, &argv[o.ind+1], &mo);
	mb_idx_destroy(idx);
	return 0;
}
