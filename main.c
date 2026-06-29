#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "kommon.h"
#include "mbpriv.h"
#include "ketopt.h"

int main_index(int argc, char *argv[]);
int main_map(int argc, char *argv[]);
int main_mem(int argc, char *argv[]);

int main_fa2bit(int argc, char *argv[]);
int main_raw2bwt(int argc, char *argv[]);
int main_genraw(int argc, char *argv[]);
int main_genbwt(int argc, char *argv[]);
int main_gensa(int argc, char *argv[]);

int main_getref(int argc, char *argv[]);
int main_fastmap(int argc, char *argv[]);
int main_bench(int argc, char *argv[]);

static int usage(FILE *fp, int is_long)
{
	fprintf(fp, "Usage: minibwt <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	if (is_long) {
		fprintf(fp, "  General:\n");
		fprintf(fp, "    index      index reference FASTA\n");
		fprintf(fp, "    map        read alignment\n");
		fprintf(fp, "    mem        legacy bwa-mem CLI (not recommended)\n");
		fprintf(fp, "    version    print the version number\n");
		fprintf(fp, "  Separate indexing routines:\n");
		fprintf(fp, "    fa2bit     convert FASTA to the long-2bit format\n");
		fprintf(fp, "    genraw     generate BWT from pac with the BWT-SW algorithm\n");
		fprintf(fp, "    raw2bwt    recode bwtgen raw BWT\n");
		fprintf(fp, "    gensa      generate sampled SA from BWT\n");
		fprintf(fp, "    genbwt     generate BWT+SSA from long-2bit with libsais\n");
		fprintf(fp, "  Debugging:\n");
		fprintf(fp, "    bench      performance evaluation\n");
		fprintf(fp, "    fastmap    test seeding strategies\n");
		fprintf(fp, "    getref     get the reference genome from .l2b\n");
		fprintf(fp, "  Help:\n");
		fprintf(fp, "    --help     print this help message\n");
	} else {
		fprintf(fp, "  index      index reference FASTA\n");
		fprintf(fp, "  map        read alignment\n");
		fprintf(fp, "  mem        legacy bwa-mem CLI (not recommended)\n");
		fprintf(fp, "  version    print the version number\n");
	}
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc == 1) return usage(stdout, 0);
	else if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "map") == 0) ret = main_map(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2bit") == 0) ret = main_fa2bit(argc-1, argv+1);
	else if (strcmp(argv[1], "genraw") == 0) ret = main_genraw(argc-1, argv+1);
	else if (strcmp(argv[1], "raw2bwt") == 0) ret = main_raw2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "genbwt") == 0) ret = main_genbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "gensa") == 0) ret = main_gensa(argc-1, argv+1);
	else if (strcmp(argv[1], "getref") == 0) ret = main_getref(argc-1, argv+1);
	else if (strcmp(argv[1], "bench") == 0) ret = main_bench(argc-1, argv+1);
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else if (strcmp(argv[1], "--help") == 0) return usage(stdout, 1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", MB_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (kom_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MB_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, kom_realtime(), kom_cputime(), kom_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

static int usage_getref(FILE *fp)
{
	fprintf(fp, "Usage: minibwa getref <ref.l2b>\n");
	return fp == stdout? 0 : 1;
}

int main_getref(int argc, char *argv[])
{
	uint64_t i, j;
	l2b_t *l2b;
	if (argc == 1) return usage_getref(stderr);
	l2b = l2b_load(argv[1]);
	assert(l2b);
	for (i = 0; i < l2b->n_ctg; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[i];
		uint8_t *seq;
		printf(">%s\n", ctg->name);
		seq = kom_malloc(uint8_t, ctg->len);
		l2b_getseq(l2b, i, 0, ctg->len, seq);
		for (j = 0; j < ctg->len; ++j) seq[j] = "ACGTN"[seq[j]];
		fwrite(seq, 1, j, stdout);
		fputc('\n', stdout);
		free(seq);
	}
	return 0;
}

typedef enum { MB_BENCH_2A, MB_BENCH_SA, MB_BENCH_MSA } mb_bench_type_t;

static int usage_bench(FILE *fp, int intv)
{
	fprintf(fp, "Usage: minibwa bench [options] <in.mbw>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -b STR         type: 2a, sa or msa [2a]\n");
	fprintf(fp, "  -n NUM         number of data points [1m]\n");
	fprintf(fp, "  -v INT         interval size for msa [%d]\n", intv);
	fprintf(fp, "  -p             print results for each data point\n");
	fprintf(fp, "  -1             use unbatched sa for msa\n");
	fprintf(fp, "  --help         print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_bench(int argc, char *argv[])
{
	mb_bench_type_t type = MB_BENCH_2A;
	uint64_t x = 11, cs = 1;
	int64_t i, n = 1000000;
	int c, print_val = 0, use_single = 0, intv = 20;
	mb_bwt_t *bwt;
	ketopt_t o = KETOPT_INIT;
	double t;
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};

	while ((c = ketopt(&o, argc, argv, 1, "pn:b:v:1", long_opts)) >= 0) {
		if (c == 'n') n = kom_parse_num(o.arg, 0);
		else if (c == 'p') print_val = 1;
		else if (c == '1') use_single = 1;
		else if (c == 'v') intv = atoi(o.arg);
		else if (c == 'b') {
			if (strcmp(o.arg, "2a") == 0) type = MB_BENCH_2A;
			else if (strcmp(o.arg, "sa") == 0) type = MB_BENCH_SA;
			else if (strcmp(o.arg, "msa") == 0) type = MB_BENCH_MSA;
			else kom_assert(0, "unknown type");
		}
		else if (c == 901) return usage_bench(stdout, intv);
	}
	if (argc - o.ind < 1) return usage_bench(stderr, intv);

	bwt = mb_bwt_load(argv[o.ind]);
	t = kom_cputime();
	if (type == MB_BENCH_2A) {
		for (i = 0; i < n; ++i) {
			uint64_t k = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t l = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t cntk[4], cntl[4];
			mb_bwt_rank2a(bwt, k, l, cntk, cntl);
			cs = cs * cntk[1] + cntl[0];
			if (print_val) printf("%ld\n", (long)cntk[1]);
		}
	} else if (type == MB_BENCH_SA) {
		for (i = 0; i < n; ++i) {
			uint64_t s, k = kom_splitmix64(&x) % bwt->seq_len;
			s = mb_bwt_sa(bwt, k);
			cs = cs * 0xbf58476d1ce4e5b9ULL ^ s;
			if (print_val) printf("%ld\n", (long)s);
		}
	} else if (type == MB_BENCH_MSA) {
		for (i = 0; i < n; ++i) {
			uint64_t j, xor = 0, k = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t l = k + intv < bwt->seq_len? k + intv : bwt->seq_len;
			if (use_single) {
				for (j = k; j < l; ++j) {
					uint64_t s = mb_bwt_sa(bwt, j);
					xor ^= s;
				}
			} else {
				uint64_t sa[intv], n_sa = l - k;
				for (j = 0; j < n_sa; ++j) sa[j] = k + j;
				mb_bwt_sa_batch(0, bwt, l - k, sa);

				for (j = 0; j < n_sa; ++j) xor ^= sa[j];
			}
			cs = cs * 0xbf58476d1ce4e5b9ULL ^ xor;
			if (print_val) printf("%ld\n", (long)xor);
		}
	}
	fprintf(stderr, "checksum = %lx\n", (unsigned long)cs);
	fprintf(stderr, "t = %.3f\n", kom_cputime() - t);
	mb_bwt_destroy(bwt);
	return 0;
}
