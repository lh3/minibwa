#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "sys.h"
#include "utils.h"
#include "ketopt.h"
#include "l2bit.h"
#include "hl.h"

#define MB_VERSION "0.0"

int main_fa2bit(int argc, char *argv[]);
int main_raw2bwt(int argc, char *argv[]);
int main_test(int argc, char *argv[]);
int main_genraw(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minibwt <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  fa2bit     convert FASTA to the long-2-bit format\n");
	fprintf(fp, "  genraw     build BWT from pac with the BWT-SW algorithm\n");
	fprintf(fp, "  raw2bwt    recode bwtgen raw BWT\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	mb_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "fa2bit") == 0) ret = main_fa2bit(argc-1, argv+1);
	else if (strcmp(argv[1], "genraw") == 0) ret = main_genraw(argc-1, argv+1);
	else if (strcmp(argv[1], "raw2bwt") == 0) ret = main_raw2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "test") == 0) ret = main_test(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", MB_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (mb_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MB_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, mb_realtime(), mb_cputime(), mb_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

int main_fa2bit(int argc, char *argv[])
{
	l2b_t *l2b;
	int out_pac = 0, both_strand = 0;
	uint64_t seed = 11;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "s:p2", 0)) >= 0) {
		if (c == 's') seed = atol(o.arg);
		else if (c == 'p') out_pac = 1;
		else if (c == '2') both_strand = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa fa2bit [options] <in.fa> <out.l2b>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s INT    random seed [%lu]\n", (unsigned long)seed);
		fprintf(stderr, "  -p        output the BWA pac format\n");
		fprintf(stderr, "  -2        output both strands (effective with -p)\n");
		return 1;
	}
	l2b = l2b_import(argv[o.ind], seed);
	if (out_pac)
		l2b_save_pac(argv[o.ind+1], l2b, both_strand);
	else
		l2b_save(argv[o.ind+1], l2b);
	l2b_destroy(l2b);
	return 0;
}

int main_genraw(int argc, char *argv[])
{
	extern void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);
	ketopt_t o = KETOPT_INIT;
	int c, block_size = 10000000;
	while ((c = ketopt(&o, argc, argv, 1, "b:", 0)) >= 0) {
		if (c == 'b') block_size = hl_parse_num(o.arg, 0);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa genraw [options] <in.pac> <out.raw-bwt>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b NUM      block size [10m]\n");
		return 1;
	}
	mb_bwtgen(argv[o.ind], argv[o.ind+1], block_size);
	return 0;
}

int main_raw2bwt(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	if (argc < 3) {
		printf("Usage: minibwa raw2bwt <raw.bwt> <recode.bwt>\n");
		return 1;
	}
	bwt = mb_bwt_load_raw(argv[1]);
	mb_bwt_save(argv[2], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

int main_test(int argc, char *argv[])
{
	uint64_t i, n = 10000, x = 11;
	mb_bwt_t *bwt;
	if (argc == 1) {
		fprintf(stderr, "Usage: minibwa test <in.bwt>\n");
		return 1;
	}
	if (argc >= 3) n = atol(argv[2]);

	bwt = mb_bwt_load(argv[1]);
	if (n > 0) {
		double t = mb_cputime();
		for (i = 0; i < n; ++i) {
			uint64_t k = hl_splitmix64(&x) % bwt->seq_len;
			#if 1
			uint64_t cnt[4];
			mb_bwt_rank1a(bwt, k, cnt);
			printf("%lld\n", cnt[1]);
			#else
			printf("%lld\n", mb_bwt_rank11(bwt, k, 1));
			#endif
		}
		fprintf(stderr, "t = %.3f\n", mb_cputime() - t);
	} else {
		uint64_t k = 10, cnt[4];
		int c;
		mb_bwt_rank1a(bwt, k, cnt);
		for (c = 0; c < 4; ++c)
			printf("%lld\n", cnt[c]);
	}
	mb_bwt_destroy(bwt);
	return 0;
}
