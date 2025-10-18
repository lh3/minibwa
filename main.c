#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "ketopt.h"
#include "kommon.h"
#include "mbpriv.h"

#define MB_VERSION "0.0"

int main_index(int argc, char *argv[]);
int main_fa2bit(int argc, char *argv[]);
int main_raw2bwt(int argc, char *argv[]);
int main_test(int argc, char *argv[]);
int main_genraw(int argc, char *argv[]);
int main_genbwt(int argc, char *argv[]);
int main_gensa(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minibwt <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  General:\n");
	fprintf(fp, "    index      index reference FASTA\n");
	fprintf(fp, "    version    print the version number\n");
	fprintf(fp, "  Separate indexing routines:\n");
	fprintf(fp, "    fa2bit     convert FASTA to the long-2bit format\n");
	fprintf(fp, "    genraw     generate BWT from pac with the BWT-SW algorithm\n");
	fprintf(fp, "    raw2bwt    recode bwtgen raw BWT\n");
	fprintf(fp, "    gensa      generate sampled SA from BWT\n");
	fprintf(fp, "    genbwt     generate BWT+SSA from long-2bit with libsais\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2bit") == 0) ret = main_fa2bit(argc-1, argv+1);
	else if (strcmp(argv[1], "genraw") == 0) ret = main_genraw(argc-1, argv+1);
	else if (strcmp(argv[1], "raw2bwt") == 0) ret = main_raw2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "genbwt") == 0) ret = main_genbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "gensa") == 0) ret = main_gensa(argc-1, argv+1);
	else if (strcmp(argv[1], "test") == 0) ret = main_test(argc-1, argv+1);
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

int main_test(int argc, char *argv[])
{
	uint64_t x = 11;
	int64_t i, n = 10000;
	mb_bwt_t *bwt;
	if (argc == 1) {
		fprintf(stderr, "Usage: minibwa test <in.bwt>\n");
		return 1;
	}
	if (argc >= 3) n = atol(argv[2]);

	bwt = mb_bwt_load(argv[1]);
	if (n > 0) {
		double t = kom_cputime();
		for (i = 0; i < n; ++i) {
			uint64_t k = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t l = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t cntk[4], cntl[4];
			mb_bwt_rank2a(bwt, k, l, cntk, cntl);
			//printf("%lld\n", cntk[1]);
		}
		fprintf(stderr, "t = %.3f\n", kom_cputime() - t);
	} else {
		uint64_t k = 10, cnt[4];
		int c;
		mb_bwt_rank1a(bwt, k, cnt);
		for (c = 0; c < 4; ++c)
			printf("%lld\n", cnt[c]);
		//uint64_t j; for (j = 0; j < bwt->n_sa; ++j) printf("%lld\n", bwt->sa[j]);
	}
	mb_bwt_destroy(bwt);
	return 0;
}
