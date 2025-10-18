#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kommon.h"
#include "mbpriv.h"
#include "ketopt.h"

#define MB_VERSION "0.0"

int main_index(int argc, char *argv[]);
int main_seed(int argc, char *argv[]);
int main_bench(int argc, char *argv[]);

int main_fa2bit(int argc, char *argv[]);
int main_raw2bwt(int argc, char *argv[]);
int main_genraw(int argc, char *argv[]);
int main_genbwt(int argc, char *argv[]);
int main_gensa(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minibwt <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  General:\n");
	fprintf(fp, "    index      index reference FASTA\n");
	fprintf(fp, "    bench      performance evaluation\n");
	fprintf(fp, "    seed       test seeding strategies\n");
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
	else if (strcmp(argv[1], "bench") == 0) ret = main_bench(argc-1, argv+1);
	else if (strcmp(argv[1], "seed") == 0) ret = main_seed(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2bit") == 0) ret = main_fa2bit(argc-1, argv+1);
	else if (strcmp(argv[1], "genraw") == 0) ret = main_genraw(argc-1, argv+1);
	else if (strcmp(argv[1], "raw2bwt") == 0) ret = main_raw2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "genbwt") == 0) ret = main_genbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "gensa") == 0) ret = main_gensa(argc-1, argv+1);
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

typedef enum { MB_BENCH_2A, MB_BENCH_SA, MB_BENCH_MSA } mb_bench_type_t;

int main_bench(int argc, char *argv[])
{
	mb_bench_type_t type = MB_BENCH_2A;
	uint64_t x = 11, cs = 1;
	int64_t i, n = 1000000;
	int c, print_val = 0, intv = 20;
	mb_bwt_t *bwt;
	ketopt_t o = KETOPT_INIT;
	double t;

	while ((c = ketopt(&o, argc, argv, 1, "pn:b:v:", 0)) >= 0) {
		if (c == 'n') n = kom_parse_num(o.arg, 0);
		else if (c == 'p') print_val = 1;
		else if (c == 'v') intv = atoi(o.arg);
		else if (c == 'b') {
			if (strcmp(o.arg, "2a") == 0) type = MB_BENCH_2A;
			else if (strcmp(o.arg, "sa") == 0) type = MB_BENCH_SA;
			else if (strcmp(o.arg, "msa") == 0) type = MB_BENCH_MSA;
			else kom_assert(0, "unknown type");
		}
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minibwa bench [options] <in.mbw>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b STR         type: 2a, sa or msa [2a]\n");
		fprintf(stderr, "  -n NUM         number of data points [10k]\n");
		fprintf(stderr, "  -v INT         interval size for msa [%d]\n", intv);
		fprintf(stderr, "  -p             print results for each data point\n");
		return 1;
	}

	bwt = mb_bwt_load(argv[o.ind]);
	t = kom_cputime();
	if (type == MB_BENCH_2A) {
		for (i = 0; i < n; ++i) {
			uint64_t k = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t l = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t cntk[4], cntl[4];
			mb_bwt_rank2a(bwt, k, l, cntk, cntl);
			cs = cs * cntk[1] + cntl[0];
			if (print_val) printf("%lld\n", cntk[1]);
		}
	} else if (type == MB_BENCH_SA) {
		for (i = 0; i < n; ++i) {
			uint64_t s, k = kom_splitmix64(&x) % bwt->seq_len;
			s = mb_bwt_sa(bwt, k);
			cs = cs * 0xbf58476d1ce4e5b9ULL ^ s;
			if (print_val) printf("%lld\n", s);
		}
	} else if (type == MB_BENCH_MSA) {
		for (i = 0; i < n; ++i) {
			uint64_t j, xor = 0, k = kom_splitmix64(&x) % bwt->seq_len;
			uint64_t l = k + intv < bwt->seq_len? k + intv : bwt->seq_len;
			#if 0
			for (j = k; j < l; ++j) {
				uint64_t s = mb_bwt_sa(bwt, j);
				xor ^= s;
			}
			#else
			uint64_t sa[intv], n_sa = l - k;
			for (j = 0; j < n_sa; ++j)
				sa[j] = k + j;
			mb_bwt_sa2(0, bwt, l - k, sa);
			for (j = 0; j < n_sa; ++j)
				xor ^= sa[j];
			#endif
			cs = cs * 0xbf58476d1ce4e5b9ULL ^ xor;
			if (print_val) printf("%lld\n", xor);
		}
	}
	fprintf(stderr, "checksum = %lx\n", (unsigned long)cs);
	fprintf(stderr, "t = %.3f\n", kom_cputime() - t);
	mb_bwt_destroy(bwt);
	return 0;
}

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

int main_seed(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	ketopt_t o = KETOPT_INIT;
	gzFile fp;
	kseq_t *ks;
	int c, min_len = 19, min_occ = 1, max_size_out = 20;
	uint64_t *sa, m_a = 0;
	mb_sai_t *a = 0;
	kstring_t out = {0};

	while ((c = ketopt(&o, argc, argv, 1, "l:s:w:", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 's') min_occ = atoi(o.arg);
		else if (c == 'w') max_size_out = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa seed [options] <in.mbw> <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT     min seed length [%d]\n", min_len);
		fprintf(stderr, "  -s INT     min interval size [%d]\n", min_occ);
		return 1;
	}

	bwt = mb_bwt_load(argv[o.ind]);
	kom_assert(bwt, "failed to open the BWT file.");
	fp = strcmp(argv[o.ind+1], "-")? gzopen(argv[o.ind+1], "rb") : gzdopen(0, "rb");
	ks = kseq_init(fp);
	sa = kom_calloc(uint64_t, max_size_out);
	while (kseq_read(ks) >= 0) {
		int64_t x = 0, i, n_a = 0;
		mb_sai_t p;
		out.l = 0;
		kom_sprintf_lite(&out, "SQ\t%s\t%ld\n", ks->name.s, ks->seq.l);
		for (i = 0; i < ks->seq.l; ++i)
			ks->seq.s[i] = kom_nt4_table[(uint8_t)ks->seq.s[i]];
		do {
			x = mb_bwt_smem(0, bwt, min_len, min_occ, ks->seq.l, (uint8_t*)ks->seq.s, x, &p);
			if (p.size > 0) {
				kom_grow(mb_sai_t, a, n_a, m_a);
				a[n_a++] = p;
			}
		} while (x < ks->seq.l);
		for (i = 0; i < n_a; ++i) {
			kom_sprintf_lite(&out, "EM\t%ld\t%ld\t%ld", a[i].info>>32, a[i].info&0xffffffffull, a[i].size);
			if (a[i].size <= max_size_out) {
				int64_t j, n_sa = a[i].size;
				#if 0
				for (j = 0; j < a[i].size; ++j)
					sa[j] = mb_bwt_sa(bwt, a[i].x[0] + j);
				#else
				for (j = 0; j < a[i].size; ++j)
					sa[j] = a[i].x[0] + j;
				mb_bwt_sa2(0, bwt, a[i].size, sa);
				#endif
				for (j = 0; j < n_sa; ++j)
					kom_sprintf_lite(&out, "\t%ld", sa[j]);
			} else kom_sprintf_lite(&out, "\t*");
			kom_sprintf_lite(&out, "\n");
		}
		kom_sprintf_lite(&out, "//\n");
		fputs(out.s, stdout);
	}
	free(sa);
	kseq_destroy(ks);
	gzclose(fp);
	mb_bwt_destroy(bwt);
	return 0;
}
