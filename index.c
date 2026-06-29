#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "libsais.h"
#include "libsais64.h"
#include "kommon.h"
#include "ketopt.h"
#include "mbpriv.h"

void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

static ko_longopt_t long_opts[] = { // common long options shared across all index-related functions
	{ "help", ko_no_argument, 901 },
	{ "meth", ko_no_argument, 902 },
	{ 0, 0, 0 }
};

static inline uint8_t l2b_c2t(uint8_t b) { return b == 1? 3 : b; } // C(1) -> T(3)
static inline uint8_t l2b_g2a(uint8_t b) { return b == 2? 0 : b; } // G(2) -> A(0)

// invert the suffix array a[] (32- or 64-bit) to the BWT in seq[], sample the SSA and drop the primary ($); return the primary
static int64_t sa_to_bwt(void *a, int use_int32, uint8_t *seq, int64_t len, int sa_bit, uint64_t *ssa)
{
	int32_t *a32 = a;
	int64_t *a64 = a, i, primary = -1;
	uint64_t mask = (1ULL<<sa_bit) - 1;
	for (i = 0; i <= len; ++i) {
		int64_t v = use_int32? a32[i] : a64[i];
		if ((i & mask) == 0) ssa[i>>sa_bit] = v;
		if (v == 0) primary = i;
		else if (use_int32) a32[i] = seq[v - 1];
		else a64[i] = seq[v - 1];
	}
	ssa[0] = (uint64_t)-1;
	for (i = 0; i < primary; ++i) seq[i] = use_int32? a32[i] : a64[i];
	for (; i < len; ++i) seq[i] = use_int32? a32[i+1] : a64[i+1];
	return primary;
}

static mb_bwt_t *mb_bwt_libsais(const l2b_t *l2b, int sa_bit, int both_strand, int is_meth, int n_thread)
{
	const int fs = 10000;
	uint8_t *seq;
	int64_t i, j, primary, len;
	mb_bwt_t *bwt;
	uint64_t *ssa, n_ssa;
	void *a;
	int use_int32;

	len = l2b->tot_len * (is_meth? 2 : 1) * (both_strand? 2 : 1);
	// use a 32-bit suffix array (half the memory) when the concatenated length fits in int32_t
	use_int32 = (len + fs + 1 <= INT32_MAX);
	seq = kom_malloc(uint8_t, len);
	if (use_int32) a = kom_malloc(int32_t, len + fs + 1);
	else a = kom_malloc(int64_t, len + fs + 1);
	if (is_meth) {
		// c2t forward
		for (i = 0, j = 0; i < l2b->tot_len; ++i, ++j)
			seq[j] = l2b_c2t(l2b_get0(l2b, i));
		// g2a forward
		for (i = 0; i < l2b->tot_len; ++i, ++j)
			seq[j] = l2b_g2a(l2b_get0(l2b, i));
		if (both_strand) {
			// g2a reverse (reverse complement of g2a converted)
			for (i = l2b->tot_len - 1; i >= 0; --i, ++j)
				seq[j] = 3 - l2b_g2a(l2b_get0(l2b, i));
			// c2t reverse (reverse complement of c2t converted)
			for (i = l2b->tot_len - 1; i >= 0; --i, ++j)
				seq[j] = 3 - l2b_c2t(l2b_get0(l2b, i));
		}
	} else {
		for (i = 0, j = 0; i < l2b->tot_len; ++i, ++j)
			seq[j] = l2b_get0(l2b, i);
		if (both_strand)
			for (i = l2b->tot_len - 1; i >= 0; --i, ++j)
				seq[j] = 3 - l2b_get0(l2b, i);
	}
	if (use_int32) {
		int32_t *a32 = a;
#ifdef LIBSAIS_OPENMP
		libsais_omp(seq, a32 + 1, (int32_t)len, fs, 0, n_thread);
#else
		libsais(seq, a32 + 1, (int32_t)len, fs, 0);
#endif
		a32[0] = (int32_t)len; // libsais doesn't write a[0], which always equals to len
	} else {
		int64_t *a64 = a;
#ifdef LIBSAIS_OPENMP
		libsais64_omp(seq, a64 + 1, len, fs, 0, n_thread);
#else
		libsais64(seq, a64 + 1, len, fs, 0);
#endif
		a64[0] = len;
	}

	n_ssa = (len + (1<<sa_bit)) >> sa_bit;
	ssa = kom_calloc(uint64_t, n_ssa);
	primary = sa_to_bwt(a, use_int32, seq, len, sa_bit, ssa);
	assert(primary != -1);
	free(a);
	bwt = mb_bwt_init_from_raw(1, seq, len, primary);
	bwt->sa_bit = sa_bit, bwt->n_sa = n_ssa, bwt->sa = ssa;
	free(seq);
	return bwt;
}

static int usage_fa2bit(FILE *fp, uint64_t seed)
{
	fprintf(fp, "Usage: minibwa fa2bit [options] <in.fa> <out.l2b>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -s INT    random seed [%lu]\n", (unsigned long)seed);
	fprintf(fp, "  -p        output the BWA pac format\n");
	fprintf(fp, "  -2        output both strands (effective with -p)\n");
	fprintf(fp, "  --help    print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_fa2bit(int argc, char *argv[])
{
	l2b_t *l2b;
	int out_pac = 0, both_strand = 0;
	uint64_t seed = 11;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "s:p2", long_opts)) >= 0) {
		if (c == 's') seed = atol(o.arg);
		else if (c == 'p') out_pac = 1;
		else if (c == '2') both_strand = 1;
		else if (c == 901) return usage_fa2bit(stdout, seed);
	}
	if (argc - o.ind < 2) return usage_fa2bit(stderr, seed);
	l2b = l2b_import(argv[o.ind], seed);
	if (out_pac)
		l2b_save_pac(argv[o.ind+1], l2b, both_strand);
	else
		l2b_save(argv[o.ind+1], l2b);
	l2b_destroy(l2b);
	return 0;
}

#ifdef USE_GPL
static int usage_genraw(FILE *fp)
{
	fprintf(fp, "Usage: minibwa genraw [options] <in.pac> <out.raw-bwt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -b NUM      block size [10m]\n");
	fprintf(fp, "  --help      print this help message\n");
	return fp == stdout? 0 : 1;
}
#endif

int main_genraw(int argc, char *argv[])
{
#ifdef USE_GPL
	ketopt_t o = KETOPT_INIT;
	int c, block_size = 10000000;
	while ((c = ketopt(&o, argc, argv, 1, "b:", long_opts)) >= 0) {
		if (c == 'b') block_size = kom_parse_num(o.arg, 0);
		else if (c == 901) return usage_genraw(stdout);
	}
	if (argc - o.ind < 2) return usage_genraw(stderr);
	mb_bwtgen(argv[o.ind], argv[o.ind+1], block_size);
	return 0;
#else
	(void)argc; (void)argv;
	if (kom_verbose >= 1) fprintf(stderr, "ERROR: genraw not compiled as it depends on GPL'd code\n");
	return 1;
#endif
}

static int usage_raw2bwt(FILE *fp)
{
	fprintf(fp, "Usage: minibwa raw2bwt <raw.bwt> <recode.bwt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  --help    print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_raw2bwt(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	int i;
	for (i = 1; i < argc; ++i)
		if (strcmp(argv[i], "--help") == 0) return usage_raw2bwt(stdout);
	if (argc < 3) return usage_raw2bwt(stderr);
	bwt = mb_bwt_load_raw(argv[1]);
	mb_bwt_save(argv[2], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

static int usage_genbwt(FILE *fp, int sa_bit, int n_thread)
{
	(void)n_thread;
	fprintf(fp, "Usage: minibwa genbwt [options] <in.l2b> <out.bwt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -u INT      SA sample rate at 1/(1<<INT) [%d]\n", sa_bit);
	fprintf(fp, "  -1          forward strand only\n");
#ifdef LIBSAIS_OPENMP
	fprintf(fp, "  -t INT      number of threads [%d]\n", n_thread);
#endif
	fprintf(fp, "  --help      print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_genbwt(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, n_thread = 4, both_strand = 1, sa_bit = 4;
	mb_bwt_t *bwt;
	l2b_t *l2b;
	while ((c = ketopt(&o, argc, argv, 1, "1u:t:", long_opts)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == '1') both_strand = 0;
		else if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 901) return usage_genbwt(stdout, sa_bit, n_thread);
	}
	if (argc - o.ind < 2) return usage_genbwt(stderr, sa_bit, n_thread);
	l2b = l2b_load(argv[o.ind]);
	kom_assert(l2b, "failed to open the input file.");
	bwt = mb_bwt_libsais(l2b, sa_bit, both_strand, 0, n_thread);
	l2b_destroy(l2b);
	mb_bwt_save(argv[o.ind+1], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

static int usage_gensa(FILE *fp, int sa_bit)
{
	fprintf(fp, "Usage: minibwa gensa [options] <in.bwt> <out.bwt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -u INT    sample rate at 1/(1<<INT) [%d]\n", sa_bit);
	fprintf(fp, "  -r        input BWT in the raw BWA format\n");
	fprintf(fp, "  --help    print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_gensa(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	int c, sa_bit = 4, is_raw = 0;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "ru:", long_opts)) >= 0) {
		if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 'r') is_raw = 1;
		else if (c == 901) return usage_gensa(stdout, sa_bit);
	}
	if (argc - o.ind < 2) return usage_gensa(stderr, sa_bit);
	bwt = is_raw? mb_bwt_load_raw(argv[o.ind]) : mb_bwt_load(argv[o.ind]);
	mb_bwt_gen_sa(bwt, sa_bit);
	mb_bwt_save(argv[o.ind+1], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

static int usage_index(FILE *fp, uint64_t seed, int sa_bit, int n_thread)
{
	(void)n_thread;
	fprintf(fp, "Usage: minibwa index [options] <in.fasta> [out.prefix]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -s INT    random seed for amibiguous bases [%ld]\n", (unsigned long)seed);
	fprintf(fp, "  -u INT    SA sample rate at 1/(1<<INT) [%d]\n", sa_bit);
	fprintf(fp, "  -l        low-memory GPL'd algorithm for BWT construction\n");
	fprintf(fp, "  -b NUM    block size (effective with -l) [10m]\n");
#ifdef LIBSAIS_OPENMP
	fprintf(fp, "  -t INT    number of threads (effective w/o -l) [%d]\n", n_thread);
#endif
	fprintf(fp, "  --meth    build FM-index for BS-seq mapping\n");
	fprintf(fp, "  --help    print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_index(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, low_mem = 0, n_thread = 4, sa_bit = 4, is_meth = 0;
	int64_t block_size = 10000000;
	uint64_t seed = 11;
	char *prefix, *fn_l2b, *fn_bwt, *fn_meth_bwt = 0;
	l2b_t *l2b;
	mb_bwt_t *bwt;

	while ((c = ketopt(&o, argc, argv, 1, "ls:u:b:t:", long_opts)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'l') low_mem = 1;
		else if (c == 'b') block_size = kom_parse_num(o.arg, 0);
		else if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 's') seed = atol(o.arg);
		else if (c == 901) return usage_index(stdout, seed, sa_bit, n_thread);
		else if (c == 902) is_meth = 1;
	}
	if (argc - o.ind == 0) return usage_index(stderr, seed, sa_bit, n_thread);

	prefix = o.ind + 1 < argc? argv[o.ind+1] : argv[o.ind];
	fn_l2b = kom_calloc(char, strlen(prefix) + 10);
	strcat(strcpy(fn_l2b, prefix), ".l2b");
	fn_bwt = kom_calloc(char, strlen(prefix) + 10);
	strcat(strcpy(fn_bwt, prefix), ".mbw");
	if (is_meth) {
		fn_meth_bwt = kom_calloc(char, strlen(prefix) + 10);
		strcat(strcpy(fn_meth_bwt, prefix), ".meth.mbw");
	}

	l2b = l2b_import(argv[o.ind], seed);
	kom_assert(l2b, "failed to read the genome FASTA.");
	if (low_mem) {
#ifdef USE_GPL
		l2b_save_pac(fn_l2b, l2b, 1);
		mb_bwtgen(fn_l2b, fn_bwt, block_size);
		l2b_save(fn_l2b, l2b);
		bwt = mb_bwt_load_raw(fn_bwt);
		mb_bwt_gen_sa(bwt, sa_bit);
		mb_bwt_save(fn_bwt, bwt);
		mb_bwt_destroy(bwt);
		if (is_meth) {
			l2b_save_pac_meth(fn_l2b, l2b, 1);
			mb_bwtgen(fn_l2b, fn_meth_bwt, block_size);
			l2b_save(fn_l2b, l2b); // restore the real .l2b; the meth pac above overwrote it (cf. the regular pass)
			bwt = mb_bwt_load_raw(fn_meth_bwt);
			mb_bwt_gen_sa(bwt, sa_bit);
			mb_bwt_save(fn_meth_bwt, bwt);
			mb_bwt_destroy(bwt);
		}
#else
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: option -l not compiled as it depends on GPL'd code\n");
		abort();
#endif
	} else {
		l2b_save(fn_l2b, l2b);
		bwt = mb_bwt_libsais(l2b, sa_bit, 1, 0, n_thread);
		mb_bwt_save(fn_bwt, bwt);
		mb_bwt_destroy(bwt);
		if (is_meth) {
			bwt = mb_bwt_libsais(l2b, sa_bit, 1, 1, n_thread);
			mb_bwt_save(fn_meth_bwt, bwt);
			mb_bwt_destroy(bwt);
		}
	}
	l2b_destroy(l2b);
	free(fn_meth_bwt); free(fn_bwt); free(fn_l2b);
	return 0;
}
