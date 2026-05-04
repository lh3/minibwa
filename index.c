#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "libsais64.h"
#include "kommon.h"
#include "ketopt.h"
#include "mbpriv.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

// --meth: write doubled c2t FASTA (>r<name> G->A, >f<name> C->T) next to <fa>.
// Format: uppercased, 100bp wrap. Skipped if <fa>.bwameth.c2t is fresher than <fa>.

static void meth_project_and_write(FILE *out, const char *prefix, const char *name,
                                   const char *seq, size_t len, char from, char to)
{
	char buf[65536];
	size_t bl = 0, i;
	fprintf(out, ">%s%s\n", prefix, name);
	for (i = 0; i < len; ++i) {
		if (bl + 2 > sizeof(buf)) { fwrite(buf, 1, bl, out); bl = 0; }
		char c = seq[i];
		buf[bl++] = (c == from) ? to : c;
		if (((i + 1) % 100) == 0) buf[bl++] = '\n';
	}
	if (len % 100 != 0) {
		if (bl + 1 > sizeof(buf)) { fwrite(buf, 1, bl, out); bl = 0; }
		buf[bl++] = '\n';
	}
	if (bl) fwrite(buf, 1, bl, out);
}

static int meth_c2t_is_fresh(const char *in_fa, const char *out_fa)
{
	struct stat a, b;
	if (stat(in_fa, &a) != 0 || stat(out_fa, &b) != 0) return 0;
	return b.st_mtime > a.st_mtime; // strict: same mtime (1s resolution) means rebuild
}

static int meth_write_c2t_fasta(const char *fa, char *out_fa, size_t out_fa_cap)
{
	int n = snprintf(out_fa, out_fa_cap, "%s.bwameth.c2t", fa);
	if (n <= 0 || (size_t)n >= out_fa_cap) {
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: reference path too long\n");
		return 1;
	}
	if (meth_c2t_is_fresh(fa, out_fa)) {
		if (kom_verbose >= 2)
			fprintf(stderr, "[index:--meth] %s is newer than %s; skipping c2t FASTA emission\n",
			        out_fa, fa);
		return 0;
	}
	gzFile in = gzopen(fa, "r");
	if (in == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: cannot open %s\n", fa);
		return 2;
	}
	FILE *out = fopen(out_fa, "w");
	if (out == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: cannot open %s for writing\n", out_fa);
		gzclose(in);
		return 3;
	}
	if (kom_verbose >= 2) fprintf(stderr, "[index:--meth] writing %s ...\n", out_fa);

	kseq_t *seq = kseq_init(in);
	int64_t total_bases = 0, n_seqs = 0;
	int kr = 0;
	while ((kr = kseq_read(seq)) >= 0) {
		size_t i;
		// upper-case before projection (soft-masked input round-trips like bwameth.py)
		for (i = 0; i < seq->seq.l; ++i) {
			char c = seq->seq.s[i];
			if (c >= 'a' && c <= 'z') seq->seq.s[i] = (char)(c - 'a' + 'A');
		}
		meth_project_and_write(out, "r", seq->name.s, seq->seq.s, seq->seq.l, 'G', 'A');
		meth_project_and_write(out, "f", seq->name.s, seq->seq.s, seq->seq.l, 'C', 'T');
		total_bases += (int64_t)seq->seq.l;
		++n_seqs;
	}
	kseq_destroy(seq);
	gzclose(in);
	// kseq_read: -1 = EOF, < -1 = parse/IO error; drop partial file so freshness check rebuilds
	if (kr < -1) {
		fclose(out);
		unlink(out_fa);
		if (kom_verbose >= 1)
			fprintf(stderr, "ERROR: failed while reading %s (kseq_read=%d)\n", fa, kr);
		return 4;
	}
	if (fclose(out) != 0) {
		unlink(out_fa);
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: failed to close %s\n", out_fa);
		return 4;
	}
	if (kom_verbose >= 2)
		fprintf(stderr, "[index:--meth] emitted %ld seqs, %ld bp (doubled to %ld bp of c2t text)\n",
		        (long)n_seqs, (long)total_bases, (long)(2 * total_bases));
	return 0;
}

static mb_bwt_t *mb_bwt_libsais(const l2b_t *l2b, int sa_bit, int both_strand, int n_thread)
{
	const int fs = 10000;
	uint8_t *seq;
	int64_t i, j, *a, primary, len;
	mb_bwt_t *bwt;
	uint64_t *ssa, n_ssa, mask;

	len = both_strand? l2b->tot_len * 2 : l2b->tot_len;
	seq = kom_malloc(uint8_t, len);
	a = kom_malloc(int64_t, len + fs + 1);
	for (i = 0, j = 0; i < l2b->tot_len; ++i, ++j)
		seq[j] = l2b_get0(l2b, i);
	if (both_strand)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++j)
			seq[j] = 3 - l2b_get0(l2b, i);
#ifdef LIBSAIS_OPENMP
    libsais64_omp(seq, a + 1, len, fs, 0, n_thread);
#else
    libsais64(seq, a + 1, len, fs, 0);
#endif
	a[0] = len; // libsais doesn't write a[0], which always equals to len

	n_ssa = (len + (1<<sa_bit)) >> sa_bit;
	ssa = kom_calloc(uint64_t, n_ssa);
	mask = (1<<sa_bit) - 1;
	for (i = 0; i <= len; ++i)
		if ((i & mask) == 0)
			ssa[i >> sa_bit] = a[i];
	ssa[0] = (uint64_t)-1;
	primary = (uint64_t)-1;
	for (i = 0; i <= len; ++i) {
		if (a[i] == 0) primary = i;
		else a[i] = seq[a[i] - 1];
	}
	assert(primary != (uint64_t)-1);
	for (i = 0; i < primary; ++i) seq[i] = a[i];
	for (; i < len; ++i) seq[i] = a[i + 1];
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
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};
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
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};
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
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};
	while ((c = ketopt(&o, argc, argv, 1, "1u:t:", long_opts)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == '1') both_strand = 0;
		else if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 901) return usage_genbwt(stdout, sa_bit, n_thread);
	}
	if (argc - o.ind < 2) return usage_genbwt(stderr, sa_bit, n_thread);
	l2b = l2b_load(argv[o.ind]);
	kom_assert(l2b, "failed to open the input file.");
	bwt = mb_bwt_libsais(l2b, sa_bit, both_strand, n_thread);
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
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ 0, 0, 0 }
	};
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
	fprintf(fp, "  --meth    build a bwameth-style doubled c2t reference (writes\n");
	fprintf(fp, "            <in.fasta>.bwameth.c2t and indexes that). Use with\n");
	fprintf(fp, "            `minibwa map --meth <in.fasta> R1.fq [R2.fq]`.\n");
	fprintf(fp, "  --help    print this help message\n");
	return fp == stdout? 0 : 1;
}

int main_index(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, low_mem = 0, n_thread = 4, sa_bit = 4, meth = 0;
	int64_t block_size = 10000000;
	uint64_t seed = 11;
	char *prefix, *fn_l2b, *fn_bwt;
	const char *fa_in;
	char meth_c2t_path[4096];
	l2b_t *l2b;
	mb_bwt_t *bwt;
	static ko_longopt_t long_opts[] = {
		{ "help", ko_no_argument, 901 },
		{ "meth", ko_no_argument, 1000 },
		{ 0, 0, 0 }
	};

	while ((c = ketopt(&o, argc, argv, 1, "ls:u:b:t:", long_opts)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'l') low_mem = 1;
		else if (c == 'b') block_size = kom_parse_num(o.arg, 0);
		else if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 's') seed = atol(o.arg);
		else if (c == 901) return usage_index(stdout, seed, sa_bit, n_thread);
		else if (c == 1000) meth = 1;
	}
	if (argc - o.ind == 0) return usage_index(stderr, seed, sa_bit, n_thread);

	fa_in = argv[o.ind];
	if (meth) {
		/* Reject -p / explicit out-prefix for --meth: the c2t flow is
		 * keyed off "<in.fasta>.bwameth.c2t" both at index time and at
		 * map time, and silently overriding it would only confuse. */
		if (o.ind + 1 < argc) {
			if (kom_verbose >= 1)
				fprintf(stderr, "ERROR: --meth does not accept an explicit out-prefix (use bare <in.fasta>)\n");
			return 1;
		}
		int rc = meth_write_c2t_fasta(fa_in, meth_c2t_path, sizeof(meth_c2t_path));
		if (rc != 0) return rc;
		fa_in = meth_c2t_path; /* index the doubled FASTA */
	}

	prefix = (!meth && o.ind + 1 < argc) ? argv[o.ind+1] : (char *)fa_in;
	fn_l2b = kom_calloc(char, strlen(prefix) + 5);
	strcat(strcpy(fn_l2b, prefix), ".l2b");
	fn_bwt = kom_calloc(char, strlen(prefix) + 5);
	strcat(strcpy(fn_bwt, prefix), ".mbw");

	l2b = l2b_import(fa_in, seed);
	kom_assert(l2b, "failed to read the genome FASTA.");
	if (low_mem) {
#ifdef USE_GPL
		l2b_save_pac(fn_l2b, l2b, 1);
		mb_bwtgen(fn_l2b, fn_bwt, block_size);
		l2b_save(fn_l2b, l2b);
		bwt = mb_bwt_load_raw(fn_bwt);
		mb_bwt_gen_sa(bwt, sa_bit);
#else
		if (kom_verbose >= 1) fprintf(stderr, "ERROR: option -l not compiled as it depends on GPL'd code\n");
		abort();
#endif
	} else {
		l2b_save(fn_l2b, l2b);
		bwt = mb_bwt_libsais(l2b, sa_bit, 1, n_thread);
	}
	mb_bwt_save(fn_bwt, bwt);
	l2b_destroy(l2b);
	mb_bwt_destroy(bwt);
	free(fn_bwt); free(fn_l2b);
	return 0;
}
