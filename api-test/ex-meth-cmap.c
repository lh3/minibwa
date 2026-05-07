/* Smoke tests for mb_meth_cmap_build: f<X>/r<X> partner detection,
 * 3+-collision poisoning, length-mismatch fallback. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l2bit.h"
#include "meth.h"

static l2b_t make_l2b(int n, const char **names, const uint64_t *lens)
{
	l2b_t l = {0};
	int i;
	l.n_ctg = n;
	l.ctg = calloc(n, sizeof(l2b_ctg_t));
	for (i = 0; i < n; ++i) {
		l.ctg[i].name = strdup(names[i]);
		l.ctg[i].len  = lens[i];
		l.ctg[i].off  = (uint64_t)i * 1000;
	}
	return l;
}

static void free_l2b(l2b_t *l)
{
	uint64_t i;
	for (i = 0; i < l->n_ctg; ++i) free(l->ctg[i].name);
	free(l->ctg);
}

static int n_pass, n_fail;
#define CHECK(cond, msg) do { \
	if (cond) { ++n_pass; } \
	else { ++n_fail; fprintf(stderr, "  FAIL: %s (line %d)\n", msg, __LINE__); } \
} while (0)

static void test_no_doubled(void)
{
	const char *names[] = {"chr1", "chr2"};
	const uint64_t lens[] = {1000, 2000};
	l2b_t l = make_l2b(2, names, lens);
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	fprintf(stderr, "test_no_doubled:\n");
	CHECK(m->n_internal == 2, "n_internal");
	CHECK(m->n_output == 2, "n_output");
	CHECK(m->canon_tid[0] == 0 && m->canon_tid[1] == 1, "canon_tid identity");
	CHECK(m->paired_tid[0] == -1 && m->paired_tid[1] == -1, "no partners");
	CHECK(m->direction[0] == 0 && m->direction[1] == 0, "no direction");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

static void test_paired_fr(void)
{
	const char *names[] = {"fchr1", "rchr1"};
	const uint64_t lens[] = {1000, 1000};
	l2b_t l = make_l2b(2, names, lens);
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	fprintf(stderr, "test_paired_fr:\n");
	CHECK(m->n_internal == 2 && m->n_output == 1, "collapse");
	CHECK(m->canon_tid[0] == 0 && m->canon_tid[1] == 0, "canon=0 for both");
	CHECK(m->paired_tid[0] == 1 && m->paired_tid[1] == 0, "partners cross-link");
	CHECK(m->direction[0] == 'f' && m->direction[1] == 'r', "direction");
	CHECK(m->out_tid[0] == 0 && m->out_tid[1] == 0, "shared out_tid");
	CHECK(strcmp(m->output_names[0], "chr1") == 0, "stripped name");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

static void test_three_collision(void)
{
	/* All three names must strip to the same string ("chr1") to share an out_tid;
	 * use two f-prefixed and one r-prefixed dupes. */
	const char *names[] = {"fchr1", "rchr1", "fchr1"};
	const uint64_t lens[] = {1000, 1000, 1000};
	l2b_t l = make_l2b(3, names, lens);
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	fprintf(stderr, "test_three_collision:\n");
	CHECK(m->n_output == 1, "single out_tid");
	CHECK(m->paired_tid[0] == -1 && m->paired_tid[1] == -1 && m->paired_tid[2] == -1,
	      "all partners -1 after 3+ collision");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

static void test_four_collision(void)
{
	const char *names[] = {"fchr1", "rchr1", "fchr1", "rchr1"};
	const uint64_t lens[] = {1000, 1000, 1000, 1000};
	l2b_t l = make_l2b(4, names, lens);
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	fprintf(stderr, "test_four_collision:\n");
	CHECK(m->n_output == 1, "single out_tid");
	/* All four must have paired_tid == -1: regression for the toggling-after-poison bug */
	CHECK(m->paired_tid[0] == -1, "partner[0] -1 after 4-collision");
	CHECK(m->paired_tid[1] == -1, "partner[1] -1 after 4-collision");
	CHECK(m->paired_tid[2] == -1, "partner[2] -1 after 4-collision");
	CHECK(m->paired_tid[3] == -1, "partner[3] -1 after 4-collision");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

static void test_length_mismatch(void)
{
	const char *names[] = {"fchr1", "rchr1"};
	const uint64_t lens[] = {1000, 1001}; /* mismatched */
	l2b_t l = make_l2b(2, names, lens);
	fprintf(stderr, "test_length_mismatch:\n");
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	CHECK(m != NULL, "build does not abort on mismatch");
	CHECK(m->paired_tid[0] == -1 && m->paired_tid[1] == -1, "mismatch -> no partner");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

static void test_mixed_doubled_and_plain(void)
{
	const char *names[] = {"fchr1", "rchr1", "chrM", "fchr2", "rchr2"};
	const uint64_t lens[] = {1000, 1000, 16569, 2000, 2000};
	l2b_t l = make_l2b(5, names, lens);
	mb_meth_cmap_t *m = mb_meth_cmap_build(&l);
	fprintf(stderr, "test_mixed_doubled_and_plain:\n");
	CHECK(m->n_output == 3, "chr1, chrM, chr2");
	CHECK(m->paired_tid[0] == 1 && m->paired_tid[1] == 0, "chr1 partners");
	CHECK(m->paired_tid[2] == -1, "chrM no partner");
	CHECK(m->paired_tid[3] == 4 && m->paired_tid[4] == 3, "chr2 partners");
	CHECK(m->direction[2] == 0, "chrM direction 0");
	mb_meth_cmap_free(m);
	free_l2b(&l);
}

int main(void)
{
	test_no_doubled();
	test_paired_fr();
	test_three_collision();
	test_four_collision();
	test_length_mismatch();
	test_mixed_doubled_and_plain();
	fprintf(stderr, "\n%d passed, %d failed\n", n_pass, n_fail);
	return n_fail == 0? 0 : 1;
}
