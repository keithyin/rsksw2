/* Deterministic fuzz driver for the C KSW2 library.
 *
 * Used by scripts/gen_fuzz_expected.sh to produce tests/expected/fuzz.txt,
 * which tests/fuzz_against_c.rs replays with the Rust port. It covers API
 * paths that the cli.c test program cannot reach: KSW_EZ_EXTZ_ONLY,
 * KSW_EZ_GENERIC_SC, KSW_EZ_EQX, splice junc arrays and X-drop values.
 *
 * Build: gcc -O2 -msse4.1 -I ksw2 -o fuzz_ref tests/c/fuzz_ref.c ksw2/ksw2_{gg,gg2,gg2_sse,extz,extz2_sse,extd,extd2_sse,extf2_sse,exts2_sse}.c
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ksw2.h"

static uint64_t rng_state;
static uint64_t rng_next(void)
{
	uint64_t x = rng_state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	rng_state = x;
	return x;
}

static const int W_LIST[] = { -1, 0, 1, 2, 5, 20, 100 };
static const int Z_LIST[] = { -1, 0, 5, 15, 100 };
/* flag list for extz/extz2/extd/extd2 */
static const int F_LIST[] = {
	0x000,            /* nothing          */
	0x001,            /* SCORE_ONLY       */
	0x040,            /* EXTZ_ONLY        */
	0x002,            /* RIGHT            */
	0x080,            /* REV_CIGAR        */
	0x082,            /* RIGHT|REV_CIGAR  */
	0x003,            /* SCORE_ONLY|RIGHT */
	0x004,            /* GENERIC_SC       */
	0x006,            /* GENERIC_SC|RIGHT */
	0x044,            /* GENERIC_SC|EXTZ_ONLY */
	0x018,            /* APPROX_MAX|APPROX_DROP */
	0x800             /* EQX (only honoured by extd2_sse in this ksw2 version) */
};
/* flag list for exts2 */
static const int FS_LIST[] = {
	0x100, 0x200, 0x104, 0x204, 0x500, 0x600, 0x101, 0x201, 0x000
};

/* verbatim copy of the helper from cli.c (it is static there) */
static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

typedef struct {
	int qlen, tlen;
	uint8_t q[64], t[64], junc[64];
	int8_t mat[25];
	int a_, b_;
	int gq, ge, gq2, ge2;
	int w, zdrop, T, xdrop, flag, sflag, noncan, jb;
} param_t;

static void gen_case(int id, param_t *P)
{
	int i;
	rng_state = (uint64_t)id * 0x9E3779B97F4A7C15ull + 1;
	P->qlen = 1 + (int)(rng_next() % 60);
	P->tlen = 1 + (int)(rng_next() % 60);
	for (i = 0; i < P->qlen; ++i) P->q[i] = (uint8_t)(rng_next() % 5);
	for (i = 0; i < P->tlen; ++i) P->t[i] = (uint8_t)(rng_next() % 5);
	P->a_ = 1 + (int)(rng_next() % 4);
	P->b_ = -(1 + (int)(rng_next() % 5));
	ksw_gen_simple_mat(5, P->mat, (int8_t)P->a_, (int8_t)P->b_);
	P->gq  = 2 + (int)(rng_next() % 9);
	P->ge  = 1 + (int)(rng_next() % 4);
	P->gq2 = 10 + (int)(rng_next() % 11);
	P->ge2 = 1 + (int)(rng_next() % 4);
	/* A band that does not cover the alignment makes the C *scalar*
	   implementations read out of bounds in ksw_backtrack() (they lack the
	   off_end[] guard of the SSE versions) -- undefined behaviour, not a
	   meaningful oracle. Keep regenerating until the band covers the ends. */
	do P->w = W_LIST[rng_next() % 7];
	while (P->w >= 0 && abs(P->qlen - P->tlen) > P->w);
	P->zdrop = Z_LIST[rng_next() % 5];
	P->T = (int)(rng_next() % 6);
	P->flag = F_LIST[rng_next() % 12];
	P->sflag = FS_LIST[rng_next() % 9];
	P->xdrop = (rng_next() % 2) ? -1 : 8 + (int)(rng_next() % 25);
	P->noncan = (int)(rng_next() % 8);
	P->jb = (int)(rng_next() % 15);
	for (i = 0; i < P->tlen; ++i) P->junc[i] = (uint8_t)rng_next();
}

static void print_cigar(int n, const uint32_t *c)
{
	int k;
	printf(" [");
	for (k = 0; k < n; ++k)
		printf("%s%d:%u", k ? " " : "", c[k] >> 4, c[k] & 0xf);
	printf("]");
}

static void print_ez(const char *name, const ksw_extz_t *ez)
{
	printf("%s sc=%d mqe=%d mqe_t=%d mte=%d mte_q=%d max=%u max_t=%d max_q=%d z=%d re=%d",
	       name, ez->score, ez->mqe, ez->mqe_t, ez->mte, ez->mte_q,
	       ez->max, ez->max_t, ez->max_q, (int)ez->zdropped, ez->reach_end);
	print_cigar(ez->n_cigar, ez->cigar);
	printf("\n");
}

/* Pre-size the CIGAR buffer. ksw_cigar2eqx() in ksw2.h drops the (possibly
   realloc'ed) pointer returned by ksw_push_cigar(), so if the EQX expansion
   needs more capacity than the backtrack buffer has, the C library writes
   through a dangling pointer. 4096 operations are far more than any 60x60
   alignment can produce, so with this preallocation no realloc ever happens
   and the C library follows its intended semantics. */
static void ez_init(ksw_extz_t *ez)
{
	memset(ez, 0, sizeof(ksw_extz_t));
	ez->m_cigar = 4096;
	ez->cigar = (uint32_t*)malloc(4096 * sizeof(uint32_t));
}

static void run_gg(const char *name,
                   int (*fn)(void *, int, const uint8_t *, int, const uint8_t *, int8_t, const int8_t *,
                             int8_t, int8_t, int, int *, int *, uint32_t **),
                   const param_t *P, int gq, int ge)
{
	int mc = 0, nc = 0, sc;
	uint32_t *ci = 0;
	sc = fn(0, P->qlen, P->q, P->tlen, P->t, 5, P->mat, (int8_t)gq, (int8_t)ge, P->w, &mc, &nc, &ci);
	printf("%s %d", name, sc);
	print_cigar(nc, ci);
	printf("\n");
	free(ci);
}

int main(void)
{
	param_t P;
	int id;
	for (id = 0; id < 2000; ++id) {
		ksw_extz_t ez;
		gen_case(id, &P);
		printf("== %d ", id);
		{
			int i;
			for (i = 0; i < P.qlen; ++i) putchar("ACGTN"[P.q[i]]);
			putchar(' ');
			for (i = 0; i < P.tlen; ++i) putchar("ACGTN"[P.t[i]]);
		}
		printf(" q=%d e=%d q2=%d e2=%d w=%d z=%d T=%d f=%d sf=%d x=%d nc=%d jb=%d\n",
		       P.gq, P.ge, P.gq2, P.ge2, P.w, P.zdrop, P.T, P.flag, P.sflag, P.xdrop, P.noncan, P.jb);

		run_gg("gg",     ksw_gg,     &P, P.gq,  P.ge);
		run_gg("gg2",    ksw_gg2,    &P, P.gq,  P.ge);
		run_gg("gg2_sse",ksw_gg2_sse,&P, P.gq,  P.ge);

		ez_init(&ez);
		ksw_extz(0, P.qlen, P.q, P.tlen, P.t, 5, P.mat, (int8_t)P.gq, (int8_t)P.ge, P.w, P.zdrop, P.flag, &ez);
		print_ez("extz", &ez);
		free(ez.cigar);

		ez_init(&ez);
		ksw_extz2_sse(0, P.qlen, P.q, P.tlen, P.t, 5, P.mat, (int8_t)P.gq, (int8_t)P.ge, P.w, P.zdrop, P.T, P.flag, &ez);
		print_ez("extz2_sse", &ez);
		free(ez.cigar);

		ez_init(&ez);
		ksw_extd(0, P.qlen, P.q, P.tlen, P.t, 5, P.mat, (int8_t)P.gq, (int8_t)P.ge, (int8_t)P.gq2, (int8_t)P.ge2,
		         P.w, P.zdrop, P.flag, &ez);
		print_ez("extd", &ez);
		free(ez.cigar);

		ez_init(&ez);
		ksw_extd2_sse(0, P.qlen, P.q, P.tlen, P.t, 5, P.mat, (int8_t)P.gq, (int8_t)P.ge, (int8_t)P.gq2, (int8_t)P.ge2,
		              P.w, P.zdrop, P.T, P.flag, &ez);
		print_ez("extd2_sse", &ez);
		free(ez.cigar);

		ez_init(&ez);
		ksw_extf2_sse(0, P.qlen, P.q, P.tlen, P.t, (int8_t)P.a_, (int8_t)P.b_, (int8_t)P.ge, P.w, P.xdrop, &ez);
		print_ez("extf2_sse", &ez);
		free(ez.cigar);

		ez_init(&ez);
		ksw_exts2_sse(0, P.qlen, P.q, P.tlen, P.t, 5, P.mat, (int8_t)P.gq, (int8_t)P.ge, (int8_t)P.gq2,
		              (int8_t)P.noncan, P.zdrop, (int8_t)P.jb, P.sflag, P.junc, &ez);
		print_ez("exts2_sse", &ez);
		free(ez.cigar);
	}
	return 0;
}
