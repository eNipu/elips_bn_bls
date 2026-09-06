/*
 * ELiPS known-answer vector runner  (issue #12)
 *
 * Reads a .vec file produced by tools/reference/gen_vectors.py and drives the
 * library through each record, comparing against the reference result.
 *
 * Exit status is the whole point: 0 only if every record matched. The existing
 * test drivers return 0 unconditionally, so nothing in this repository was
 * CI-gradeable before this file existed.
 *
 * The runner takes the prime from the vector file rather than from one of the
 * library's init_* routines. That keeps it curve-agnostic and, usefully, lets
 * the field layer be exercised on BLS12-381 even though the library has no
 * BLS12-381 curve support yet.
 *
 * Canonicality: the library has routines that return non-canonical residues.
 * Fp_neg computes p - a, so Fp_neg(0) yields p rather than 0, and Fp_set_neg
 * uses mpz_neg and can yield a negative value. Those are congruent to the right
 * answer but are not reduced, and Fp_cmp compares integers, not residues. The
 * runner therefore compares modulo p and counts non-canonical outputs
 * separately, so a representation problem is visible without being conflated
 * with a wrong value.
 */

/* getline, strtok_r and ssize_t are POSIX, and the build compiles with
 * -std=c11 (no GNU extensions), which hides them on glibc unless asked for.
 * Without this the two functions are implicitly declared as returning int,
 * which truncates their pointer results on any 64-bit target. */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include <ELiPS_bn_bls/bn_fp.h>
#include <ELiPS_bn_bls/bn_fp2.h>
#include <ELiPS_bn_bls/bn_fp6.h>
#include <ELiPS_bn_bls/bn_fp12.h>
#include <ELiPS_bn_bls/bn_efp.h>
#include <ELiPS_bn_bls/bn_efp2.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>
#include <ELiPS_bn_bls/curve_settings.h>

#define MAX_TOK 40

static mpz_t PRIME;
static long n_pass, n_fail, n_noncanon;
static char cur_op[64];

/* ------------------------------------------------------------ token helpers */

typedef struct { mpz_t v[MAX_TOK]; int n; } Toks;

static void toks_init(Toks *t)  { for (int i = 0; i < MAX_TOK; i++) mpz_init(t->v[i]); t->n = 0; }
static void toks_clear(Toks *t) { for (int i = 0; i < MAX_TOK; i++) mpz_clear(t->v[i]); }

/* Report one mismatch with enough context to locate it. */
static void fail(long line, int idx, const mpz_t got, const mpz_t want)
{
    n_fail++;
    if (n_fail <= 12) {
        fprintf(stderr, "MISMATCH line %ld op %s token %d\n", line, cur_op, idx);
        gmp_fprintf(stderr, "   expected %Zx\n   actual   %Zx\n", want, got);
    }
}

/* Compare one produced value against the reference, modulo p. */
static void cmp1(long line, int idx, mpz_t got, const mpz_t want)
{
    if (mpz_sgn(got) < 0 || mpz_cmp(got, PRIME) >= 0)
        n_noncanon++;
    mpz_mod(got, got, PRIME);
    if (mpz_cmp(got, want) != 0) fail(line, idx, got, want);
}

/* ---------------------------------------------------- struct <-> token glue */

static void ld_fp (Fp  *a, mpz_t *t) { mpz_set(a->x0, t[0]); }
static void ld_fp2(Fp2 *a, mpz_t *t) { ld_fp(&a->x0, t); ld_fp(&a->x1, t+1); }
static void ld_fp6(Fp6 *a, mpz_t *t) { ld_fp2(&a->x0, t); ld_fp2(&a->x1, t+2); ld_fp2(&a->x2, t+4); }
static void ld_fp12(Fp12 *a, mpz_t *t){ ld_fp6(&a->x0, t); ld_fp6(&a->x1, t+6); }

static void st_fp  (Fp  *a, mpz_t *o) { mpz_set(o[0], a->x0); }
static void st_fp2 (Fp2 *a, mpz_t *o) { st_fp(&a->x0, o); st_fp(&a->x1, o+1); }
static void st_fp6 (Fp6 *a, mpz_t *o) { st_fp2(&a->x0, o); st_fp2(&a->x1, o+2); st_fp2(&a->x2, o+4); }
static void st_fp12(Fp12 *a, mpz_t *o){ st_fp6(&a->x0, o); st_fp6(&a->x1, o+6); }

/* EFp point tokens are (infinity flag, x, y); EFp2 is (flag, x0, x1, y0, y1). */
static void ld_efp (EFp  *P, mpz_t *t) { P->infinity = mpz_get_ui(t[0]); ld_fp(&P->x, t+1); ld_fp(&P->y, t+2); }
static void ld_efp2(EFp2 *P, mpz_t *t) { P->infinity = mpz_get_ui(t[0]); ld_fp2(&P->x, t+1); ld_fp2(&P->y, t+3); }

/* ------------------------------------------------------------------ dispatch */

/* in_n / out_n are token counts; the runner uses them to slice each record. */
typedef struct {
    const char *name;
    int in_n, out_n;
    void (*run)(mpz_t *in, mpz_t *out);
} Op;

#define FP_BIN(nm, fn)                                                        \
static void nm(mpz_t *in, mpz_t *out) {                                       \
    Fp a, b, r; Fp_init(&a); Fp_init(&b); Fp_init(&r);                        \
    ld_fp(&a, in); ld_fp(&b, in+1); fn(&r, &a, &b); st_fp(&r, out);           \
    Fp_clear(&a); Fp_clear(&b); Fp_clear(&r); }

#define FP_UN(nm, fn)                                                         \
static void nm(mpz_t *in, mpz_t *out) {                                       \
    Fp a, r; Fp_init(&a); Fp_init(&r);                                        \
    ld_fp(&a, in); fn(&r, &a); st_fp(&r, out);                                \
    Fp_clear(&a); Fp_clear(&r); }

#define FPX_BIN(nm, T, fn, init, clr, ld, st)                                 \
static void nm(mpz_t *in, mpz_t *out) {                                       \
    T a, b, r; init(&a); init(&b); init(&r);                                  \
    ld(&a, in); ld(&b, in + (sizeof(#T)?0:0) + FPX_W_##T); fn(&r, &a, &b);    \
    st(&r, out); clr(&a); clr(&b); clr(&r); }

#define FPX_UN(nm, T, fn, init, clr, ld, st)                                  \
static void nm(mpz_t *in, mpz_t *out) {                                       \
    T a, r; init(&a); init(&r);                                               \
    ld(&a, in); fn(&r, &a); st(&r, out); clr(&a); clr(&r); }

#define FPX_W_Fp2  2
#define FPX_W_Fp6  6
#define FPX_W_Fp12 12

FP_BIN(op_fp_add, Fp_add)
FP_BIN(op_fp_sub, Fp_sub)
FP_BIN(op_fp_mul, Fp_mul)
FP_UN (op_fp_inv, Fp_inv)
FP_UN (op_fp_neg, Fp_neg)

FPX_BIN(op_fp2_add, Fp2, Fp2_add, Fp2_init, Fp2_clear, ld_fp2, st_fp2)
FPX_BIN(op_fp2_sub, Fp2, Fp2_sub, Fp2_init, Fp2_clear, ld_fp2, st_fp2)
FPX_BIN(op_fp2_mul, Fp2, Fp2_mul, Fp2_init, Fp2_clear, ld_fp2, st_fp2)
FPX_UN (op_fp2_sqr, Fp2, Fp2_sqr, Fp2_init, Fp2_clear, ld_fp2, st_fp2)
FPX_UN (op_fp2_inv, Fp2, Fp2_inv, Fp2_init, Fp2_clear, ld_fp2, st_fp2)
FPX_UN (op_fp2_mulbasis, Fp2, Fp2_mul_basis, Fp2_init, Fp2_clear, ld_fp2, st_fp2)

FPX_BIN(op_fp6_add, Fp6, Fp6_add, Fp6_init, Fp6_clear, ld_fp6, st_fp6)
FPX_BIN(op_fp6_sub, Fp6, Fp6_sub, Fp6_init, Fp6_clear, ld_fp6, st_fp6)
FPX_BIN(op_fp6_mul, Fp6, Fp6_mul, Fp6_init, Fp6_clear, ld_fp6, st_fp6)
FPX_UN (op_fp6_sqr, Fp6, Fp6_sqr, Fp6_init, Fp6_clear, ld_fp6, st_fp6)
FPX_UN (op_fp6_inv, Fp6, Fp6_inv, Fp6_init, Fp6_clear, ld_fp6, st_fp6)
FPX_UN (op_fp6_mulbasis, Fp6, Fp6_mul_basis, Fp6_init, Fp6_clear, ld_fp6, st_fp6)

FPX_BIN(op_fp12_add, Fp12, Fp12_add, Fp12_init, Fp12_clear, ld_fp12, st_fp12)
FPX_BIN(op_fp12_sub, Fp12, Fp12_sub, Fp12_init, Fp12_clear, ld_fp12, st_fp12)
FPX_BIN(op_fp12_mul, Fp12, Fp12_mul, Fp12_init, Fp12_clear, ld_fp12, st_fp12)
FPX_UN (op_fp12_sqr, Fp12, Fp12_sqr, Fp12_init, Fp12_clear, ld_fp12, st_fp12)
FPX_UN (op_fp12_inv, Fp12, Fp12_inv, Fp12_init, Fp12_clear, ld_fp12, st_fp12)

/* Points: the infinity flag is compared as a plain integer, coordinates mod p. */
static void st_efp(EFp *P, mpz_t *o)
{
    mpz_set_ui(o[0], (unsigned long)(P->infinity != 0));
    if (P->infinity) { mpz_set_ui(o[1], 0); mpz_set_ui(o[2], 0); }
    else             { st_fp(&P->x, o+1); st_fp(&P->y, o+2); }
}
static void st_efp2(EFp2 *P, mpz_t *o)
{
    mpz_set_ui(o[0], (unsigned long)(P->infinity != 0));
    if (P->infinity) { for (int i = 1; i < 5; i++) mpz_set_ui(o[i], 0); }
    else             { st_fp2(&P->x, o+1); st_fp2(&P->y, o+3); }
}

static void op_efp_dbl(mpz_t *in, mpz_t *out)
{ EFp P, R; EFp_init(&P); EFp_init(&R); ld_efp(&P, in);
  if (P.infinity) R.infinity = 1; else EFp_ECD(&R, &P);
  st_efp(&R, out); EFp_clear(&P); EFp_clear(&R); }

static void op_efp_add(mpz_t *in, mpz_t *out)
{ EFp P, Q, R; EFp_init(&P); EFp_init(&Q); EFp_init(&R);
  ld_efp(&P, in); ld_efp(&Q, in+3); EFp_ECA(&R, &P, &Q);
  st_efp(&R, out); EFp_clear(&P); EFp_clear(&Q); EFp_clear(&R); }

static void op_efp_mul(mpz_t *in, mpz_t *out)
{ EFp P, R; EFp_init(&P); EFp_init(&R); ld_efp(&P, in);
  EFp_SCM(&R, &P, in[3]); st_efp(&R, out); EFp_clear(&P); EFp_clear(&R); }

static void op_efp2_dbl(mpz_t *in, mpz_t *out)
{ EFp2 P, R; EFp2_init(&P); EFp2_init(&R); ld_efp2(&P, in);
  if (P.infinity) R.infinity = 1; else EFp2_ECD(&R, &P);
  st_efp2(&R, out); EFp2_clear(&P); EFp2_clear(&R); }

static void op_efp2_add(mpz_t *in, mpz_t *out)
{ EFp2 P, Q, R; EFp2_init(&P); EFp2_init(&Q); EFp2_init(&R);
  ld_efp2(&P, in); ld_efp2(&Q, in+5); EFp2_ECA(&R, &P, &Q);
  st_efp2(&R, out); EFp2_clear(&P); EFp2_clear(&Q); EFp2_clear(&R); }

static void op_efp2_mul(mpz_t *in, mpz_t *out)
{ EFp2 P, R; EFp2_init(&P); EFp2_init(&R); ld_efp2(&P, in);
  EFp2_SCM(&R, &P, in[5]); st_efp2(&R, out); EFp2_clear(&P); EFp2_clear(&R); }

static const Op OPS[] = {
    {"fp_add",1+1,1,op_fp_add},  {"fp_sub",2,1,op_fp_sub},  {"fp_mul",2,1,op_fp_mul},
    {"fp_inv",1,1,op_fp_inv},    {"fp_neg",1,1,op_fp_neg},
    {"fp2_add",4,2,op_fp2_add},  {"fp2_sub",4,2,op_fp2_sub},{"fp2_mul",4,2,op_fp2_mul},
    {"fp2_sqr",2,2,op_fp2_sqr},  {"fp2_inv",2,2,op_fp2_inv},
    {"fp2_mulbasis",2,2,op_fp2_mulbasis},
    {"fp6_add",12,6,op_fp6_add}, {"fp6_sub",12,6,op_fp6_sub},{"fp6_mul",12,6,op_fp6_mul},
    {"fp6_sqr",6,6,op_fp6_sqr},  {"fp6_inv",6,6,op_fp6_inv},
    {"fp6_mulbasis",6,6,op_fp6_mulbasis},
    {"fp12_add",24,12,op_fp12_add},{"fp12_sub",24,12,op_fp12_sub},
    {"fp12_mul",24,12,op_fp12_mul},{"fp12_sqr",12,12,op_fp12_sqr},
    {"fp12_inv",12,12,op_fp12_inv},
    {"efp_dbl",3,3,op_efp_dbl},  {"efp_add",6,3,op_efp_add},{"efp_mul",4,3,op_efp_mul},
    {"efp2_dbl",5,5,op_efp2_dbl},{"efp2_add",10,5,op_efp2_add},{"efp2_mul",6,5,op_efp2_mul},
    {NULL,0,0,NULL}
};

/* ---------------------------------------------------------------------- main */

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: kat_runner <file.vec>\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }

    char *line = NULL; size_t cap = 0; ssize_t len;
    long lineno = 0;
    int primed = 0;
    Toks in, out, want;
    toks_init(&in); toks_init(&out); toks_init(&want);
    mpz_init(PRIME);

    while ((len = getline(&line, &cap, f)) > 0) {
        lineno++;
        if (line[0] == '#' || line[0] == '\n') continue;

        char *save = NULL;
        char *tok = strtok_r(line, " \t\n", &save);
        if (!tok) continue;

        if (strcmp(tok, "prime") == 0) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h || mpz_set_str(PRIME, h, 16) != 0) {
                fprintf(stderr, "bad prime on line %ld\n", lineno); return 2;
            }
            /* Set up just enough library state: the field layer needs only the
             * modulus, and init_precoms derives epsilon and the Frobenius
             * constants from it. Deliberately not calling init_bn_settings(),
             * which would overwrite the prime with its own curve. */
            mpz_init(curve_parameters.prime);
            mpz_set(curve_parameters.prime, PRIME);
            init_precoms(1);
            primed = 1;
            continue;
        }
        if (!primed) { fprintf(stderr, "vectors before 'prime' on line %ld\n", lineno); return 2; }

        const Op *op = NULL;
        for (int i = 0; OPS[i].name; i++)
            if (strcmp(tok, OPS[i].name) == 0) { op = &OPS[i]; break; }
        if (!op) { fprintf(stderr, "unknown op '%s' line %ld\n", tok, lineno); return 2; }
        snprintf(cur_op, sizeof cur_op, "%s", op->name);

        int i = 0;
        for (; i < op->in_n; i++) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h || mpz_set_str(in.v[i], h, 16) != 0) {
                fprintf(stderr, "bad operand %d line %ld\n", i, lineno); return 2;
            }
        }
        char *eq = strtok_r(NULL, " \t\n", &save);
        if (!eq || strcmp(eq, "=") != 0) {
            fprintf(stderr, "missing '=' line %ld\n", lineno); return 2;
        }
        for (i = 0; i < op->out_n; i++) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h || mpz_set_str(want.v[i], h, 16) != 0) {
                fprintf(stderr, "bad result %d line %ld\n", i, lineno); return 2;
            }
        }

        for (i = 0; i < op->out_n; i++) mpz_set_ui(out.v[i], 0);
        op->run(in.v, out.v);

        long before = n_fail;
        for (i = 0; i < op->out_n; i++) cmp1(lineno, i, out.v[i], want.v[i]);
        if (n_fail == before) n_pass++;
    }

    free(line); fclose(f);
    toks_clear(&in); toks_clear(&out); toks_clear(&want);

    printf("%s: %ld passed, %ld failed", argv[1], n_pass, n_fail);
    if (n_noncanon) printf(", %ld non-canonical output value(s)", n_noncanon);
    printf("\n");
    return n_fail ? 1 : 0;
}
