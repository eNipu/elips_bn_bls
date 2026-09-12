/*
 * Differential test for the new Montgomery Fp layer (Phase 3, issue #4).
 *
 * The new arithmetic is checked against GMP's mpz on random inputs. GMP is an
 * independent implementation, so agreement across a large sample is real
 * evidence rather than a restatement of our own code.
 *
 * Deliberately includes the boundary values that a Montgomery implementation
 * gets wrong: 0, 1, p-1, and values that force the conditional subtraction.
 *
 * Since Phase 6 there can be two fp_mul implementations in the binary, the
 * portable C and an assembly backend the CPU may or may not have. Both are run
 * against GMP here, side by side, on every check. Which one fp_mul itself
 * dispatched to is printed at the end, so a green run states the path it
 * exercised instead of leaving it to be inferred from the build flags.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "elips/fp.h"

static mpz_t P;
static long failures, checks;

static void limbs_of(limb_t *out, const mpz_t v)
{
    memset(out, 0, FP_BYTES);
    mpz_export(out, NULL, -1, sizeof(limb_t), 0, 0, v);
}
static void mpz_of(mpz_t out, const limb_t *in)
{
    mpz_import(out, FP_LIMBS, -1, sizeof(limb_t), 0, 0, in);
}

/* Round-trip through Montgomery form and compare against the mpz answer. */
static void check(const char *what, const mpz_t want, const fp_t got)
{
    limb_t plain[FP_LIMBS];
    fp_to_limbs(plain, got);
    mpz_t g; mpz_init(g); mpz_of(g, plain);
    mpz_t w; mpz_init(w); mpz_mod(w, want, P);
    checks++;
    if (mpz_cmp(g, w) != 0) {
        failures++;
        if (failures <= 6)
            gmp_fprintf(stderr, "MISMATCH %s\n  expected %Zx\n  actual   %Zx\n", what, w, g);
    }
    mpz_clear(g); mpz_clear(w);
}

int main(int argc, char **argv)
{
    long iters = (argc > 1) ? atol(argv[1]) : 20000;

    mpz_init(P);
    { limb_t m[FP_LIMBS]; memcpy(m, FP_MODULUS, sizeof m); mpz_of(P, m); }
    gmp_printf("%s: p = %Zd\n  %d limbs, %ld iterations\n",
               ELIPS_CURVE_NAME, P, (int)FP_LIMBS, iters);

    /* p must actually be the prime the reference implementation uses. */
    if (mpz_probab_prime_p(P, 25) == 0) {
        fprintf(stderr, "FATAL: generated modulus is not prime\n");
        return 2;
    }

    gmp_randstate_t st; gmp_randinit_default(st); gmp_randseed_ui(st, 20260906UL);

    mpz_t a, b, w; mpz_inits(a, b, w, NULL);
    fp_t fa, fb, fr;
    limb_t la[FP_LIMBS], lb[FP_LIMBS];

    /* Boundary cases first: these are where Montgomery code actually breaks. */
    const int NEDGE = 5;
    for (int i = 0; i < NEDGE; i++) {
        for (int j = 0; j < NEDGE; j++) {
            switch (i) {
                case 0: mpz_set_ui(a, 0); break;
                case 1: mpz_set_ui(a, 1); break;
                case 2: mpz_sub_ui(a, P, 1); break;
                case 3: mpz_sub_ui(a, P, 2); break;
                default: mpz_tdiv_q_ui(a, P, 2); break;
            }
            switch (j) {
                case 0: mpz_set_ui(b, 0); break;
                case 1: mpz_set_ui(b, 1); break;
                case 2: mpz_sub_ui(b, P, 1); break;
                case 3: mpz_sub_ui(b, P, 2); break;
                default: mpz_tdiv_q_ui(b, P, 2); break;
            }
            limbs_of(la, a); limbs_of(lb, b);
            fp_from_limbs(fa, la); fp_from_limbs(fb, lb);

            mpz_add(w, a, b);  fp_add(fr, fa, fb); check("edge add", w, fr);
            mpz_sub(w, a, b);  fp_sub(fr, fa, fb); check("edge sub", w, fr);
            mpz_mul(w, a, b);  fp_mul(fr, fa, fb); check("edge mul", w, fr);
            mpz_mul(w, a, b);  fp_mul_portable(fr, fa, fb);
                                                  check("edge mul portable", w, fr);
            mpz_mul(w, a, a);  fp_sqr(fr, fa);     check("edge sqr", w, fr);
            mpz_neg(w, a);     fp_neg(fr, fa);     check("edge neg", w, fr);
        }
    }

    for (long k = 0; k < iters; k++) {
        mpz_urandomm(a, st, P);
        mpz_urandomm(b, st, P);
        limbs_of(la, a); limbs_of(lb, b);
        fp_from_limbs(fa, la); fp_from_limbs(fb, lb);

        mpz_add(w, a, b); fp_add(fr, fa, fb); check("add", w, fr);
        mpz_sub(w, a, b); fp_sub(fr, fa, fb); check("sub", w, fr);
        mpz_mul(w, a, b); fp_mul(fr, fa, fb); check("mul", w, fr);
        mpz_mul(w, a, b); fp_mul_portable(fr, fa, fb); check("mul portable", w, fr);
        mpz_mul(w, a, a); fp_sqr(fr, fa);     check("sqr", w, fr);
        mpz_neg(w, a);    fp_neg(fr, fa);     check("neg", w, fr);

        /* Montgomery round-trip must be the identity. */
        limb_t back[FP_LIMBS];
        fp_to_limbs(back, fa);
        if (memcmp(back, la, FP_BYTES) != 0) { failures++; checks++;
            fprintf(stderr, "MISMATCH montgomery round-trip\n"); }
        else checks++;

        /* Inversion, on a sample so the test stays quick. */
        if ((k % 64) == 0 && mpz_sgn(a) != 0) {
            mpz_invert(w, a, P);
            fp_inv(fr, fa);
            check("inv", w, fr);
            fp_t prod; fp_mul(prod, fr, fa);
            fp_t one;  fp_set_one(one);
            checks++;
            if (!fp_eq(prod, one)) { failures++; fprintf(stderr, "MISMATCH a*a^-1 != 1\n"); }
        }

        /* Predicates must agree with the obvious answer. */
        checks += 2;
        if (fp_is_zero(fa) != (mpz_sgn(a) == 0)) { failures++; fprintf(stderr, "MISMATCH is_zero\n"); }
        if (fp_eq(fa, fb) != (mpz_cmp(a, b) == 0)) { failures++; fprintf(stderr, "MISMATCH eq\n"); }
    }

    /* Zero is excluded from the loop above because it has no inverse, which
     * makes it exactly the case a constant-time inversion is most likely to get
     * wrong: the natural implementation early-returns on it, and an early
     * return is a branch on the operand. fp_inv is defined to map 0 to 0 and to
     * do it without branching, so pin the value here. */
    {
        fp_t z, iz;
        fp_set_zero(z);
        fp_inv(iz, z);
        checks++;
        if (!fp_is_zero(iz)) { failures++; fprintf(stderr, "MISMATCH fp_inv(0) != 0\n"); }
    }

    /* The two backends must agree bit for bit, including when the result
     * aliases an input. Vectors never alias, and the assembly writes its
     * output only after its last read of a and b, which is a property no
     * value-based test would catch if it were broken. */
    for (long k = 0; k < iters; k++) {
        mpz_urandomm(a, st, P);
        mpz_urandomm(b, st, P);
        limbs_of(la, a); limbs_of(lb, b);
        fp_from_limbs(fa, la); fp_from_limbs(fb, lb);

        fp_t want, u, v;
        fp_mul_portable(want, fa, fb);

        fp_mul(fr, fa, fb);
        checks++;
        if (!fp_eq(fr, want)) { failures++; fprintf(stderr, "MISMATCH backends disagree\n"); }

        fp_copy(u, fa); fp_mul(u, u, fb);
        checks++;
        if (!fp_eq(u, want)) { failures++; fprintf(stderr, "MISMATCH alias r==a\n"); }

        fp_copy(v, fb); fp_mul(v, fa, v);
        checks++;
        if (!fp_eq(v, want)) { failures++; fprintf(stderr, "MISMATCH alias r==b\n"); }

        fp_mul_portable(want, fa, fa);
        fp_copy(u, fa); fp_mul(u, u, u);
        checks++;
        if (!fp_eq(u, want)) { failures++; fprintf(stderr, "MISMATCH alias r==a==b\n"); }
    }

    mpz_clears(a, b, w, NULL); gmp_randclear(st); mpz_clear(P);
    printf("  %ld checks, %ld failures, fp_mul backend: %s\n",
           checks, failures, fp_mul_backend_name());
    return failures ? 1 : 0;
}
