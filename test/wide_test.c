/*
 * elips_mod_wide against GMP.
 *
 * This test is the reason GMP is still a dependency of the test suite after it
 * stopped being one of the library. elips_mod_wide replaced mpn_sec_div_r, and
 * a replacement checked only against itself is not checked at all. GMP is an
 * independent implementation, so agreement across a large sample is evidence.
 *
 * Boundary cases come first and are not decoration. Modular reduction gets
 * wrong exactly the values around the modulus: m-1, m, m+1, 2m-1, 2m, 2m+1,
 * exact multiples, zero, and the largest value the width can hold. A random
 * sweep alone would hit none of them with meaningful probability.
 *
 * Both moduli are covered, at the widths the library actually calls with: the
 * field prime for hash_to_field and fp_rand, and the group order for
 * elips_random_scalar.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "elips/fp.h"
#include "elips/random.h"
#include "arith/wide.h"

#define MAXW 64

static long checks, failures;

static void from_mpz(limb_t *o, int n, const mpz_t v)
{
    memset(o, 0, (size_t)n * sizeof(limb_t));
    mpz_export(o, NULL, -1, sizeof(limb_t), 0, 0, v);
}

static void one(const limb_t *w, int nn, const limb_t *m, int dn, const char *what)
{
    mpz_t W, M, R;
    mpz_inits(W, M, R, NULL);
    mpz_import(W, (size_t)nn, -1, sizeof(limb_t), 0, 0, w);
    mpz_import(M, (size_t)dn, -1, sizeof(limb_t), 0, 0, m);
    mpz_mod(R, W, M);

    limb_t got[MAXW];
    memcpy(got, w, (size_t)nn * sizeof(limb_t));
    elips_mod_wide(got, nn, m, dn);

    limb_t want[MAXW];
    memset(want, 0, sizeof want);
    from_mpz(want, dn, R);

    checks++;
    if (memcmp(got, want, (size_t)dn * sizeof(limb_t)) != 0) {
        failures++;
        if (failures <= 5)
            gmp_fprintf(stderr, "MISMATCH %s\n  w = %Zx\n  m = %Zx\n  want %Zx\n",
                        what, W, M, R);
    }
    /* The contract says the limbs above the remainder are cleared. A caller
     * reading nn limbs, as fp_from_limbs does not but a future one might,
     * would otherwise see leftover dividend. */
    checks++;
    for (int i = dn; i < nn; i++) {
        if (got[i] != 0) {
            failures++;
            fprintf(stderr, "MISMATCH %s: limb %d above the remainder is not zero\n", what, i);
            break;
        }
    }
    mpz_clears(W, M, R, NULL);
}

static void sweep(const limb_t *m, int dn, int nn, long iters, const char *label)
{
    mpz_t M, w, t;
    mpz_inits(M, w, t, NULL);
    mpz_import(M, (size_t)dn, -1, sizeof(limb_t), 0, 0, m);
    limb_t buf[MAXW];

    static const char *names[] = { "0", "1", "m-1", "m", "m+1", "2m-1", "2m", "2m+1" };
    for (int i = 0; i < 8; i++) {
        switch (i) {
        case 0: mpz_set_ui(w, 0); break;
        case 1: mpz_set_ui(w, 1); break;
        case 2: mpz_sub_ui(w, M, 1); break;
        case 3: mpz_set(w, M); break;
        case 4: mpz_add_ui(w, M, 1); break;
        case 5: mpz_mul_ui(t, M, 2); mpz_sub_ui(w, t, 1); break;
        case 6: mpz_mul_ui(w, M, 2); break;
        case 7: mpz_mul_ui(t, M, 2); mpz_add_ui(w, t, 1); break;
        }
        from_mpz(buf, nn, w);
        one(buf, nn, m, dn, names[i]);
    }

    /* The widest value the buffer can hold, which is where a missing carry in
     * the shift shows up. */
    mpz_ui_pow_ui(w, 2, (unsigned long)nn * 64);
    mpz_sub_ui(w, w, 1);
    from_mpz(buf, nn, w);
    one(buf, nn, m, dn, "2^(64*nn) - 1");

    /* Exact multiples must give exactly zero. */
    for (unsigned long k = 3; k < 48; k++) {
        mpz_mul_ui(w, M, k);
        if (mpz_sizeinbase(w, 2) > (size_t)nn * 64) break;
        from_mpz(buf, nn, w);
        one(buf, nn, m, dn, "k*m");
    }

    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, 20260912);          /* fixed: a failure must reproduce */
    mpz_t hi; mpz_init(hi);
    mpz_ui_pow_ui(hi, 2, (unsigned long)nn * 64);
    for (long i = 0; i < iters; i++) {
        mpz_urandomm(w, st, hi);
        from_mpz(buf, nn, w);
        one(buf, nn, m, dn, "random");
    }
    mpz_clear(hi);
    gmp_randclear(st);
    mpz_clears(M, w, t, NULL);
    printf("  %-42s ok\n", label);
}

int main(int argc, char **argv)
{
    long iters = (argc > 1) ? atol(argv[1]) : 4000;
    printf("%s, elips_mod_wide against GMP\n", ELIPS_CURVE_NAME);

    sweep(FP_MODULUS, (int)FP_LIMBS, (int)FP_LIMBS + 2, iters,
          "wide mod p, the fp_rand width");
    sweep(FP_MODULUS, (int)FP_LIMBS, 2 * (int)FP_LIMBS, iters,
          "wide mod p, double width");
    sweep(ELIPS_ORDER, (int)ELIPS_ORDER_LIMBS, (int)ELIPS_ORDER_LIMBS + 2, iters,
          "wide mod r, the random-scalar width");

    printf("  %ld checks, %ld failures\n", checks, failures);
    return failures ? 1 : 0;
}
