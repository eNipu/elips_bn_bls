/*
 * Does the pairing compute the right value?
 *
 * The reference is test/kat/pairing_<curve>.vec, produced by
 * tools/reference/pairing_ref.py -- an independent optimal ate written from the
 * defining equations, with the final exponentiation done as one exponentiation
 * by (p^12-1)/r rather than by any addition chain.
 *
 * Until issue #17 the reference was the legacy mpz layer, so this test said
 * "the new pairing agrees with the old one". That is the weaker of the two
 * available statements: two implementations can share a misreading of the twist
 * conventions and agree on a wrong answer. Against the oracle they cannot, and
 * BLS12-381 gets a pairing reference for the first time -- the legacy layer
 * never supported it.
 *
 * Exit codes: 0 all passed, 1 a check failed, 2 the vector file is unusable.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "elips/pairing.h"

static int fails, checks;
static long cur_line;

static void ok(int c, const char *what)
{
    checks++;
    if (!c) { fails++; printf("  [FAIL] %s (line %ld)\n", what, cur_line); }
}

static void ld(fp_t r, const mpz_t v)
{
    limb_t l[FP_LIMBS];
    memset(l, 0, sizeof l);
    mpz_export(l, NULL, -1, sizeof(limb_t), 0, 0, v);
    fp_from_limbs(r, l);
}

static void st(mpz_t o, const fp_t a)
{
    limb_t l[FP_LIMBS];
    fp_to_limbs(l, a);
    mpz_import(o, FP_LIMBS, -1, sizeof(limb_t), 0, 0, l);
}

/* The vector file lists the twelve coefficients as d0.c0.a, d0.c0.b, d0.c1.a,
 * ... d1.c2.b, which is the order tools/reference/elips_ref.py's fp12_to_list
 * produces. */
static void fp12_coeffs(const limb_t *out[12], const fp12_t z)
{
    int k = 0;
    for (int d = 0; d < 2; d++)
        for (int c = 0; c < 3; c++)
            for (int b = 0; b < 2; b++)
                out[k++] = z[d][c][b];
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <pairing.vec>\n", argv[0]); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }

    printf("pairing vectors [%s]\n", ELIPS_CURVE_NAME);

    mpz_t prime, tmp, got, want;
    mpz_inits(prime, tmp, got, want, NULL);

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    long records = 0;

    while ((len = getline(&line, &cap, f)) > 0) {
        cur_line++;
        if (line[0] == '#' || line[0] == '\n') continue;

        char *save = NULL;
        char *op = strtok_r(line, " \t\n", &save);
        if (!op) continue;

        if (strcmp(op, "prime") == 0) {
            mpz_set_str(prime, strtok_r(NULL, " \t\n", &save), 16);
            /* The vector file and the compiled curve must be the same one, or
             * every value below is quietly meaningless. */
            mpz_import(tmp, FP_LIMBS, -1, sizeof(limb_t), 0, 0, FP_MODULUS);
            if (mpz_cmp(tmp, prime) != 0) {
                fprintf(stderr, "vector file is for a different curve than %s\n",
                        ELIPS_CURVE_NAME);
                free(line); fclose(f);
                mpz_clears(prime, tmp, got, want, NULL);
                return 2;
            }
            continue;
        }

        if (strcmp(op, "pair") != 0) {
            fprintf(stderr, "unknown op '%s' line %ld\n", op, cur_line);
            free(line); fclose(f);
            mpz_clears(prime, tmp, got, want, NULL);
            return 2;
        }

        /* px py qx0 qx1 qy0 qy1 = c0 .. c11 */
        char *tok[19];
        int missing = 0;
        for (int i = 0; i < 19; i++) {
            tok[i] = strtok_r(NULL, " \t\n", &save);
            if (!tok[i]) missing = 1;
        }
        if (missing || strcmp(tok[6], "=") != 0) {
            fprintf(stderr, "malformed record at line %ld\n", cur_line);
            free(line); fclose(f);
            mpz_clears(prime, tmp, got, want, NULL);
            return 2;
        }

        fp_t px, py;
        fp2_t qx, qy;
        mpz_set_str(tmp, tok[0], 16); ld(px, tmp);
        mpz_set_str(tmp, tok[1], 16); ld(py, tmp);
        mpz_set_str(tmp, tok[2], 16); ld(qx[0], tmp);
        mpz_set_str(tmp, tok[3], 16); ld(qx[1], tmp);
        mpz_set_str(tmp, tok[4], 16); ld(qy[0], tmp);
        mpz_set_str(tmp, tok[5], 16); ld(qy[1], tmp);

        /* The inputs must be what the file says they are before the output is
         * worth comparing. */
        ep_t P;
        ep2_t Q;
        ep_from_affine(&P, px, py);
        ep2_from_affine(&Q, qx, qy);
        ok(ep_in_subgroup(&P),  "P from the vector is in G1");
        ok(ep2_in_subgroup(&Q), "Q from the vector is in G2");

        fp12_t raw, exact;
        pairing_miller(raw, qx, qy, px, py);
        pairing_final_exp_plain(exact, raw);

        const limb_t *c[12];
        fp12_coeffs(c, exact);
        int diff = 0;
        for (int i = 0; i < 12; i++) {
            st(got, c[i]);
            mpz_set_str(want, tok[7 + i], 16);
            if (mpz_cmp(got, want) != 0) {
                if (diff == 0)
                    gmp_printf("  coord %2d  got  %Zx\n            want %Zx\n",
                               i, got, want);
                diff++;
            }
        }
        ok(diff == 0, "exact pairing matches the reference value");

        /* Properties the vector cannot fake: a reference that was itself wrong
         * would still have to satisfy these. */
        fp12_t one, chk;
        fp12_set_one(one);
        ok(!fp12_eq(exact, one), "e(P,Q) != 1");
        fp12_exp(chk, exact, ELIPS_ORDER, ELIPS_ORDER_BITS);
        ok(fp12_eq(chk, one), "e(P,Q)^r == 1");

        /* And the relationship between the two final exponentiations, which is
         * internal to this library and has no reference outside it. */
        fp12_t fast;
        pairing_final_exp_fast(fast, raw);
#ifdef ELIPS_FAMILY_BLS12
        /* The fast chain computes e^3 by design: 3*lambda has a short
         * evaluation and lambda alone does not. See issue #16. */
        fp12_t cube;
        fp12_sqr(cube, exact);
        fp12_mul(cube, cube, exact);
        ok(fp12_eq(fast, cube), "fast final exponentiation == exact^3");
#else
        /* BN's hard part decomposes exactly, so no stray factor. */
        ok(fp12_eq(fast, exact), "fast final exponentiation == exact e");
#endif

        /* elips_pairing must agree with the pieces it is built from. */
        fp12_t whole;
        ok(elips_pairing(whole, &P, &Q) == 1, "elips_pairing accepts the inputs");
        ok(fp12_eq(whole, fast), "elips_pairing == miller + fast final exp");

        records++;
    }

    free(line);
    fclose(f);
    mpz_clears(prime, tmp, got, want, NULL);

    if (records == 0) {
        fprintf(stderr, "no pairing records in %s\n", argv[1]);
        return 2;
    }

    printf("  %ld pairings, %d checks, %d failed\n", records, checks, fails);
    return fails ? 1 : 0;
}
