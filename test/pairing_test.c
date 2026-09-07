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

    /* Kept for the multi-pairing checks at the end. */
#define MAXREC 16
    ep_t   allP[MAXREC];
    ep2_t  allQ[MAXREC];
    fp12_t allExact[MAXREC];

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

        /* Fixed-argument precomputation. The claim is bit-identity with the
         * direct Miller loop, not agreement after the final exponentiation --
         * the exponentiation would hide any per-line Fp2 scaling, and hiding
         * the difference is exactly what would let a wrong table pass. */
        ep2_prec_t pc;
        ok(ep2_precompute(&pc, &Q) == 1, "ep2_precompute accepts a G2 point");
        fp12_t fprec;
        pairing_miller_prec(fprec, &pc, px, py);
        ok(fp12_eq(fprec, raw), "precomputed Miller == direct Miller, bit for bit");

        fp12_t eprec;
        ok(elips_pairing_prec(eprec, &P, &pc) == 1, "elips_pairing_prec accepts P");
        ok(fp12_eq(eprec, whole), "elips_pairing_prec == elips_pairing");

        if (records < MAXREC) {
            ep_copy(&allP[records], &P);
            ep2_copy(&allQ[records], &Q);
            fp12_copy(allExact[records], exact);
        }
        records++;
    }

    free(line);
    fclose(f);
    mpz_clears(prime, tmp, got, want, NULL);

    if (records == 0) {
        fprintf(stderr, "no pairing records in %s\n", argv[1]);
        return 2;
    }

    /* ---------------------------------------------------- multi-pairing ---
     *
     * Two independent statements, because each can pass while the other fails.
     *
     * Against the oracle: the product of the recorded exact values is a number
     * this library never computed, so agreeing with it says the shared loop
     * really evaluated every pair. On BLS12 the fast chain contributes a cube
     * to each factor, so the product picks up one cube overall.
     *
     * Against itself: raising to a fixed exponent is a homomorphism, so the
     * multi-pairing must equal the product of the individual elips_pairing
     * results bit for bit. That catches a wrong accumulator or a dropped term
     * even if the oracle comparison were somehow satisfied.
     */
    long nrec = records < MAXREC ? records : MAXREC;
    cur_line = 0;

    if (nrec >= 2) {
        fp12_t multi, prod_exact, prod_single, one;
        fp12_set_one(one);

        ok(elips_pairing_multi(multi, allP, allQ, (size_t)nrec) == 1,
           "elips_pairing_multi accepts the vector points");

        fp12_set_one(prod_exact);
        for (long i = 0; i < nrec; i++) fp12_mul(prod_exact, prod_exact, allExact[i]);
#ifdef ELIPS_FAMILY_BLS12
        {   /* the fast chain cubes each factor, so the product picks up one cube */
            fp12_t sq;
            fp12_sqr(sq, prod_exact);
            fp12_mul(prod_exact, sq, prod_exact);
        }
#endif
        ok(fp12_eq(multi, prod_exact),
           "multi-pairing == product of the reference values");

        fp12_set_one(prod_single);
        for (long i = 0; i < nrec; i++) {
            fp12_t e1;
            if (elips_pairing(e1, &allP[i], &allQ[i]) != 1) { fails++; checks++; continue; }
            fp12_mul(prod_single, prod_single, e1);
        }
        ok(fp12_eq(multi, prod_single),
           "multi-pairing == product of individual elips_pairing");

        /* Cross the internal chunk boundary. Ten terms cycled from the
         * available records, so the run is longer than one chunk however the
         * chunk size is tuned. */
        ep_t  longP[10];
        ep2_t longQ[10];
        for (int i = 0; i < 10; i++) {
            ep_copy(&longP[i],  &allP[i % nrec]);
            ep2_copy(&longQ[i], &allQ[i % nrec]);
        }
        fp12_t mlong, plong;
        ok(elips_pairing_multi(mlong, longP, longQ, 10) == 1,
           "multi-pairing accepts a run longer than one chunk");
        fp12_set_one(plong);
        for (int i = 0; i < 10; i++) {
            fp12_t e1;
            if (elips_pairing(e1, &longP[i], &longQ[i]) != 1) { fails++; checks++; continue; }
            fp12_mul(plong, plong, e1);
        }
        ok(fp12_eq(mlong, plong),
           "multi-pairing across a chunk boundary == product of singles");

        /* The empty product is one, and it is not an error. */
        fp12_t empty;
        ok(elips_pairing_multi(empty, allP, allQ, 0) == 1,
           "multi-pairing of nothing succeeds");
        ok(fp12_eq(empty, one), "multi-pairing of nothing is one");

        /* Precomputation and multi-pairing together, which is the shape a
         * verification actually has: fixed G2 arguments, one product. */
        ep2_prec_t pcs[MAXREC];
        int pc_ok = 1;
        for (long i = 0; i < nrec; i++)
            if (ep2_precompute(&pcs[i], &allQ[i]) != 1) pc_ok = 0;
        ok(pc_ok, "every vector G2 point precomputes");

        fp12_t mprec;
        ok(elips_pairing_multi_prec(mprec, allP, pcs, (size_t)nrec) == 1,
           "elips_pairing_multi_prec accepts the vector points");
        ok(fp12_eq(mprec, multi),
           "multi-pairing with precomputation == multi-pairing without");

        /* The count the table is sized by must match the loop that fills it.
         * If these disagree, ep2_precompute writes past the end of the table,
         * so it is checked here rather than discovered there. */
        {
            int lines = 0;
            for (int i = ELIPS_LOOP_TOP - 1; i >= 0; i--)
                lines += (ELIPS_LOOP[i] != 0) ? 2 : 1;
#ifdef ELIPS_FAMILY_BN
            lines += 2;
#endif
            ok(lines == ELIPS_MILLER_LINES,
               "ELIPS_MILLER_LINES matches the loop that fills the table");
        }

        /* A precomputed table must refuse a Q outside G2, since nothing
         * downstream can check it any more. */
        {
            ep2_prec_t bad_pc;
            ep2_t inf2;
            ep2_set_infinity(&inf2);
            ok(ep2_precompute(&bad_pc, &inf2) == 0,
               "ep2_precompute rejects the identity");
        }

        /* And a bad input must be refused, with the result left at one rather
         * than a product over whichever pairs were checked first. */
        ep_t  badP[2];
        ep2_t badQ[2];
        ep_copy(&badP[0], &allP[0]);  ep2_copy(&badQ[0], &allQ[0]);
        ep_set_infinity(&badP[1]);    ep2_copy(&badQ[1], &allQ[0]);
        fp12_t bad;
        ok(elips_pairing_multi(bad, badP, badQ, 2) == 0,
           "multi-pairing rejects an identity input");
        ok(fp12_eq(bad, one), "a rejected multi-pairing leaves the result at one");
    }

    {   /* Constant-time G_T exponentiation against the reference ladder.
         *
         * fp12_exp_gt splits the exponent with the Frobenius and squares in
         * the cyclotomic subgroup; fp12_exp is a plain square-and-multiply
         * over the whole exponent. They must agree on every scalar, and a
         * decomposition that is wrong for one scalar in a million would pass
         * a casual test, so this includes 0, 1, r-1 and values straddling the
         * digit base as well as random ones.
         *
         * Both are exercised on a real pairing output, which is the only kind
         * of element fp12_exp_gt is defined on. */
        fp12_t base;
        ep_t Pg; ep2_t Qg;
        ep_generator(&Pg); ep2_generator(&Qg);
        ok(elips_pairing(base, &Pg, &Qg) == 1, "pairing of the generators for the G_T test");

        mpz_t R_, kk;
        mpz_inits(R_, kk, NULL);
        mpz_import(R_, (ELIPS_ORDER_BITS + 63) / 64, -1, sizeof(limb_t), 0, 0, ELIPS_ORDER);
        gmp_randstate_t rs;
        gmp_randinit_default(rs);
        gmp_randseed_ui(rs, 20260907UL);

        int bad_gt = 0;
        for (int i = 0; i < 40; i++) {
            if (i == 0)      mpz_set_ui(kk, 0);
            else if (i == 1) mpz_set_ui(kk, 1);
            else if (i == 2) mpz_sub_ui(kk, R_, 1);
            else if (i == 3) mpz_set_ui(kk, 2);
            else if (i < 12) {
                int nb = ELIPS_ORDER_BITS / 2 + (i - 4) - 3;
                mpz_ui_pow_ui(kk, 2, (unsigned)(nb > 1 ? nb : 1));
                if (i % 3 == 1) mpz_sub_ui(kk, kk, 1);
                if (i % 3 == 2) mpz_add_ui(kk, kk, 1);
                mpz_mod(kk, kk, R_);
            } else {
                mpz_urandomm(kk, rs, R_);
            }
            limb_t kb[16];
            memset(kb, 0, sizeof kb);
            mpz_export(kb, NULL, -1, sizeof(limb_t), 0, 0, kk);
            int kbits = mpz_sgn(kk) ? (int)mpz_sizeinbase(kk, 2) : 1;

            fp12_t slow, fast_gt;
            fp12_exp(slow, base, kb, kbits);
            fp12_exp_gt(fast_gt, base, kb, kbits);
            if (!fp12_eq(slow, fast_gt)) bad_gt++;
        }
        ok(bad_gt == 0, "fp12_exp_gt matches fp12_exp on 40 scalars");

        gmp_randclear(rs);
        mpz_clears(R_, kk, NULL);
    }

    printf("  %ld pairings, %d checks, %d failed\n", records, checks, fails);
    return fails ? 1 : 0;
}
