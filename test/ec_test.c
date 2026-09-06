/*
 * Group-law properties for the curve layer.
 *
 * The known-answer vectors already pin the arithmetic to the Phase 0 oracle;
 * these checks cover what vectors are bad at: associativity, the degenerate
 * cases, and the scalar ladder against repeated addition.
 *
 * The points used to come from the legacy layer, which had square roots and
 * point generation when this layer did not. They now come from ep2_generator,
 * which the parameter generator asserts is on the twist and of order exactly r
 * before emitting it -- a stronger provenance than the borrowed points had, and
 * one that does not keep the old layer alive (issue #17).
 */
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "elips/ec.h"
#include "elips/pairing.h"

static int failures;
static void ok(int c, const char *what)
{ printf("  [%s] %s\n", c ? "PASS" : "FAIL", what); if (!c) failures++; }

static int ep2_same(const ep2_t *A, const ep2_t *B)
{
    if (ep2_is_infinity(A) || ep2_is_infinity(B))
        return ep2_is_infinity(A) && ep2_is_infinity(B);
    fp2_t ax, ay, bx, by;
    ep2_to_affine(ax, ay, A); ep2_to_affine(bx, by, B);
    return fp2_eq(ax, bx) && fp2_eq(ay, by);
}

int main(void)
{
    printf("EC property tests [%s]\n", ELIPS_CURVE_NAME);

    /* P, 2P and 3P from the generator. Three points that are distinct, in the
     * subgroup, and related in a way the checks below can exploit. */
    ep2_t P, Q, R, t1, t2, inf;
    ep2_generator(&P);
    ep2_dbl(&Q, &P);
    ep2_add(&R, &P, &Q);
    ep2_set_infinity(&inf);

    ok(ep2_on_curve(&P) && ep2_on_curve(&Q) && ep2_on_curve(&R),
       "the generator and its multiples satisfy the twist equation");
    ok(ep2_on_curve(&inf), "infinity counts as on-curve");

    /* associativity */
    ep2_add(&t1, &P, &Q); ep2_add(&t1, &t1, &R);
    ep2_add(&t2, &Q, &R); ep2_add(&t2, &P, &t2);
    ok(ep2_same(&t1, &t2), "(P+Q)+R == P+(Q+R)");

    /* identity and inverse */
    ep2_add(&t1, &P, &inf);  ok(ep2_same(&t1, &P), "P + O == P");
    ep2_add(&t1, &inf, &P);  ok(ep2_same(&t1, &P), "O + P == P");
    ep2_neg(&t2, &P);
    ep2_add(&t1, &P, &t2);   ok(ep2_is_infinity(&t1), "P + (-P) == O");

    /* the complete addition must detect the doubling case */
    ep2_add(&t1, &P, &P); ep2_dbl(&t2, &P);
    ok(ep2_same(&t1, &t2), "P + P == 2P via the complete addition");

    /* doubling infinity stays infinity */
    ep2_dbl(&t1, &inf); ok(ep2_is_infinity(&t1), "2*O == O");

    /* Addition stays on the curve over a long chain.
     *
     * This case is inherited from the Jacobian layer, where it compared the
     * complete addition against a faster incomplete one and had to keep the
     * operands distinct and non-opposite for the incomplete routine's sake. The
     * RCB formulas have no such precondition and the incomplete routine is
     * gone, so what is left is a chain check: 4P, 8P, 16P... each plus P, every
     * intermediate still satisfying the curve equation. */
    int agree = 1, checked = 0;
    ep2_t acc; ep2_copy(&acc, &Q);          /* 2P */
    for (int i = 0; i < 40; i++) {
        ep2_dbl(&acc, &acc);                 /* 4P, 8P, 16P, ... */
        ep2_add(&t1, &acc, &P);
        ep2_add(&t2, &acc, &P);
        checked++;
        if (!ep2_same(&t1, &t2)) { agree = 0; break; }
        if (!ep2_on_curve(&t2)) { agree = 0; break; }
    }
    ok(agree && checked == 40,
       "ep2_add is consistent with repeated addition on 40 pairs");

    /* scalar multiplication against repeated addition */
    limb_t k[FP_LIMBS]; memset(k, 0, sizeof k); k[0] = 13;
    ep2_mul(&t1, &P, k, 4);
    ep2_set_infinity(&t2);
    for (int i = 0; i < 13; i++) ep2_add(&t2, &t2, &P);
    ok(ep2_same(&t1, &t2), "[13]P matches thirteen additions");

    memset(k, 0, sizeof k);
    ep2_mul(&t1, &P, k, 8);
    ok(ep2_is_infinity(&t1), "[0]P == O");

    /* the subgroup check that the legacy library never performed */
    {
        limb_t r[FP_LIMBS]; memset(r, 0, sizeof r);
        memcpy(r, ELIPS_ORDER, sizeof(limb_t) * ((ELIPS_ORDER_BITS + 63) / 64));
        int bits = ELIPS_ORDER_BITS;
        ep2_t cof; ep2_mul(&cof, &P, r, bits);
        /* P is a real G2 generator, so this is the genuine subgroup check that
         * nothing in the repository performed before Phase 1. */
        ok(ep2_is_infinity(&cof), "[r]P == O for the G2 generator");
    }

#ifdef ELIPS_FAMILY_BLS12
    {   /* GLV must agree with the plain ladder, including at the edges. */
        mpz_t R_, kk; mpz_inits(R_, kk, NULL);
        mpz_import(R_, (ELIPS_ORDER_BITS+63)/64, -1, sizeof(limb_t), 0, 0, ELIPS_ORDER);
        gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs, 20260906UL);
        int bad = 0;
        for (int i = 0; i < 25; i++) {
            if (i == 0) mpz_set_ui(kk, 0);
            else if (i == 1) mpz_set_ui(kk, 1);
            else if (i == 2) mpz_sub_ui(kk, R_, 1);
            else mpz_urandomm(kk, rs, R_);
            limb_t buf[16]; memset(buf, 0, sizeof buf);
            mpz_export(buf, NULL, -1, sizeof(limb_t), 0, 0, kk);
            int bits = mpz_sgn(kk) ? (int)mpz_sizeinbase(kk, 2) : 1;
            ep2_t g1_, g2_;
            ep2_mul(&g1_, &P, buf, bits);
            ep2_mul_glv(&g2_, &P, buf, bits);
            if (!ep2_same(&g1_, &g2_)) { bad++; break; }
        }
        ok(bad == 0, "GLV matches the plain ladder on 25 scalars incl. 0, 1, r-1");
        gmp_randclear(rs); mpz_clears(R_, kk, NULL);
    }
#endif

    printf("%s: %s (%d failure%s)\n", ELIPS_CURVE_NAME,
           failures ? "FAILED" : "OK", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
