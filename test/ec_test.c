/*
 * Group-law properties for the new Jacobian layer.
 *
 * The known-answer vectors already pin the arithmetic to the Phase 0 oracle;
 * these checks cover the things vectors are bad at: associativity, the
 * degenerate cases, and the agreement between the complete addition and the
 * faster incomplete one the Miller loop will use.
 *
 * Valid curve points come from the legacy layer, which has square roots and
 * point generation. Borrowing them keeps this test self-contained without
 * adding API surface to the new layer before anything needs it.
 */
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "elips/ec.h"
#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bn_efp.h>
#include <ELiPS_bn_bls/bn_efp2.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>
#include <ELiPS_bn_bls/curve_settings.h>

static int failures;
static void ok(int c, const char *what)
{ printf("  [%s] %s\n", c ? "PASS" : "FAIL", what); if (!c) failures++; }

static void ld(fp_t r, const mpz_t v)
{ limb_t l[FP_LIMBS]; memset(l,0,sizeof l);
  mpz_export(l,NULL,-1,sizeof(limb_t),0,0,v); fp_from_limbs(r,l); }

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
    bls12_inits();
    printf("EC property tests [%s]\n", ELIPS_CURVE_NAME);

    /* A genuine G2 point, taken through the legacy untwisting map.
     *
     * Note: EFp2_rational_point cannot be used here. It builds points on
     * y^2 = x^3 - b over Fp2, which is NOT the sextic twist y^2 = x^3 + b*xi
     * that carries G2 -- measured residual p-4 rather than 4+4u. The legacy
     * tests never noticed because the group law does not reference b. */
    EFp12 G2; EFp12_init(&G2); bls12_generate_G2_point(&G2);
    EFp2 L1, L2, L3;
    EFp2_init(&L1); EFp2_init(&L2); EFp2_init(&L3);
    bls12_EFp12_to_EFp2(&L1, &G2);
    EFp2_ECD(&L2, &L1);
    EFp2_ECA(&L3, &L1, &L2);

    ep2_t P, Q, R, t1, t2, inf;
    fp2_t x, y;
    #define CONV(dst, src) do { \
        ld(x[0], (src).x.x0.x0); ld(x[1], (src).x.x1.x0); \
        ld(y[0], (src).y.x0.x0); ld(y[1], (src).y.x1.x0); \
        ep2_from_affine(&(dst), x, y); } while (0)
    CONV(P, L1); CONV(Q, L2); CONV(R, L3);
    ep2_set_infinity(&inf);

    ok(ep2_on_curve(&P) && ep2_on_curve(&Q) && ep2_on_curve(&R),
       "borrowed points satisfy the twist equation");
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

    /* The incomplete addition must agree wherever its precondition holds.
     *
     * Choosing the inputs matters: accumulating 2P, 4P, 8P... and adding P
     * keeps the operands distinct and non-opposite throughout. An earlier
     * version of this test accumulated from P and added 2P, which hit
     * acc == Q on the very first iteration -- the doubling case, exactly what
     * ep2_add_generic is documented not to handle. The test was wrong, not the
     * code, but it is a fair demonstration that the precondition is real. */
    int agree = 1, checked = 0;
    ep2_t acc; ep2_copy(&acc, &Q);          /* 2P */
    for (int i = 0; i < 40; i++) {
        ep2_dbl(&acc, &acc);                 /* 4P, 8P, 16P, ... */
        ep2_add(&t1, &acc, &P);
        ep2_add_generic(&t2, &acc, &P);
        checked++;
        if (!ep2_same(&t1, &t2)) { agree = 0; break; }
        if (!ep2_on_curve(&t2)) { agree = 0; break; }
    }
    ok(agree && checked == 40,
       "ep2_add_generic matches ep2_add on 40 admissible pairs");

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
        mpz_export(r, NULL, -1, sizeof(limb_t), 0, 0, curve_parameters.order);
        int bits = (int)mpz_sizeinbase(curve_parameters.order, 2);
        ep2_t cof; ep2_mul(&cof, &P, r, bits);
        /* P is a real G2 generator, so this is the genuine subgroup check that
         * nothing in the repository performed before Phase 1. */
        ok(ep2_is_infinity(&cof), "[r]P == O for the G2 generator");
    }

    printf("%s: %s (%d failure%s)\n", ELIPS_CURVE_NAME,
           failures ? "FAILED" : "OK", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
