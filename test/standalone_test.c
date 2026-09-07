/*
 * The new layer on its own: no legacy library linked at all.
 *
 * Until this existed every test borrowed points from the old implementation,
 * which is the last thing keeping that layer alive. Everything here comes from
 * generated constants.
 */
#include <stdio.h>
#include <string.h>
#include "elips/pairing.h"

static int fails;
static void ok(int c, const char *w){ printf("  [%s] %s\n", c?"PASS":"FAIL", w); if(!c) fails++; }

/* Small deterministic scalars; no RNG dependency, so failures reproduce. */
static void scalar(limb_t *k, int *bits, uint64_t v)
{ memset(k, 0, sizeof(limb_t)*FP_LIMBS); k[0] = (limb_t)v; *bits = 64; }

int main(void)
{
    printf("standalone test [%s], no legacy library\n", ELIPS_CURVE_NAME);

    ep_t P; ep2_t Q;
    ep_generator(&P);
    ep2_generator(&Q);

    ok(ep_on_curve(&P),      "G1 generator is on the curve");
    ok(ep2_on_curve(&Q),     "G2 generator is on the twist");
    ok(ep_in_subgroup(&P),   "[r]G1 == O");
    ok(ep2_in_subgroup(&Q),  "[r]G2 == O");

    fp12_t z, one;
    fp12_set_one(one);
    ok(elips_pairing(z, &P, &Q) == 1, "pairing accepts the generators");
    ok(!fp12_eq(z, one),              "e(G1,G2) != 1");

    {   /* e(P,Q)^r == 1 */
        fp12_t chk;
        fp12_exp(chk, z, ELIPS_ORDER, ELIPS_ORDER_BITS);
        ok(fp12_eq(chk, one), "e(G1,G2)^r == 1");
    }

    {   /* bilinearity in both arguments at once */
        limb_t ka[FP_LIMBS], kb[FP_LIMBS]; int ba, bb;
        scalar(ka, &ba, 0x9E3779B97F4A7C15ULL);
        scalar(kb, &bb, 0xC2B2AE3D27D4EB4FULL);
        ep_t aP; ep2_t bQ;
        ep_mul(&aP, &P, ka, ba);
        ep2_mul(&bQ, &Q, kb, bb);
        ok(ep_in_subgroup(&aP) && ep2_in_subgroup(&bQ),
           "scalar multiples stay in the subgroup");

        fp12_t lhs, rhs;
        ok(elips_pairing(lhs, &aP, &bQ) == 1, "pairing accepts the multiples");
        /* e(P,Q)^(a*b): exponentiate twice rather than multiplying a*b, which
         * would need reduction mod r that this test has no need to implement */
        fp12_exp(rhs, z, ka, ba);
        fp12_exp(rhs, rhs, kb, bb);
        ok(fp12_eq(lhs, rhs), "e([a]G1,[b]G2) == e(G1,G2)^(ab)");
    }

    {   /* the identity must be rejected, not silently paired */
        ep_t inf; ep_set_infinity(&inf);
        fp12_t r;
        ok(elips_pairing(r, &inf, &Q) == 0, "pairing rejects the identity in G1");
        ok(fp12_eq(r, one), "rejected pairing leaves the result at one");
    }

    printf("%s: %s (%d failure%s)\n", ELIPS_CURVE_NAME,
           fails ? "FAILED" : "OK", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}
