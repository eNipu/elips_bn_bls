/* Does the dedicated doubling agree with the complete one, everywhere?
 *
 * ep2_dbl uses the dedicated formulas for a = 0, which are 4M + 5S against the
 * complete formulas' 7M + 2S. The argument that they are exception-free on the
 * curve is on the function; this is the check, because that routine is what
 * validates attacker-chosen points and an argument is not a check.
 *
 * The reference below is Renes-Costello-Batina Algorithm 9 written out here,
 * so the two implementations share no code. Comparison is by ep2_eq and not by
 * memcmp: the formulas produce different projective representatives of the
 * same point, which is exactly what should be allowed.
 *
 * Inputs covered: random points on the twist, NOT restricted to G2, each in a
 * random projective representative; the point at infinity; and r aliasing p.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elips/pairing.h"
#include "elips/ec.h"

static int fails, checks;
static void ok(int c, const char *what)
{ checks++; if (!c) { fails++; printf("  [FAIL] %s\n", what); } }

/* RCB Algorithm 9, independent of the library's copy. */
static void rcb_dbl(ep2_t *r, const ep2_t *p)
{
    fp2_t t0, t1, t2, x3, y3, z3, b3, b;
    ep2_curve_b(b);
    fp2_add(b3, b, b); fp2_add(b3, b3, b);

    fp2_mul(t0, p->y, p->y);
    fp2_add(z3, t0, t0); fp2_add(z3, z3, z3); fp2_add(z3, z3, z3);
    fp2_mul(t1, p->y, p->z);
    fp2_mul(t2, p->z, p->z);
    fp2_mul(t2, b3, t2);
    fp2_mul(x3, t2, z3);
    fp2_add(y3, t0, t2);
    fp2_mul(z3, t1, z3);
    fp2_add(t1, t2, t2); fp2_add(t2, t1, t2);
    fp2_sub(t0, t0, t2);
    fp2_mul(y3, t0, y3);
    fp2_add(y3, x3, y3);
    fp2_mul(t1, p->x, p->y);
    fp2_mul(x3, t0, t1);
    fp2_add(x3, x3, x3);
    fp2_copy(r->x, x3); fp2_copy(r->y, y3); fp2_copy(r->z, z3);
}

static unsigned long long st = 0x853c49e6748fea9bull;
static unsigned long long xs(void)
{ st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }

static void rand_fp(fp_t r)
{
    limb_t k[FP_LIMBS];
    for (int i = 0; i < FP_LIMBS; i++) k[i] = xs();
    k[FP_LIMBS - 1] >>= 8;              /* comfortably below p */
    fp_from_limbs(r, k);
}
static void rand_fp2(fp2_t r) { rand_fp(r[0]); rand_fp(r[1]); }

/* A random point on the twist. Almost never in G2: the cofactor is enormous,
 * which is the point -- the doubling has to be right off the subgroup too. */
static int rand_point(ep2_t *P)
{
    fp2_t x, y, t, b;
    rand_fp2(x);
    fp2_sqr(t, x); fp2_mul(t, t, x);
    ep2_curve_b(b); fp2_add(t, t, b);
    if (!fp2_sqrt(y, t)) return 0;
    ep2_from_affine(P, x, y);
    return 1;
}

/* Same point, different representative: (X:Y:Z) -> (lX:lY:lZ). */
static void rescale(ep2_t *P)
{
    fp2_t l;
    do { rand_fp2(l); } while (fp2_is_zero(l));
    fp2_mul(P->x, P->x, l);
    fp2_mul(P->y, P->y, l);
    fp2_mul(P->z, P->z, l);
}

int main(int argc, char **argv)
{
    long want = (argc > 1) ? atol(argv[1]) : 200000;
    long tried = 0, used = 0, offsub = 0;

    printf("dedicated ep2_dbl against RCB [%s]\n", ELIPS_CURVE_NAME);

    while (used < want && tried < want * 8) {
        ep2_t P, a, b;
        tried++;
        if (!rand_point(&P)) continue;      /* x^3 + b' was not a square */
        used++;
        if (!ep2_in_subgroup(&P)) offsub++;
        rescale(&P);

        ep2_dbl(&a, &P);
        rcb_dbl(&b, &P);
        if (!ep2_eq(&a, &b)) {
            if (fails < 3) printf("  [FAIL] dedicated != RCB at %ld\n", used);
            fails++; checks++;
        } else checks++;
        if (!ep2_on_curve(&a)) { printf("  [FAIL] output off the curve\n"); fails++; }
        checks++;

        /* aliasing: r and p the same object */
        ep2_t c; ep2_copy(&c, &P);
        ep2_dbl(&c, &c);
        if (!ep2_eq(&c, &b)) { printf("  [FAIL] aliased dbl differs\n"); fails++; }
        checks++;
    }

    /* the point at infinity, which is the one exceptional case that can reach
     * the ladder: acc goes to O whenever the order of Q divides the partial
     * scalar, and every later doubling sees it */
    {
        ep2_t O, d;
        ep2_set_infinity(&O);
        ep2_dbl(&d, &O);
        ok(ep2_is_infinity(&d), "dbl(O) is O");
        ok(!fp2_is_zero(d.y), "dbl(O) is (0:non-zero:0), not (0:0:0)");
        rcb_dbl(&d, &O);
        ok(ep2_is_infinity(&d), "RCB agrees that dbl(O) is O");
    }

    /* and O reached the way the ladder reaches it, from a real point */
    {
        ep2_t P, acc;
        ep2_generator(&P);
        ep2_mul(&acc, &P, ELIPS_ORDER, ELIPS_ORDER_BITS);   /* [r]P = O */
        ok(ep2_is_infinity(&acc), "[r]G is O");
        ep2_t d1, d2;
        ep2_dbl(&d1, &acc);
        rcb_dbl(&d2, &acc);
        ok(ep2_is_infinity(&d1) && ep2_eq(&d1, &d2), "dbl of a computed O agrees");
    }

    printf("  %ld points (%ld outside G2), %d checks, %d failed\n",
           used, offsub, checks, fails);
    if (used == 0) { printf("  no points generated\n"); return 2; }
    if (offsub == 0) { printf("  every point was in G2, test is vacuous\n"); return 2; }
    return fails ? 1 : 0;
}
