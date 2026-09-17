/*
 * Does the projective group law agree with affine arithmetic, everywhere?
 *
 * ep2 is in Jacobian coordinates (x = X/Z^2, y = Y/Z^3), where the doubling is
 * 2M + 5S and the addition is the unified add-or-double form. Neither is
 * complete the way the RCB formulas in ec_homog.h are, so both carry an
 * exception argument on their definition. This is the check on those
 * arguments, because ep2_dbl and ep2_add are what ep2_in_subgroup runs on
 * points an attacker chose, and an argument is not a check.
 *
 * THE ORACLE IS AFFINE, with a real field inversion and explicit special
 * cases. That is deliberately not a second projective formula: it shares no
 * algebra with the code under test, it is slow and obviously correct, and it
 * cannot be wrong in the same direction.
 *
 * Inputs: random points on the twist NOT restricted to G2, each in a random
 * projective representative, and every degenerate pair -- P+P, P+(-P), P+O,
 * O+P, O+O, dbl(O) -- plus the aliasing cases.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elips/pairing.h"
#include "elips/ec.h"

static int fails, checks;
static void ok(int c, const char *what)
{ checks++; if (!c) { fails++; printf("  [FAIL] %s\n", what); } }

/* ---- the affine oracle. inf is carried as a flag, not as a coordinate ---- */
typedef struct { fp2_t x, y; int inf; } aff_t;

static void aff_dbl(aff_t *r, const aff_t *p)
{
    if (p->inf || fp2_is_zero(p->y)) { r->inf = 1; return; }
    fp2_t l, t, x3;
    fp2_sqr(l, p->x);
    fp2_add(t, l, l); fp2_add(l, t, l);          /* 3x^2 */
    fp2_add(t, p->y, p->y);                      /* 2y   */
    fp2_inv(t, t);
    fp2_mul(l, l, t);                            /* lambda */
    fp2_sqr(x3, l);
    fp2_sub(x3, x3, p->x); fp2_sub(x3, x3, p->x);
    fp2_sub(t, p->x, x3);
    fp2_mul(t, t, l);
    fp2_sub(r->y, t, p->y);
    fp2_copy(r->x, x3);
    r->inf = 0;
}

static void aff_add(aff_t *r, const aff_t *p, const aff_t *q)
{
    if (p->inf) { *r = *q; return; }
    if (q->inf) { *r = *p; return; }
    if (fp2_eq(p->x, q->x)) {
        if (fp2_eq(p->y, q->y)) { aff_dbl(r, p); return; }
        r->inf = 1; return;                      /* q = -p */
    }
    fp2_t l, t, x3;
    fp2_sub(l, q->y, p->y);
    fp2_sub(t, q->x, p->x);
    fp2_inv(t, t);
    fp2_mul(l, l, t);
    fp2_sqr(x3, l);
    fp2_sub(x3, x3, p->x); fp2_sub(x3, x3, q->x);
    fp2_sub(t, p->x, x3);
    fp2_mul(t, t, l);
    fp2_sub(r->y, t, p->y);
    fp2_copy(r->x, x3);
    r->inf = 0;
}

static void to_aff(aff_t *a, const ep2_t *P)
{
    a->inf = !ep2_to_affine(a->x, a->y, P);
}

static int aff_eq(const aff_t *a, const aff_t *b)
{
    if (a->inf || b->inf) return a->inf && b->inf;
    return fp2_eq(a->x, b->x) && fp2_eq(a->y, b->y);
}

/* ---- random on-curve points ---- */
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

/* Almost never in G2: the cofactor is enormous, which is the point -- the
 * group law has to be right off the subgroup too. */
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

/* Same point, a different representative: (X:Y:Z) -> (l^2 X : l^3 Y : l Z). */
static void rescale(ep2_t *P)
{
    fp2_t l, l2, l3;
    do { rand_fp2(l); } while (fp2_is_zero(l));
    fp2_sqr(l2, l); fp2_mul(l3, l2, l);
    fp2_mul(P->x, P->x, l2);
    fp2_mul(P->y, P->y, l3);
    fp2_mul(P->z, P->z, l);
}

static void check_dbl(const ep2_t *P, const char *what)
{
    ep2_t d; aff_t a, b, pa;
    ep2_dbl(&d, P);
    to_aff(&a, &d);
    to_aff(&pa, P);
    aff_dbl(&b, &pa);
    ok(aff_eq(&a, &b), what);
    ok(ep2_on_curve(&d), "dbl output is on the curve");
    /* aliasing */
    ep2_t c; ep2_copy(&c, P); ep2_dbl(&c, &c);
    ok(ep2_eq(&c, &d), "aliased dbl agrees");
}

static void check_add(const ep2_t *P, const ep2_t *Q, const char *what)
{
    ep2_t s; aff_t a, b, pa, qa;
    ep2_add(&s, P, Q);
    to_aff(&a, &s);
    to_aff(&pa, P); to_aff(&qa, Q);
    aff_add(&b, &pa, &qa);
    ok(aff_eq(&a, &b), what);
    ok(ep2_on_curve(&s), "add output is on the curve");
    /* commutative, and both aliasing directions */
    ep2_t t; ep2_add(&t, Q, P);
    ok(ep2_eq(&t, &s), "add is commutative");
    ep2_copy(&t, P); ep2_add(&t, &t, Q);
    ok(ep2_eq(&t, &s), "add aliasing the first argument agrees");
    ep2_copy(&t, Q); ep2_add(&t, P, &t);
    ok(ep2_eq(&t, &s), "add aliasing the second argument agrees");
}

int main(int argc, char **argv)
{
    long want = (argc > 1) ? atol(argv[1]) : 20000;
    long used = 0, offsub = 0;

    printf("projective group law against affine [%s]\n", ELIPS_CURVE_NAME);

    ep2_t O; ep2_set_infinity(&O);

    for (long tried = 0; used < want && tried < want * 8; tried++) {
        ep2_t P, Q, N;
        if (!rand_point(&P)) continue;
        used++;
        if (!ep2_in_subgroup(&P)) offsub++;
        rescale(&P);

        check_dbl(&P, "dbl matches affine");

        if (rand_point(&Q)) {
            rescale(&Q);
            check_add(&P, &Q, "add matches affine");
        }

        /* the degenerate pairs, on every point rather than once */
        check_add(&P, &P, "P + P matches affine (the doubling path)");
        ep2_neg(&N, &P);
        check_add(&P, &N, "P + (-P) is the identity");
        check_add(&P, &O, "P + O is P");
        check_add(&O, &P, "O + P is P");
    }

    check_dbl(&O, "dbl(O) is O");
    check_add(&O, &O, "O + O is O");
    {
        ep2_t d; ep2_dbl(&d, &O);
        ok(ep2_is_infinity(&d), "dbl(O) reports infinity");
        ep2_add(&d, &O, &O);
        ok(ep2_is_infinity(&d), "O + O reports infinity");
    }
    /* O reached the way a ladder reaches it, from a real point */
    {
        ep2_t G, acc, d;
        ep2_generator(&G);
        ep2_mul(&acc, &G, ELIPS_ORDER, ELIPS_ORDER_BITS);
        ok(ep2_is_infinity(&acc), "[r]G is O");
        ep2_dbl(&d, &acc);
        ok(ep2_is_infinity(&d), "dbl of a computed O is O");
        ep2_add(&d, &acc, &G);
        ok(ep2_eq(&d, &G), "a computed O plus G is G");
    }

    printf("  %ld points (%ld outside G2), %d checks, %d failed\n",
           used, offsub, checks, fails);
    if (used == 0) { printf("  no points generated\n"); return 2; }
    if (offsub == 0) { printf("  every point was in G2, test is vacuous\n"); return 2; }
    return fails ? 1 : 0;
}
