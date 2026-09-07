/*
 * Elliptic curve points in homogeneous projective coordinates, with the
 * complete addition formulas of Renes, Costello and Batina.
 *
 * Pulled forward from Phase 4 into Phase 3 because the two are coupled: the
 * affine group law needs a field inversion per point operation, and constant-
 * time inversion costs 26.8 us against 1.4 us for the variable-time kind. With
 * roughly 79 inversions in a BLS12 Miller loop that is 2.1 ms of inversion
 * against a 2.49 ms loop, so constant-time arithmetic is unaffordable until the
 * inversions are gone. Jacobian coordinates remove all but one of them.
 *
 * A point is (X : Y : Z) with x = X/Z and y = Y/Z. The identity has Z = 0.
 * Both curves have a = 0, so the RCB formulas specialise to their cheapest form.
 *
 * ep_add and ep2_add are correct for every input pair -- equal, opposite,
 * identity -- with no branch and no fixup. There is deliberately no faster
 * incomplete variant: the Jacobian version had one, guarded by a precondition
 * the caller had to honour, and a routine that is wrong for inputs a caller can
 * plausibly supply is a defect waiting for its first careless call site.
 */
#ifndef ELIPS_EC_H
#define ELIPS_EC_H

#include "elips/fpx.h"

typedef struct { fp_t  x, y, z; } ep_t;    /* E  over Fp  */
typedef struct { fp2_t x, y, z; } ep2_t;   /* E' over Fp2, the sextic twist */

/* --- E(Fp) --- */
void ep_set_infinity(ep_t *r);
int  ep_is_infinity(const ep_t *p);
void ep_copy(ep_t *r, const ep_t *p);
void ep_neg(ep_t *r, const ep_t *p);
void ep_dbl(ep_t *r, const ep_t *p);
void ep_add(ep_t *r, const ep_t *p, const ep_t *q);
void ep_from_affine(ep_t *r, const fp_t x, const fp_t y);
int  ep_to_affine(fp_t x, fp_t y, const ep_t *p);   /* 0 if p is infinity */
void ep_mul(ep_t *r, const ep_t *p, const limb_t *k, int kbits);
int  ep_on_curve(const ep_t *p);
/* Projective equality: no inversion, and correct for infinity. */
int  ep_eq(const ep_t *a, const ep_t *b);
#ifdef ELIPS_FAMILY_BLS12
/* The GLV endomorphism (x, y) -> (beta*x, y). On G1 it acts as [-x^2];
 * on the rest of E(Fp) it does not, which is what makes it a subgroup test. */
void ep_phi(ep_t *r, const ep_t *p);
#endif

/* --- E'(Fp2) --- */
void ep2_set_infinity(ep2_t *r);
int  ep2_is_infinity(const ep2_t *p);
void ep2_copy(ep2_t *r, const ep2_t *p);
void ep2_neg(ep2_t *r, const ep2_t *p);
void ep2_dbl(ep2_t *r, const ep2_t *p);
void ep2_add(ep2_t *r, const ep2_t *p, const ep2_t *q);
void ep2_from_affine(ep2_t *r, const fp2_t x, const fp2_t y);
int  ep2_to_affine(fp2_t x, fp2_t y, const ep2_t *p);
void ep2_mul(ep2_t *r, const ep2_t *p, const limb_t *k, int kbits);
/* The skew Frobenius on the twist. On G2 it acts as multiplication by a fixed
 * eigenvalue, which is what makes GLV possible. */
void ep2_psi(ep2_t *r, const ep2_t *p);

/* GLV scalar multiplication on G2. Four-dimensional via psi on BLS12,
 * two-dimensional on BN.
 *
 *
 * On BLS12, psi acts on G2 as multiplication by the mother parameter x, which
 * is only 64 to 77 bits against a 255 to 308 bit group order. Writing the
 * scalar in base |x| therefore gives four short digits and cuts the ladder to a
 * quarter of its length. On BN, psi acts as 6x^2 -- 231 bits against a 462-bit
 * order, exactly sqrt(r) -- so the split is base 6x^2 and two-dimensional,
 * halving the ladder rather than quartering it.
 *
 * Constant time throughout, decomposition included: the division is restoring,
 * one bit at a time, with the conditional subtraction done by mask. Safe for
 * secret scalars.
 *
 * PRECONDITION: q must lie in G2, the order-r subgroup. psi acts as
 * multiplication by x only on that eigenspace; on an arbitrary point of the
 * twist it does not, and this routine then returns a wrong answer rather than
 * failing. Use ep2_in_subgroup if the caller cannot guarantee it, or ep2_mul,
 * which is general-purpose and has no such requirement. This is why ep2_mul
 * does not simply dispatch here. */
void ep2_mul_glv(ep2_t *r, const ep2_t *q, const limb_t *k, int kbits);

#ifdef ELIPS_FAMILY_BLS12
/* Two-dimensional GLV on G1, using phi. Same precondition: p must be in G1.
 *
 * BLS12 only. phi acts as [-x^2] there and |x^2| is sqrt(r), so the scalar
 * splits in base x^2 with no lattice reduction. BN's lambda is 348 bits
 * against sqrt(r) = 231, so no such split exists and ep_mul stays the routine
 * to use there. */
void ep_mul_glv(ep_t *r, const ep_t *p, const limb_t *k, int kbits);
#endif
int  ep2_on_curve(const ep2_t *p);
int  ep2_eq(const ep2_t *a, const ep2_t *b);

/* Curve constant b in Montgomery form, and the twist constant b*xi. */
void ep_curve_b(fp_t b);
void ep2_curve_b(fp2_t b);

#define ELIPS_HAVE_EC 1
#endif /* ELIPS_EC_H */
