/*
 * Elliptic curve points in Jacobian coordinates.
 *
 * Pulled forward from Phase 4 into Phase 3 because the two are coupled: the
 * affine group law needs a field inversion per point operation, and constant-
 * time inversion costs 26.8 us against 1.4 us for the variable-time kind. With
 * roughly 79 inversions in a BLS12 Miller loop that is 2.1 ms of inversion
 * against a 2.49 ms loop, so constant-time arithmetic is unaffordable until the
 * inversions are gone. Jacobian coordinates remove all but one of them.
 *
 * A point is (X : Y : Z) with x = X/Z^2 and y = Y/Z^3. The point at infinity is
 * Z = 0. Both curves have a = 0 (y^2 = x^3 + b), so the doubling and addition
 * formulas below never reference b.
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
/* Faster, but wrong for p==q, p==-q or infinity. Miller loop only. */
void ep_add_generic(ep_t *r, const ep_t *p, const ep_t *q);
void ep_from_affine(ep_t *r, const fp_t x, const fp_t y);
int  ep_to_affine(fp_t x, fp_t y, const ep_t *p);   /* 0 if p is infinity */
void ep_mul(ep_t *r, const ep_t *p, const limb_t *k, int kbits);
int  ep_on_curve(const ep_t *p);

/* --- E'(Fp2) --- */
void ep2_set_infinity(ep2_t *r);
int  ep2_is_infinity(const ep2_t *p);
void ep2_copy(ep2_t *r, const ep2_t *p);
void ep2_neg(ep2_t *r, const ep2_t *p);
void ep2_dbl(ep2_t *r, const ep2_t *p);
void ep2_add(ep2_t *r, const ep2_t *p, const ep2_t *q);
/* Faster, but wrong for p==q, p==-q or infinity. Miller loop only. */
void ep2_add_generic(ep2_t *r, const ep2_t *p, const ep2_t *q);
void ep2_from_affine(ep2_t *r, const fp2_t x, const fp2_t y);
int  ep2_to_affine(fp2_t x, fp2_t y, const ep2_t *p);
void ep2_mul(ep2_t *r, const ep2_t *p, const limb_t *k, int kbits);
int  ep2_on_curve(const ep2_t *p);

/* Curve constant b in Montgomery form, and the twist constant b*xi. */
void ep_curve_b(fp_t b);
void ep2_curve_b(fp2_t b);

#define ELIPS_HAVE_EC 1
#endif /* ELIPS_EC_H */
