/*
 * Extension field tower. See include/elips/fpx.h.
 *
 * Karatsuba throughout, matching the formulas the Phase 0 reference uses, so
 * the committed known-answer vectors validate this code directly.
 *
 * ponytail: no lazy reduction yet. Proper lazy reduction needs double-width
 * accumulators threaded through every routine, and the measured projection said
 * plain Karatsuba over the new fp_t already gets close to the gate. Add it only
 * if the measurement at the end of the phase falls short -- the ceiling is one
 * Montgomery reduction per fp_t operation instead of one per fp2/fp6 operation.
 */
#include "elips/fpx.h"

/* ============================== fp2 ==================================== */

void fp2_set_zero(fp2_t r) { fp_set_zero(r[0]); fp_set_zero(r[1]); }
void fp2_set_one(fp2_t r)  { fp_set_one(r[0]);  fp_set_zero(r[1]); }
void fp2_copy(fp2_t r, const fp2_t a) { fp_copy(r[0], a[0]); fp_copy(r[1], a[1]); }
void fp2_add(fp2_t r, const fp2_t a, const fp2_t b)
{ fp_add(r[0], a[0], b[0]); fp_add(r[1], a[1], b[1]); }
void fp2_sub(fp2_t r, const fp2_t a, const fp2_t b)
{ fp_sub(r[0], a[0], b[0]); fp_sub(r[1], a[1], b[1]); }
void fp2_neg(fp2_t r, const fp2_t a) { fp_neg(r[0], a[0]); fp_neg(r[1], a[1]); }
void fp2_conj(fp2_t r, const fp2_t a) { fp_copy(r[0], a[0]); fp_neg(r[1], a[1]); }
void fp2_mul_fp(fp2_t r, const fp2_t a, const fp_t b)
{ fp_mul(r[0], a[0], b); fp_mul(r[1], a[1], b); }
void fp2_cselect(fp2_t r, const fp2_t a, const fp2_t b, limb_t mask)
{ fp_cselect(r[0], a[0], b[0], mask); fp_cselect(r[1], a[1], b[1], mask); }
int fp2_is_zero(const fp2_t a) { return fp_is_zero(a[0]) & fp_is_zero(a[1]); }
int fp2_eq(const fp2_t a, const fp2_t b) { return fp_eq(a[0], b[0]) & fp_eq(a[1], b[1]); }

void fp2_mul(fp2_t r, const fp2_t a, const fp2_t b)
{
    /* (a0 + a1 u)(b0 + b1 u) = (a0b0 - a1b1) + ((a0+a1)(b0+b1) - a0b0 - a1b1) u,
     * using u^2 = -1. Three multiplies rather than four. */
    fp_t t0, t1, s, t;
    fp_mul(t0, a[0], b[0]);
    fp_mul(t1, a[1], b[1]);
    fp_add(s, a[0], a[1]);
    fp_add(t, b[0], b[1]);
    fp_mul(s, s, t);
    fp_sub(s, s, t0);
    fp_sub(s, s, t1);
    fp_sub(r[0], t0, t1);
    fp_copy(r[1], s);
}

void fp2_sqr(fp2_t r, const fp2_t a)
{
    /* (a0 + a1 u)^2 = (a0+a1)(a0-a1) + 2 a0 a1 u */
    fp_t s, d, m;
    fp_add(s, a[0], a[1]);
    fp_sub(d, a[0], a[1]);
    fp_mul(m, a[0], a[1]);
    fp_mul(r[0], s, d);
    fp_add(r[1], m, m);
}

void fp2_mul_xi(fp2_t r, const fp2_t a)
{
    /* multiply by xi = 1 + u: (a0 + a1 u)(1 + u) = (a0 - a1) + (a0 + a1) u */
    fp_t s, d;
    fp_sub(d, a[0], a[1]);
    fp_add(s, a[0], a[1]);
    fp_copy(r[0], d);
    fp_copy(r[1], s);
}

void fp2_inv(fp2_t r, const fp2_t a)
{
    /* conjugate over the norm; the norm a0^2 + a1^2 lives in Fp */
    fp_t n, t;
    fp_sqr(n, a[0]);
    fp_sqr(t, a[1]);
    fp_add(n, n, t);
    fp_inv(n, n);
    fp_mul(r[0], a[0], n);
    fp_mul(t, a[1], n);
    fp_neg(r[1], t);
}

/* ============================== fp6 ==================================== */

void fp6_set_zero(fp6_t r) { fp2_set_zero(r[0]); fp2_set_zero(r[1]); fp2_set_zero(r[2]); }
void fp6_set_one(fp6_t r)  { fp2_set_one(r[0]);  fp2_set_zero(r[1]); fp2_set_zero(r[2]); }
void fp6_copy(fp6_t r, const fp6_t a)
{ fp2_copy(r[0], a[0]); fp2_copy(r[1], a[1]); fp2_copy(r[2], a[2]); }
void fp6_add(fp6_t r, const fp6_t a, const fp6_t b)
{ fp2_add(r[0],a[0],b[0]); fp2_add(r[1],a[1],b[1]); fp2_add(r[2],a[2],b[2]); }
void fp6_sub(fp6_t r, const fp6_t a, const fp6_t b)
{ fp2_sub(r[0],a[0],b[0]); fp2_sub(r[1],a[1],b[1]); fp2_sub(r[2],a[2],b[2]); }
void fp6_neg(fp6_t r, const fp6_t a)
{ fp2_neg(r[0],a[0]); fp2_neg(r[1],a[1]); fp2_neg(r[2],a[2]); }
int fp6_is_zero(const fp6_t a)
{ return fp2_is_zero(a[0]) & fp2_is_zero(a[1]) & fp2_is_zero(a[2]); }
int fp6_eq(const fp6_t a, const fp6_t b)
{ return fp2_eq(a[0],b[0]) & fp2_eq(a[1],b[1]) & fp2_eq(a[2],b[2]); }

void fp6_mul_v(fp6_t r, const fp6_t a)
{
    /* v*(c0 + c1 v + c2 v^2) = c2*xi + c0 v + c1 v^2, since v^3 = xi */
    fp2_t t;
    fp2_mul_xi(t, a[2]);
    fp2_t c0, c1;
    fp2_copy(c0, a[0]);
    fp2_copy(c1, a[1]);
    fp2_copy(r[0], t);
    fp2_copy(r[1], c0);
    fp2_copy(r[2], c1);
}

void fp6_mul(fp6_t r, const fp6_t a, const fp6_t b)
{
    fp2_t t0, t1, t2, s, t, e0, e1, e2;
    fp2_mul(t0, a[0], b[0]);
    fp2_mul(t1, a[1], b[1]);
    fp2_mul(t2, a[2], b[2]);

    fp2_add(s, a[1], a[2]); fp2_add(t, b[1], b[2]); fp2_mul(s, s, t);
    fp2_sub(s, s, t1); fp2_sub(s, s, t2); fp2_mul_xi(s, s); fp2_add(e0, s, t0);

    fp2_add(s, a[0], a[1]); fp2_add(t, b[0], b[1]); fp2_mul(s, s, t);
    fp2_sub(s, s, t0); fp2_sub(s, s, t1); fp2_mul_xi(t, t2); fp2_add(e1, s, t);

    fp2_add(s, a[0], a[2]); fp2_add(t, b[0], b[2]); fp2_mul(s, s, t);
    fp2_sub(s, s, t0); fp2_sub(s, s, t2); fp2_add(e2, s, t1);

    fp2_copy(r[0], e0); fp2_copy(r[1], e1); fp2_copy(r[2], e2);
}

void fp6_sqr(fp6_t r, const fp6_t a) { fp6_mul(r, a, a); }

void fp6_inv(fp6_t r, const fp6_t a)
{
    /* A = c0^2 - xi c1 c2, B = xi c2^2 - c0 c1, C = c1^2 - c0 c2,
     * F = xi (c2 B + c1 C) + c0 A, inverse = (A,B,C) / F */
    fp2_t A, B, C, F, t;
    fp2_sqr(A, a[0]); fp2_mul(t, a[1], a[2]); fp2_mul_xi(t, t); fp2_sub(A, A, t);
    fp2_sqr(B, a[2]); fp2_mul_xi(B, B); fp2_mul(t, a[0], a[1]); fp2_sub(B, B, t);
    fp2_sqr(C, a[1]); fp2_mul(t, a[0], a[2]); fp2_sub(C, C, t);

    fp2_mul(F, a[2], B);
    fp2_mul(t, a[1], C);
    fp2_add(F, F, t);
    fp2_mul_xi(F, F);
    fp2_mul(t, a[0], A);
    fp2_add(F, F, t);
    fp2_inv(F, F);

    fp2_mul(r[0], A, F);
    fp2_mul(r[1], B, F);
    fp2_mul(r[2], C, F);
}

/* ============================== fp12 =================================== */

void fp12_set_zero(fp12_t r) { fp6_set_zero(r[0]); fp6_set_zero(r[1]); }
void fp12_set_one(fp12_t r)  { fp6_set_one(r[0]);  fp6_set_zero(r[1]); }
void fp12_copy(fp12_t r, const fp12_t a) { fp6_copy(r[0], a[0]); fp6_copy(r[1], a[1]); }
void fp12_add(fp12_t r, const fp12_t a, const fp12_t b)
{ fp6_add(r[0],a[0],b[0]); fp6_add(r[1],a[1],b[1]); }
void fp12_sub(fp12_t r, const fp12_t a, const fp12_t b)
{ fp6_sub(r[0],a[0],b[0]); fp6_sub(r[1],a[1],b[1]); }
void fp12_neg(fp12_t r, const fp12_t a) { fp6_neg(r[0],a[0]); fp6_neg(r[1],a[1]); }
void fp12_conj(fp12_t r, const fp12_t a) { fp6_copy(r[0],a[0]); fp6_neg(r[1],a[1]); }
int fp12_is_zero(const fp12_t a) { return fp6_is_zero(a[0]) & fp6_is_zero(a[1]); }
int fp12_eq(const fp12_t a, const fp12_t b)
{ return fp6_eq(a[0],b[0]) & fp6_eq(a[1],b[1]); }

void fp12_mul(fp12_t r, const fp12_t a, const fp12_t b)
{
    /* (d0 + d1 w)(e0 + e1 w) = (d0e0 + v d1e1)
     *                        + ((d0+d1)(e0+e1) - d0e0 - d1e1) w,  w^2 = v */
    fp6_t t0, t1, s, t;
    fp6_mul(t0, a[0], b[0]);
    fp6_mul(t1, a[1], b[1]);
    fp6_add(s, a[0], a[1]);
    fp6_add(t, b[0], b[1]);
    fp6_mul(s, s, t);
    fp6_sub(s, s, t0);
    fp6_sub(s, s, t1);
    fp6_mul_v(t, t1);
    fp6_add(r[0], t0, t);
    fp6_copy(r[1], s);
}

void fp12_sqr(fp12_t r, const fp12_t a) { fp12_mul(r, a, a); }

void fp12_inv(fp12_t r, const fp12_t a)
{
    /* (d0 + d1 w)^-1 = (d0 - d1 w) / (d0^2 - v d1^2) */
    fp6_t f, t;
    fp6_sqr(f, a[0]);
    fp6_sqr(t, a[1]);
    fp6_mul_v(t, t);
    fp6_sub(f, f, t);
    fp6_inv(f, f);
    fp6_mul(r[0], a[0], f);
    fp6_mul(t, a[1], f);
    fp6_neg(r[1], t);
}
