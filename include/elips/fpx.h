/*
 * Extension field tower on top of the fixed-width Montgomery Fp layer.
 *
 *   fp2  = Fp[u]  / (u^2 + 1)          a0 + a1*u
 *   fp6  = fp2[v] / (v^3 - (1+u))      c0 + c1*v + c2*v^2
 *   fp12 = fp6[w] / (w^2 - v)          d0 + d1*w
 *
 * Same tower the legacy layer uses, so the Phase 0 vectors apply unchanged.
 * Elements are plain nested arrays: an fp12_t is 12 contiguous fp_t with no
 * indirection and no allocation, where the old Fp12 was twelve separate heap
 * mpz_t.
 *
 * Every routine tolerates aliasing of output with input.
 */
#ifndef ELIPS_FPX_H
#define ELIPS_FPX_H

#include "elips/fp.h"

typedef fp_t  fp2_t[2];
typedef fp2_t fp6_t[3];
typedef fp6_t fp12_t[2];

/* ---- fp2 ---- */
void fp2_set_zero(fp2_t r);
void fp2_set_one(fp2_t r);
void fp2_copy(fp2_t r, const fp2_t a);
void fp2_add(fp2_t r, const fp2_t a, const fp2_t b);
void fp2_sub(fp2_t r, const fp2_t a, const fp2_t b);
void fp2_neg(fp2_t r, const fp2_t a);
void fp2_mul(fp2_t r, const fp2_t a, const fp2_t b);
void fp2_sqr(fp2_t r, const fp2_t a);
void fp2_inv(fp2_t r, const fp2_t a);
void fp2_mul_xi(fp2_t r, const fp2_t a);      /* multiply by 1+u */
void fp2_conj(fp2_t r, const fp2_t a);        /* the p-power Frobenius */
void fp2_mul_fp(fp2_t r, const fp2_t a, const fp_t b);
int  fp2_is_zero(const fp2_t a);
int  fp2_eq(const fp2_t a, const fp2_t b);
void fp2_cselect(fp2_t r, const fp2_t a, const fp2_t b, limb_t mask);

/* ---- fp6 ---- */
void fp6_set_zero(fp6_t r);
void fp6_set_one(fp6_t r);
void fp6_copy(fp6_t r, const fp6_t a);
void fp6_add(fp6_t r, const fp6_t a, const fp6_t b);
void fp6_sub(fp6_t r, const fp6_t a, const fp6_t b);
void fp6_neg(fp6_t r, const fp6_t a);
void fp6_mul(fp6_t r, const fp6_t a, const fp6_t b);
void fp6_sqr(fp6_t r, const fp6_t a);
void fp6_inv(fp6_t r, const fp6_t a);
void fp6_mul_v(fp6_t r, const fp6_t a);       /* multiply by v */
int  fp6_is_zero(const fp6_t a);
int  fp6_eq(const fp6_t a, const fp6_t b);

/* ---- fp12 ---- */
void fp12_set_zero(fp12_t r);
void fp12_set_one(fp12_t r);
void fp12_copy(fp12_t r, const fp12_t a);
void fp12_add(fp12_t r, const fp12_t a, const fp12_t b);
void fp12_sub(fp12_t r, const fp12_t a, const fp12_t b);
void fp12_neg(fp12_t r, const fp12_t a);
void fp12_mul(fp12_t r, const fp12_t a, const fp12_t b);
void fp12_sqr(fp12_t r, const fp12_t a);
void fp12_inv(fp12_t r, const fp12_t a);
void fp12_conj(fp12_t r, const fp12_t a);     /* the p^6 Frobenius */
int  fp12_is_zero(const fp12_t a);
int  fp12_eq(const fp12_t a, const fp12_t b);

#endif /* ELIPS_FPX_H */
