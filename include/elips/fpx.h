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
/* r = a * u, i.e. (c0, c1) -> (-c1, c0). */
void fp2_mul_u(fp2_t r, const fp2_t a);
/* Public exponent, same contract as fp_exp. */
void fp2_exp(fp2_t r, const fp2_t a, const limb_t *e, int ebits);
/* Square root over Fp2 for p = 3 (mod 4), by the Adj and Rodriguez-Henriquez
 * method. Returns 1 and a root when a is a square, 0 and zero when it is not.
 * The sign of the returned root is arbitrary; the point decompressor fixes it.
 *
 * Control flow depends only on p, and the one value-dependent choice inside the
 * algorithm is made by masked select rather than a branch. */
int  fp2_sqrt(fp2_t r, const fp2_t a);
/* "Lexicographically largest": c1 > (p-1)/2, or c1 == 0 and c0 > (p-1)/2.
 * The Fp2 analogue of fp_is_lex_largest, and the sign rule the compressed G2
 * encoding uses. Constant time. */
int  fp2_is_lex_largest(const fp2_t a);
/* RFC 9380 4.1 sgn0 for m = 2: sgn0(c0), unless c0 is zero, in which case
 * sgn0(c1). See the note on fp_sgn0. Constant time. */
int  fp2_sgn0(const fp2_t a);
/* Is a a square in Fp2? Equivalent to its norm being a square in Fp, which is
 * one Fp exponentiation rather than an Fp2 one. Constant time. */
int  fp2_is_square(const fp2_t a);

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
/* r = a^(p^k) for k in {1,2,3}. k=6 is fp12_conj. */
void fp12_frobenius(fp12_t r, const fp12_t a, int k);
/* Squaring for elements of the cyclotomic subgroup, where conj(a) = a^-1.
 * Cheaper than fp12_sqr, and WRONG outside that subgroup. Everything the final
 * exponentiation touches after the easy part qualifies. */
void fp12_sqr_cyc(fp12_t r, const fp12_t a);

/* n squarings in a row on a cyclotomic element, via Karabina's compressed
 * form. Cheaper per squaring than fp12_sqr_cyc but costs one Fp2 inversion to
 * come back, so it wins only for a long enough n: ELIPS_KARABINA_MIN_RUN is the
 * measured break-even. Same result as n calls to fp12_sqr_cyc. */
void fp12_sqr_cyc_run(fp12_t r, const fp12_t a, int n);
/* r = a^e, e given as little-endian limbs of ebits bits. Not constant time in
 * the exponent, which is fine: every exponent used here is a public curve
 * parameter. */
void fp12_exp(fp12_t r, const fp12_t a, const limb_t *e, int ebits);

/* r = a^k in G_T, CONSTANT TIME in k, and faster than fp12_exp as well.
 *
 * fp12_exp branches on the exponent, which is fine for the public curve
 * parameters it exists for and unusable for a secret scalar. This is the one
 * to use when k is secret. It also uses cyclotomic squaring and splits k with
 * the Frobenius, which on G_T acts as the same multiplier psi has on G2, so it
 * is roughly 2.6x to 2.8x faster than fp12_exp measured on the same element.
 *
 * PRECONDITION: a must be in G_T, the order-r cyclotomic subgroup. The output
 * of elips_pairing always is. On a general fp12 element this returns a wrong
 * answer rather than failing, exactly as ep_mul_glv does off G1. */
void fp12_exp_gt(fp12_t r, const fp12_t a, const limb_t *k, int kbits);
int  fp12_is_zero(const fp12_t a);
int  fp12_eq(const fp12_t a, const fp12_t b);

#endif /* ELIPS_FPX_H */
