/*
 * Prime field arithmetic in Montgomery form over fixed-width limb arrays.
 *
 * This is the Phase 3 replacement for the mpz_t-based Fp layer. The two
 * differences that matter:
 *
 *   1. An element is a fixed array of limbs living wherever the caller put it,
 *      usually the stack. The old Fp held one mpz_t, so an Fp12 temporary meant
 *      twelve separate heap allocations and the Miller loop performed roughly
 *      four thousand malloc/free pairs per pairing.
 *
 *   2. Values are kept in Montgomery form, so reduction after a multiply is a
 *      shift-and-multiply-add rather than mpz_mod's division.
 *
 * Every routine here is constant time with respect to its operands: no branch
 * and no memory index depends on a value. Conditional subtraction is done by
 * masked select. fp_inv exponentiates by the public constant p-2, so its
 * control flow depends only on the modulus, which is public.
 *
 * The curve is selected at compile time (plan section 1.5) because FP_LIMBS
 * must be a constant for the loops to unroll and for the Phase 6 assembly to
 * be written against a known width.
 */
#ifndef ELIPS_FP_H
#define ELIPS_FP_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "elips/fp_params.h"

typedef limb_t fp_t[FP_LIMBS];

/* Number of bytes in the canonical little-endian encoding of a field element. */
#define FP_BYTES (FP_LIMBS * 8)

/* --- conversion ------------------------------------------------------------
 * "plain" means an ordinary residue in [0,p); "mont" means a*R mod p.
 * Callers work in Montgomery form throughout and convert only at the edges. */
void fp_from_limbs(fp_t r, const limb_t *plain);   /* plain -> Montgomery */
void fp_to_limbs(limb_t *plain, const fp_t a);     /* Montgomery -> plain */

void fp_set_zero(fp_t r);
void fp_set_one(fp_t r);                           /* the Montgomery 1, i.e. R mod p */
void fp_copy(fp_t r, const fp_t a);

/* --- arithmetic, all constant time ---------------------------------------- */
void fp_add(fp_t r, const fp_t a, const fp_t b);
void fp_sub(fp_t r, const fp_t a, const fp_t b);
void fp_neg(fp_t r, const fp_t a);
void fp_mul(fp_t r, const fp_t a, const fp_t b);   /* Montgomery product */
void fp_sqr(fp_t r, const fp_t a);
void fp_inv(fp_t r, const fp_t a);         /* constant time; 0 maps to 0 */
/* Variable time. Only for values that are already public -- never a secret. */
/* The previous constant-time inversion, mpn_sec_invert based. Kept as the
 * reference fp_inv is checked against, and as the fallback if the batched
 * divstep version ever has to be backed out. Same contract, ~6x slower. */
void fp_inv_sec(fp_t r, const fp_t a);

void fp_inv_vartime(fp_t r, const fp_t a);

/* --- predicates ------------------------------------------------------------
 * Return 1 or 0 without branching on the operands. */
int fp_is_zero(const fp_t a);
int fp_eq(const fp_t a, const fp_t b);

/* Constant-time select: r = mask ? a : b, where mask must be 0 or all ones. */
void fp_cselect(fp_t r, const fp_t a, const fp_t b, limb_t mask);

/* --- exponentiation and square roots ---------------------------------------
 * The exponent is read as ebits little-endian bits and is assumed PUBLIC: the
 * square-and-multiply schedule depends on it. Every exponent used inside the
 * library is derived from p, which is public. Never pass a secret. */
void fp_exp(fp_t r, const fp_t a, const limb_t *e, int ebits);

/* The three exponents derived from the modulus, each FP_LIMBS limbs wide and
 * FP_BITS bits long: (p+1)/4, (p-1)/2 and (p-3)/4. Any argument may be NULL.
 * Exposed because the Fp2 square root needs the same values. */
void fp_exp_constants(limb_t *sqrt_e, limb_t *half_e, limb_t *quarter_e);

/* Square root, for the p = 3 (mod 4) primes this library supports.
 *
 * Returns 1 and writes a root to r when a is a quadratic residue; returns 0 and
 * zeroes r when it is not. The result is the root with the exponentiation's own
 * sign convention, which is arbitrary -- callers that need a specific one (the
 * point decompressor does) must fix it themselves.
 *
 * Constant time with respect to a: one exponentiation by the public constant
 * (p+1)/4, then a comparison. Which of the two roots comes back does depend on
 * a, but that is inherent to the function's output, not a side channel. */
int fp_sqrt(fp_t r, const fp_t a);

/* Is a strictly greater than (p-1)/2 in the canonical representation? This is
 * the "lexicographically largest" predicate the compressed point encodings use
 * to recover the sign of y from a single bit. Constant time. */
int fp_is_lex_largest(const fp_t a);

/* RFC 9380 4.1 sgn0: the low bit of the canonical representation.
 *
 * A different predicate from fp_is_lex_largest, and the two are easy to
 * confuse. Both answer "which of the two square roots is this", but the
 * encodings use "greater than (p-1)/2" and hash-to-curve uses "odd", and
 * substituting one for the other produces points that are on the curve, in the
 * group, and wrong. Constant time. */
int fp_sgn0(const fp_t a);

#endif /* ELIPS_FP_H */
