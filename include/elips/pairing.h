/*
 * Optimal ate pairing.
 *
 * G1 is E(Fp), G2 is represented on the sextic twist E'(Fp2), and the result
 * lives in the order-r subgroup of Fp12.
 */
#ifndef ELIPS_PAIRING_H
#define ELIPS_PAIRING_H

#include <stddef.h>

#include "elips/ec.h"

/* One line function of the Miller loop, with the G1 argument factored out.
 *
 * The line through the running point T evaluated at P has only its w^0, w^3
 * and w^5 coefficients non-zero, and P enters only as
 *
 *     c0 = yP * a      c3 = c      c5 = xP * b
 *
 * so a, b and c depend on the G2 argument alone. That is what makes
 * fixed-argument precomputation possible: everything below is computed once per
 * Q and replayed against as many P as you like. */
typedef struct { fp2_t a, b, c; } ep2_line_t;

/* Every line one Miller loop evaluates, for one fixed Q. ELIPS_MILLER_LINES is
 * derived by the parameter generator from the loop's own digits, not bounded by
 * a guess. */
typedef struct { ep2_line_t l[ELIPS_MILLER_LINES]; } ep2_prec_t;

/* Miller loop only, no final exponentiation. Q must be affine on the twist. */
void pairing_miller(fp12_t f, const fp2_t qx, const fp2_t qy,
                    const fp_t px, const fp_t py);

/* The same loop shared across n pairs: the accumulator is squared once per
 * iteration rather than once per pair. Returns the product of the Miller
 * values, which still needs a final exponentiation. */
void pairing_miller_multi(fp12_t f, const fp2_t *qx, const fp2_t *qy,
                          const fp_t *px, const fp_t *py, size_t n);

/* The Miller loop replayed from a precomputed table. Bit-identical to
 * pairing_miller on the Q the table was built from. */
void pairing_miller_prec(fp12_t f, const ep2_prec_t *pc,
                         const fp_t px, const fp_t py);

/* f^((p^12-1)/r), computed exactly. Slow but definitionally correct; used to
 * validate the fast chain. */
void pairing_final_exp_plain(fp12_t r, const fp12_t f);

/* The fast final exponentiation.
 *
 * IMPORTANT: this yields e^3, not e. That is a property of the standard BLS12
 * addition chain, not a defect. The chain exists because
 *
 *     3*lambda = (x-1)^2 * (x+p) * (x^2+p^2-1) + 3
 *
 * has a short evaluation while lambda itself does not; recovering e would need
 * a cube root in mu_r, which costs a full-width exponentiation. Verified
 * exactly on BLS12-381 and BLS12-461. The legacy library and RELIC both do the
 * same thing. Since gcd(3, r) = 1, e^3 is still a bilinear, non-degenerate
 * pairing, and any protocol that only compares pairings is unaffected --
 * but raw values will not match an implementation that outputs e.
 *
 * See issue #16. Use pairing_final_exp_plain when the exact value is required. */
/* On BLS12 this returns e^3 (see the note above). On BN it returns e exactly:
 * the BN hard part decomposes with no stray factor. */
void pairing_final_exp_fast(fp12_t r, const fp12_t f);

/* f^x with x the curve's mother parameter, for f in the cyclotomic subgroup
 * (where conjugation is inversion). */
void fp12_exp_param(fp12_t r, const fp12_t f);

/* --- public entry points --------------------------------------------------
 * These make the new layer usable on its own. Until they existed the tests had
 * to borrow points from the legacy library, which is the last thing keeping
 * that layer alive. */

/* The group generators, verified at generation time to be on their curves and
 * of order exactly r.
 *
 * On BLS12-381 these are the generators the specification publishes, taken from
 * it rather than searched for: a library that hands out a different generator
 * of the same subgroup does not interoperate even with a byte-exact encoding,
 * because signatures verify against the standard generator or not at all.
 * test/serialize_test.c pins both against the published encodings.
 *
 * BLS12-461 and BN-462 have no specification to follow, so their generators are
 * the first ones the search in tools/reference/gen_params.py finds. */
void ep_generator(ep_t *g);
void ep2_generator(ep2_t *g);

/* Subgroup membership. A point can sit on the curve and still be outside the
 * order-r subgroup, which is what small-subgroup attacks exploit. Nothing in
 * the legacy library ever checked this.
 *
 * These are endomorphism tests, not [r]P: see the derivation in
 * tools/reference/subgroup_ref.py and the note above their definitions. */
int ep_in_subgroup(const ep_t *p);
int ep2_in_subgroup(const ep2_t *q);

/* The full pairing: Miller loop then the fast final exponentiation.
 *
 * Returns 0 and leaves the result at one if either input is the identity or
 * fails its subgroup check, so a caller that ignores the return value gets a
 * useless answer rather than a subtly wrong one.
 *
 * On BLS12 the value is e^3 (see the note on pairing_final_exp_fast); on BN it
 * is e exactly. */
int elips_pairing(fp12_t out, const ep_t *P, const ep2_t *Q);

/* The product of n pairings, prod_j e(P_j, Q_j), with one shared Miller loop
 * and a single final exponentiation.
 *
 * This is what a verification equation actually needs. Checking
 * e(H(m), pk) * e(-sig, G2) == 1 as two separate pairings runs the final
 * exponentiation twice and throws away about half the work; here it runs once
 * however many terms there are, and the per-iteration squaring is shared
 * across the terms too.
 *
 * Same contract as elips_pairing: returns 0 and leaves the result at one if
 * any input is the identity or fails its subgroup check, and every input is
 * validated before any of them is used. On BLS12 the value is
 * (prod_j e(P_j, Q_j))^3, on BN it is the product exactly -- in both cases
 * identical to multiplying together what elips_pairing returns for each pair,
 * because raising to a fixed exponent is a homomorphism.
 *
 * n may be any length; the implementation uses a fixed-size buffer and
 * processes longer inputs in chunks, which never costs more than one
 * final exponentiation. */
int elips_pairing_multi(fp12_t out, const ep_t *P, const ep2_t *Q, size_t n);

/* Fixed-argument precomputation (Costello and Stebila, ePrint 2010/342).
 *
 * In a verification the G2 arguments are fixed -- the generator, a long-lived
 * public key -- so the entire G2 side of the Miller loop can be computed once
 * and replayed against as many G1 points as you like. What is saved is the
 * point arithmetic, which is most of an iteration's work.
 *
 * ep2_precompute checks Q's subgroup membership and returns 0 if it fails,
 * because a table has no point left to check afterwards. The pairing calls
 * then check only P.
 *
 * A precomputed pairing returns exactly what the direct one does, bit for bit,
 * not merely the same value after the final exponentiation. */
int  ep2_precompute(ep2_prec_t *pc, const ep2_t *Q);
int  elips_pairing_prec(fp12_t out, const ep_t *P, const ep2_prec_t *pc);
int  elips_pairing_multi_prec(fp12_t out, const ep_t *P,
                              const ep2_prec_t *pc, size_t n);

#endif /* ELIPS_PAIRING_H */
