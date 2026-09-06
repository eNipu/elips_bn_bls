/*
 * Optimal ate pairing.
 *
 * G1 is E(Fp), G2 is represented on the sextic twist E'(Fp2), and the result
 * lives in the order-r subgroup of Fp12.
 */
#ifndef ELIPS_PAIRING_H
#define ELIPS_PAIRING_H

#include "elips/ec.h"

/* Miller loop only, no final exponentiation. Q must be affine on the twist. */
void pairing_miller(fp12_t f, const fp2_t qx, const fp2_t qy,
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
void pairing_final_exp_fast(fp12_t r, const fp12_t f);

/* f^x with x the curve's mother parameter, for f in the cyclotomic subgroup
 * (where conjugation is inversion). */
void fp12_exp_param(fp12_t r, const fp12_t f);

#endif /* ELIPS_PAIRING_H */
