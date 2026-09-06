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

#endif /* ELIPS_PAIRING_H */
