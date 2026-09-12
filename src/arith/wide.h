/*
 * Constant-time remainder of a wide value by a public modulus.
 *
 * Internal to the library; not part of the public API. Two callers need it and
 * they need it at different widths and different moduli, which is why it lives
 * here rather than inside either of them:
 *
 *   src/util/random.c        a wide random draw mod p, and mod the group order
 *   src/hash/hash_to_curve.c an L-byte hash mod p, RFC 9380 section 5.2
 *
 * This replaces GMP's mpn_sec_div_r and keeps its contract: the running time
 * depends on the operand SIZES, which are compile-time constants at every call
 * site, and never on the values. That matters most in hash_to_curve, where the
 * value being reduced comes from a message that is usually the thing the
 * caller is trying to hide.
 */
#ifndef ELIPS_ARITH_WIDE_H
#define ELIPS_ARITH_WIDE_H

#include "elips/fp.h"

/* w (nn limbs, little-endian) <- w mod m (dn limbs, little-endian).
 *
 * The remainder is left in the low dn limbs of w and the rest is zeroed.
 * Requires nn >= dn >= 1 and the top limb of m non-zero, which every call site
 * satisfies with compile-time constants. */
void elips_mod_wide(limb_t *w, int nn, const limb_t *m, int dn);

#endif /* ELIPS_ARITH_WIDE_H */
