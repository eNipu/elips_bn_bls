/*
 * Uniform field elements and scalars from the system CSPRNG.
 *
 * The legacy layer seeded every gmp_randstate_t from time(NULL), which gives an
 * adversary a search space of a few million seeds and therefore no secrecy at
 * all for anything derived from it. Plan section 6, Phase 5, item 4.
 *
 * elips_random_bytes, the source itself, lives in elips/sysrand.h so that the
 * legacy runtime-curve layer can use it too.
 *
 * Every routine returns 0 on success and -1 on failure, and a failure leaves
 * the output zeroed rather than partly filled. Callers must check: a silently
 * non-random key is worse than an error.
 */
#ifndef ELIPS_RANDOM_H
#define ELIPS_RANDOM_H

#include <stddef.h>
#include "elips/fp.h"
#include "elips/sysrand.h"

/* A uniform field element, returned in Montgomery form.
 *
 * Sampled by reducing FP_BYTES + 16 random bytes modulo p, so the deviation
 * from uniform is below 2^-128 -- the standard "extra 128 bits then reduce"
 * construction. Rejection sampling would be exactly uniform but its running
 * time depends on the sample, and this is called with secret output. */
int fp_rand(fp_t r);

/* A uniform scalar in [0, r), written as ELIPS_ORDER_LIMBS little-endian limbs.
 * Same construction and the same bias bound. */
#define ELIPS_ORDER_LIMBS ((ELIPS_ORDER_BITS + 63) / 64)
int elips_random_scalar(limb_t *k);

/* Is k, read as ELIPS_ORDER_LIMBS limbs, strictly less than the group order?
 * Constant time. The trust boundary wants this before a scalar is used: a
 * scalar at or above r is not a distinct group element, and accepting one
 * silently maps two encodings to the same point. */
int elips_scalar_is_reduced(const limb_t *k);

#endif /* ELIPS_RANDOM_H */
