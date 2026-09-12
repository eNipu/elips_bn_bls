/*
 * Uniform field elements and scalars. See include/elips/random.h.
 *
 * The bytes come from src/util/sysrand.c; what is here is the part that turns
 * them into elements of the right set without introducing bias or a timing
 * dependence on the sample.
 */
#include "elips/random.h"

#include <string.h>

#include "arith/wide.h"

/* The reduction is elips_mod_wide: constant time in the values, with the
 * running time set by the operand sizes, which are compile-time constants at
 * both call sites below. It replaced GMP's mpn_sec_div_r, which had the same
 * contract; test/wide_test.c checks the replacement against GMP, which is
 * still linked into the tests for exactly this reason.
 */

/* The extra 128 bits that make the modular reduction statistically uniform. */
#define RAND_EXTRA_LIMBS 2

int fp_rand(fp_t r)
{
    limb_t wide[FP_LIMBS + RAND_EXTRA_LIMBS];

    fp_set_zero(r);
    if (elips_random_bytes(wide, sizeof wide) != 0) return -1;
    elips_mod_wide(wide, FP_LIMBS + RAND_EXTRA_LIMBS, FP_MODULUS, FP_LIMBS);

    fp_from_limbs(r, wide);          /* low FP_LIMBS limbs hold the remainder */
    return 0;
}

int elips_random_scalar(limb_t *k)
{
    limb_t wide[ELIPS_ORDER_LIMBS + RAND_EXTRA_LIMBS];

    memset(k, 0, ELIPS_ORDER_LIMBS * sizeof(limb_t));
    if (elips_random_bytes(wide, sizeof wide) != 0) return -1;
    elips_mod_wide(wide, ELIPS_ORDER_LIMBS + RAND_EXTRA_LIMBS,
                   ELIPS_ORDER, ELIPS_ORDER_LIMBS);

    memcpy(k, wide, ELIPS_ORDER_LIMBS * sizeof(limb_t));
    return 0;
}

int elips_scalar_is_reduced(const limb_t *k)
{
    /* Borrow out of k - r is 1 exactly when k < r. No branch on k. */
    limb_t borrow = 0;
    for (int i = 0; i < ELIPS_ORDER_LIMBS; i++) {
        limb_t ki = k[i], ri = ELIPS_ORDER[i];
        limb_t d  = ki - ri;
        limb_t b1 = (ki < ri);
        limb_t d2 = d - borrow;
        limb_t b2 = (d < borrow);
        (void)d2;
        borrow = b1 | b2;
    }
    return (int)borrow;
}
