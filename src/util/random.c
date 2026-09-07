/*
 * Uniform field elements and scalars. See include/elips/random.h.
 *
 * The bytes come from src/util/sysrand.c; what is here is the part that turns
 * them into elements of the right set without introducing bias or a timing
 * dependence on the sample.
 */
#include "elips/random.h"

#include <string.h>
#include <gmp.h>

/* Reduce a wide random value modulo a public modulus, in constant time.
 *
 * mpn_sec_div_r is GMP's side-channel-silent remainder: its running time
 * depends on the operand sizes, which are compile-time constants here, and not
 * on the values. mpn_tdiv_qr would be faster and would leak.
 *
 * Returns 0 on success. The scratch bound is checked rather than assumed: an
 * itch larger than the buffer is a wrong answer waiting to happen, and this
 * codebase has already paid for one fixed-size buffer that was not checked
 * (defect M7).
 */
static int reduce_wide(limb_t *wide, int nn, const limb_t *mod, int dn)
{
    mp_limb_t scratch[512];
    mp_size_t itch = mpn_sec_div_r_itch(nn, dn);
    if (itch > (mp_size_t)(sizeof scratch / sizeof scratch[0])) return -1;
    mpn_sec_div_r(wide, nn, mod, dn, scratch);
    return 0;
}

/* The extra 128 bits that make the modular reduction statistically uniform. */
#define RAND_EXTRA_LIMBS 2

int fp_rand(fp_t r)
{
    limb_t wide[FP_LIMBS + RAND_EXTRA_LIMBS];

    fp_set_zero(r);
    if (elips_random_bytes(wide, sizeof wide) != 0) return -1;
    if (reduce_wide(wide, FP_LIMBS + RAND_EXTRA_LIMBS, FP_MODULUS, FP_LIMBS) != 0)
        return -1;

    fp_from_limbs(r, wide);          /* low FP_LIMBS limbs hold the remainder */
    return 0;
}

int elips_random_scalar(limb_t *k)
{
    limb_t wide[ELIPS_ORDER_LIMBS + RAND_EXTRA_LIMBS];

    memset(k, 0, ELIPS_ORDER_LIMBS * sizeof(limb_t));
    if (elips_random_bytes(wide, sizeof wide) != 0) return -1;
    if (reduce_wide(wide, ELIPS_ORDER_LIMBS + RAND_EXTRA_LIMBS,
                    ELIPS_ORDER, ELIPS_ORDER_LIMBS) != 0)
        return -1;

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
