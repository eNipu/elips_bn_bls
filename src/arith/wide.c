/*
 * Schoolbook binary long division, one bit per step. See src/arith/wide.h.
 *
 * Chosen over Barrett or a Montgomery-based reduction for two reasons, in this
 * order. It is obviously constant time: the loop bound is a compile-time
 * constant, every step does identical work, and the only conditional is a
 * masked select. And it needs nothing precomputed, which matters because the
 * call sites use two different moduli -- the field prime and the group order --
 * and only the prime has Montgomery constants generated for it.
 *
 * WHAT IT COSTS, measured rather than guessed, on BLS12-381 at the
 * hash_to_field width (an 8-limb value mod a 6-limb p), all in one binary:
 *
 *   this routine            2607 ns
 *   GMP mpn_sec_div_r         85 ns      what it replaced, 30x faster
 *
 * Thirty times slower is a real cost and is stated here rather than buried.
 * What it is worth end to end is small, because the callers are not hot loops:
 * hash_to_g1 does two of these inside ~190 us, so about 3%, and hash_to_g2
 * does four inside ~1100 us, so about 1%. Nothing else in the library calls it
 * -- the pairing, the ladders and the final exponentiation are unchanged. That
 * was the trade: three percent of one operation to delete a hard dependency on
 * GMP, which is what stood between this library and a browser build, a wheel
 * and an npm package.
 *
 * Two things were tried and measured before settling here. Specialising the
 * widths as compile-time literals so the compiler could unroll gained 3%: the
 * cost is the algorithm, not the code generation. Seeding the accumulator with
 * the top dn-1 limbs instead of shifting them in one bit at a time gained 2.6x
 * (6673 ns -> 2607 ns) and is kept.
 *
 * If a caller ever makes this hot, the next step is Barrett with a per-modulus
 * precomputed mu, or Horner over limbs using one new constant 2^64 * R mod p,
 * which would land near 250 ns. Neither is worth a generated constant and a
 * second code path today. Measure first.
 */
#include "elips/fp.h"
#include "wide.h"

#include <string.h>

/* Hide from the optimiser that this value is one of exactly two.
 *
 * The masked select at the bottom of the loop is written branch-free, and at
 * -O2 clang compiled it into a branch anyway:
 *
 *     neg  %r11          ; r11 is the 0-or-1 borrow, so CF = (borrow != 0)
 *     jae  ...           ; jump on the mask, into one of two copy loops
 *
 * It could do that because `0 - borrow` is provably 0 or ~0, so the two arms
 * of the select are provably "take one" or "take the other". dudect caught it:
 * the same source read max|t| = 0.99 under gcc and 113 -> 284, growing with
 * the sample, under clang. Growth with n is the signature of a real systematic
 * bias rather than noise.
 *
 * An empty asm with the value as a read-write operand makes the mask opaque:
 * the compiler must assume it could be anything, so the select stays a select.
 * Costs one register move. The fallback path for a compiler without GNU asm is
 * a volatile round trip, which is weaker but still defeats the folding.
 */
static inline limb_t ct_opaque(limb_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ ("" : "+r"(x));
    return x;
#else
    volatile limb_t v = x;
    return v;
#endif
}

#define WIDE_MAX (2 * (int)FP_LIMBS + 2)

void elips_mod_wide(limb_t *w, int nn, const limb_t *m, int dn)
{
    /* acc needs one limb more than the modulus: the shift below can push it to
     * just under 2m before the subtraction pulls it back, and 2m does not
     * always fit in dn limbs. */
    limb_t acc[WIDE_MAX + 1];
    limb_t sub[WIDE_MAX + 1];

    if (nn > WIDE_MAX || dn < 1 || dn > nn) {   /* unreachable at the call sites */
        memset(w, 0, (size_t)nn * sizeof(limb_t));
        return;
    }
    memset(acc, 0, sizeof acc);

    /* Seed with the top dn-1 limbs instead of shifting them in a bit at a
     * time. They cannot exceed the modulus: dn-1 limbs are below 2^(64(dn-1)),
     * and m has a non-zero top limb so m >= 2^(64(dn-1)). That leaves only
     * (nn-dn+1)*64 bits to walk. At the widths used here it is the difference
     * between 512 iterations and 192, which is most of the cost. */
    int seed = dn - 1;
    for (int j = 0; j < seed; j++) acc[j] = w[nn - seed + j];

    for (int bit = (nn - seed) * 64 - 1; bit >= 0; bit--) {
        /* acc = acc*2 + the next bit of w, most significant first. */
        limb_t carry = (w[bit / 64] >> (bit % 64)) & 1;
        for (int j = 0; j <= dn; j++) {
            limb_t t = acc[j];
            acc[j] = (t << 1) | carry;
            carry  = t >> 63;
        }

        /* acc < 2m here, so ONE conditional subtraction restores acc < m. */
        limb_t borrow = 0;
        for (int j = 0; j <= dn; j++) {
            limb_t a  = acc[j];
            limb_t b  = (j < dn) ? m[j] : (limb_t)0;
            limb_t d  = a - b;
            limb_t b1 = (a < b);
            limb_t d2 = d - borrow;
            limb_t b2 = (d < borrow);
            sub[j] = d2;
            borrow = b1 | b2;
        }
        /* borrow == 1 means acc < m: keep acc. Masked, never branched --
         * see ct_opaque above for why the mask has to be laundered. */
        limb_t keep = ct_opaque((limb_t)0 - borrow);
        for (int j = 0; j <= dn; j++)
            acc[j] = (acc[j] & keep) | (sub[j] & ~keep);
    }

    memcpy(w, acc, (size_t)dn * sizeof(limb_t));
    memset(w + dn, 0, (size_t)(nn - dn) * sizeof(limb_t));
}
