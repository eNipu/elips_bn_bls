/*
 * One constant-time primitive that the compiler must not see through.
 *
 * THE PROBLEM. A masked select is the standard way to choose between two
 * values without branching:
 *
 *     mask = 0 - condition;            // condition is 0 or 1
 *     r    = (a & mask) | (b & ~mask);
 *
 * It is branch-free in C, and a compiler is entitled to notice that `mask` is
 * provably 0 or ~0, conclude that the expression is just "a or b", and emit a
 * BRANCH on the condition. Both gcc and clang do this, at different sites and
 * in different shapes, and neither warns. The source still looks constant
 * time; the object file is not.
 *
 * This is not hypothetical here. It has cost this repository twice:
 *
 *   src/arith/wide.c  clang turned the reduction's select into
 *                     `neg` then `jae`, and dudect read 113 -> 284 against
 *                     0.99 for the same source under gcc.
 *
 *   glv_divrem        clang turned the divider's select into
 *                     `setb` then `test $1` then `je`, a branch on a borrow
 *                     derived from the secret scalar. That is issue #30: it
 *                     leaked in ep_mul_glv, ep2_mul_glv and fp12_exp_gt,
 *                     reading 66 where gcc read 2, and was blamed on Apple
 *                     silicon for weeks because macOS is where clang runs by
 *                     default.
 *
 * The second one hid from a search for the first, because the same bug wears
 * different instructions depending on how the condition was computed.
 *
 * THE FIX. Pass the mask through an empty asm with a read-write operand. The
 * compiler must then assume the value could be anything, so the select stays a
 * select. It costs one register move and no memory traffic.
 *
 * USE IT ON EVERY MASK whose value the compiler can prove is one of two,
 * which in practice means every mask built as `0 - condition`. The cost is
 * nil and the failure mode is silent.
 */
#ifndef ELIPS_ARITH_CT_H
#define ELIPS_ARITH_CT_H

#include "elips/fp.h"

static inline limb_t ct_mask(limb_t m)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ ("" : "+r"(m));
    return m;
#else
    /* No GNU asm. A volatile round trip is weaker -- it survives most
     * optimisers but is not guaranteed -- and is better than nothing. A port
     * to such a compiler should check the disassembly rather than trust this. */
    volatile limb_t v = m;
    return v;
#endif
}

#endif /* ELIPS_ARITH_CT_H */
