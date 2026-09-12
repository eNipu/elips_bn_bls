/*
 * Constant-time scalar decomposition, shared by the GLV routines.
 *
 * Restoring division one bit at a time, with the conditional subtraction done
 * by mask rather than by branch. The loop count comes from a public bit length,
 * never from the value, so the timing reveals nothing about the scalar.
 *
 * This replaces an earlier version that used GMP's mpz_tdiv_qr. That was
 * correct but data dependent, which made the whole GLV path unusable on secret
 * scalars and forced it to be opt-in. It is now the default.
 *
 * Lives in a header rather than in ec.c because the G_T exponentiation in
 * fpx.c needs exactly the same split, and a second copy of a constant-time
 * divider is a copy that can drift.
 *
 * Every helper takes its width n as an argument. n is a compile-time constant
 * at every call site and is never derived from a scalar, so widening changes
 * how much work is done but not whether that amount depends on a secret.
 *
 * All are static inline: a translation unit that uses only some of them does
 * not warn about the rest.
 */
#ifndef ELIPS_GLV_SCALAR_H
#define ELIPS_GLV_SCALAR_H

#include <string.h>
#include "elips/fp.h"
#include "ct.h"

#define GLV_W  ((int)FP_LIMBS)          /* wide enough for r and for k mod r */
#define GLV_WW (2 * (int)FP_LIMBS)      /* wide enough for k times a basis entry */

static inline limb_t glv_sub(limb_t *r, const limb_t *a, const limb_t *b, int n)
{
    limb_t borrow = 0;
    for (int i = 0; i < n; i++) {
        limb_t ai = a[i], bi = b[i];
        limb_t d  = ai - bi;
        limb_t b1 = (ai < bi);
        limb_t d2 = d - borrow;
        limb_t b2 = (d < borrow);
        r[i] = d2;
        borrow = b1 | b2;
    }
    return borrow;
}

static inline void glv_shl1(limb_t *a, limb_t in, int n)
{
    limb_t carry = in;
    for (int i = 0; i < n; i++) {
        limb_t next = a[i] >> 63;
        a[i] = (a[i] << 1) | carry;
        carry = next;
    }
}

/* q = k / d, rem = k mod d. kbits is public. */
static inline void glv_divrem(limb_t *q, limb_t *rem, const limb_t *k, int kbits,
                       const limb_t *d, int n)
{
    limb_t t[GLV_WW];
    memset(q, 0, sizeof(limb_t) * n);
    memset(rem, 0, sizeof(limb_t) * n);

    for (int i = kbits - 1; i >= 0; i--) {
        glv_shl1(rem, (k[i / 64] >> (i % 64)) & 1, n);
        limb_t borrow = glv_sub(t, rem, d, n);
        /* ct_mask, not a bare `0 - (1 - borrow)`. Without it clang proves the
         * mask is one of two values and emits `setb; test $1; je` -- a branch
         * on a borrow derived from the secret scalar. See src/arith/ct.h and
         * issue #30. */
        limb_t mask   = ct_mask((limb_t)0 - (1 - borrow));  /* ones if rem >= d */
        for (int j = 0; j < n; j++)
            rem[j] = (t[j] & mask) | (rem[j] & ~mask);
        q[i / 64] |= (mask & 1) << (i % 64);
    }
}


/* r = a + b mod 2^(64n), returning the carry out. */
static inline limb_t glv_add(limb_t *r, const limb_t *a, const limb_t *b, int n)
{
    limb_t carry = 0;
    for (int i = 0; i < n; i++) {
        limb_t s  = a[i] + b[i];
        limb_t c1 = (s < a[i]);
        limb_t s2 = s + carry;
        limb_t c2 = (s2 < s);
        r[i] = s2;
        carry = c1 | c2;
    }
    return carry;
}

/* r = a * b mod 2^(64n). Schoolbook; the truncation is what makes it modular,
 * and dropping the carry out of the top limb is deliberate. Every operand is
 * read unconditionally, so nothing here branches on a value. */
static inline void glv_mul(limb_t *r, const limb_t *a, const limb_t *b, int n)
{
    limb_t t[GLV_WW];
    memset(t, 0, sizeof(limb_t) * n);
    for (int i = 0; i < n; i++) {
        limb_t carry = 0;
        for (int j = 0; i + j < n; j++) {
            unsigned __int128 s = (unsigned __int128)a[i] * (unsigned __int128)b[j]
                                + (unsigned __int128)t[i + j]
                                + (unsigned __int128)carry;
            t[i + j] = (limb_t)s;
            carry    = (limb_t)(s >> 64);
        }
    }
    memcpy(r, t, sizeof(limb_t) * n);
}

/* r = mask ? -a : a, in two's complement. mask is 0 or all ones. */
static inline void glv_cneg(limb_t *r, const limb_t *a, limb_t mask, int n)
{
    limb_t zero[GLV_WW], t[GLV_WW];
    memset(zero, 0, sizeof(limb_t) * n);
    (void)glv_sub(t, zero, a, n);
    /* The mask arrives as a parameter, but this is a static inline and every
     * caller builds it as `0 - (x >> 63)`, so the optimiser can see straight
     * through to a provable 0 or ~0. Launder it here rather than at each
     * call. */
    mask = ct_mask(mask);
    for (int i = 0; i < n; i++)
        r[i] = (t[i] & mask) | (a[i] & ~mask);
}

#endif /* ELIPS_GLV_SCALAR_H */
