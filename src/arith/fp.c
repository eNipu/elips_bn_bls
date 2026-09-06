/*
 * Montgomery arithmetic over fixed-width limb arrays. See include/elips/fp.h.
 *
 * The multiply uses CIOS (Coarsely Integrated Operand Scanning) from Koc, Acar
 * and Kaliski, "Analyzing and Comparing Montgomery Multiplication Algorithms".
 * CIOS interleaves the product and the reduction, so it needs only s+2 words of
 * scratch and has no data-dependent control flow, which is what makes the
 * constant-time requirement cheap to satisfy here.
 *
 * Deliberately no GMP: the inner loop is plain C over unsigned __int128, which
 * compiles to mul/umulh on AArch64 and mulx on x86-64. That is also exactly the
 * shape Phase 6 replaces with hand-written assembly.
 */
#include "elips/fp.h"
#include <gmp.h>

typedef unsigned __int128 dlimb_t;

#define NLIMB ((int)FP_LIMBS)

/* All-ones when x is non-zero, else zero. Branch-free. */
static inline limb_t mask_nonzero(limb_t x)
{
    return (limb_t)0 - (limb_t)((x | (~x + 1)) >> 63);
}

void fp_cselect(fp_t r, const fp_t a, const fp_t b, limb_t mask)
{
    for (int i = 0; i < NLIMB; i++)
        r[i] = (a[i] & mask) | (b[i] & ~mask);
}

void fp_copy(fp_t r, const fp_t a)    { memcpy(r, a, sizeof(fp_t)); }
void fp_set_zero(fp_t r)              { memset(r, 0, sizeof(fp_t)); }
void fp_set_one(fp_t r)               { memcpy(r, FP_ONE, sizeof(fp_t)); }

/* r = a - b, returning the borrow. Straight-line. */
static inline limb_t sub_borrow(limb_t *r, const limb_t *a, const limb_t *b)
{
    limb_t borrow = 0;
    for (int i = 0; i < NLIMB; i++) {
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

/* r = a + b, returning the carry. Straight-line. */
static inline limb_t add_carry(limb_t *r, const limb_t *a, const limb_t *b)
{
    limb_t carry = 0;
    for (int i = 0; i < NLIMB; i++) {
        limb_t s  = a[i] + b[i];
        limb_t c1 = (s < a[i]);
        limb_t s2 = s + carry;
        limb_t c2 = (s2 < s);
        r[i] = s2;
        carry = c1 | c2;
    }
    return carry;
}

void fp_add(fp_t r, const fp_t a, const fp_t b)
{
    fp_t t;
    limb_t carry  = add_carry(t, a, b);
    fp_t   red;
    limb_t borrow = sub_borrow(red, t, FP_MODULUS);
    /* Take the reduced value when the sum overflowed the limb array, or when
     * subtracting p did not borrow (meaning t >= p). Mask, never branch. */
    limb_t take = (limb_t)0 - (carry | (borrow ^ 1));
    fp_cselect(r, red, t, take);
}

void fp_sub(fp_t r, const fp_t a, const fp_t b)
{
    fp_t t;
    limb_t borrow = sub_borrow(t, a, b);
    fp_t   fixed;
    add_carry(fixed, t, FP_MODULUS);
    fp_cselect(r, fixed, t, (limb_t)0 - borrow);
}

void fp_neg(fp_t r, const fp_t a)
{
    fp_t z;
    fp_set_zero(z);
    /* 0 - 0 must stay 0, and fp_sub already gives that: no borrow, no fixup. */
    fp_sub(r, z, a);
}

/* Montgomery product: r = a*b*R^-1 mod p. */
void fp_mul(fp_t r, const fp_t a, const fp_t b)
{
    limb_t t[FP_LIMBS + 2];
    memset(t, 0, sizeof t);

    for (int i = 0; i < NLIMB; i++) {
        /* t += a * b[i] */
        limb_t c = 0;
        for (int j = 0; j < NLIMB; j++) {
            dlimb_t s = (dlimb_t)a[j] * b[i] + t[j] + c;
            t[j] = (limb_t)s;
            c    = (limb_t)(s >> 64);
        }
        dlimb_t s = (dlimb_t)t[NLIMB] + c;
        t[NLIMB]     = (limb_t)s;
        t[NLIMB + 1] = (limb_t)(s >> 64);

        /* t = (t + m*p) / 2^64, chosen so the low word cancels */
        limb_t m = (limb_t)(t[0] * FP_MONT_N0);
        dlimb_t u = (dlimb_t)m * FP_MODULUS[0] + t[0];
        c = (limb_t)(u >> 64);
        for (int j = 1; j < NLIMB; j++) {
            dlimb_t v = (dlimb_t)m * FP_MODULUS[j] + t[j] + c;
            t[j - 1] = (limb_t)v;
            c        = (limb_t)(v >> 64);
        }
        dlimb_t w = (dlimb_t)t[NLIMB] + c;
        t[NLIMB - 1] = (limb_t)w;
        t[NLIMB]     = t[NLIMB + 1] + (limb_t)(w >> 64);
    }

    /* One conditional subtraction brings the result below p. */
    fp_t red;
    limb_t borrow = sub_borrow(red, t, FP_MODULUS);
    limb_t take   = (limb_t)0 - (t[NLIMB] | (borrow ^ 1));
    fp_cselect(r, red, (const limb_t *)t, take);
}

void fp_sqr(fp_t r, const fp_t a)
{
    /* A dedicated squaring saves roughly a third of the partial products.
     * Deferred: Phase 6 writes it in assembly, and doing it twice is waste.
     * ponytail: reuse fp_mul until a measurement says squaring is the ceiling. */
    fp_mul(r, a, a);
}

void fp_from_limbs(fp_t r, const limb_t *plain)
{
    fp_t t;
    memcpy(t, plain, sizeof(fp_t));
    fp_mul(r, t, FP_R2);            /* a * R^2 * R^-1 = a*R */
}

void fp_to_limbs(limb_t *plain, const fp_t a)
{
    fp_t one;
    fp_set_zero(one);
    one[0] = 1;
    fp_t t;
    fp_mul(t, a, one);              /* a*R * 1 * R^-1 = a */
    memcpy(plain, t, sizeof(fp_t));
}

int fp_is_zero(const fp_t a)
{
    limb_t acc = 0;
    for (int i = 0; i < NLIMB; i++) acc |= a[i];
    return (int)(1 - (mask_nonzero(acc) & 1));
}

int fp_eq(const fp_t a, const fp_t b)
{
    limb_t acc = 0;
    for (int i = 0; i < NLIMB; i++) acc |= (a[i] ^ b[i]);
    return (int)(1 - (mask_nonzero(acc) & 1));
}

/* Inversion.
 *
 * Measured on BLS12-461, 8 limbs, Apple Silicon:
 *
 *   mpz_invert        (variable time)   1.4 us
 *   mpn_sec_invert    (constant time)  26.4 us
 *   a^(p-2) Fermat    (constant time)  76.3 us
 *
 * So constant-time inversion costs roughly 19x variable-time inversion no
 * matter which of the two we pick, and the Fermat chain the plan suggested
 * starting with is the worst of them. GMP's mpn_sec_invert wins and is a
 * tenth of the code, so use it.
 *
 * The consequence matters more than the numbers: at 26 us an inversion, the
 * affine Miller loop's one-inversion-per-iteration would cost about 2 ms per
 * pairing on its own. Inversions have to leave the loop (projective
 * coordinates, Phase 4) before constant-time inversion is affordable at all.
 *
 * Montgomery bookkeeping: for a_mont = A*R, mpn_sec_invert returns
 * (A*R)^-1 = A^-1 * R^-1, but the Montgomery form of A^-1 is A^-1 * R. Two
 * Montgomery multiplications by R^2 make up the missing factor of R^2, since
 * montmul(montmul(x, R2), R2) = x * R^2.
 */
void fp_inv(fp_t r, const fp_t a)
{
    mp_limb_t scratch[64];                            /* itch is 4n; 64 covers n<=16 */
    fp_t t, out, zero;

    fp_copy(t, a);                                    /* mpn_sec_invert clobbers its input */
    (void)mpn_sec_invert(out, t, FP_MODULUS, FP_LIMBS,
                         2 * FP_LIMBS * GMP_NUMB_BITS, scratch);
    fp_mul(out, out, FP_R2);
    fp_mul(out, out, FP_R2);

    /* Zero has no inverse, and fp_inv(0) is defined to be 0. The obvious way to
     * write that is an early return, which is a branch on the operand: the one
     * thing this routine is not allowed to do. mpn_sec_invert runs on zero
     * without complaint -- it just reports failure and leaves a meaningless
     * result -- so run it unconditionally and select. */
    fp_set_zero(zero);
    fp_cselect(r, zero, out, (limb_t)0 - (limb_t)fp_is_zero(a));
}

/* Variable-time inversion, for values that are already public: a point being
 * deserialised, a verification-only path, or a benchmark. Roughly 19x faster
 * than fp_inv. Never call this on a secret. */
void fp_inv_vartime(fp_t r, const fp_t a)
{
    if (fp_is_zero(a)) { fp_set_zero(r); return; }

    mpz_t A, P, R;
    mpz_inits(A, P, R, NULL);
    mpz_import(A, FP_LIMBS, -1, sizeof(limb_t), 0, 0, a);
    mpz_import(P, FP_LIMBS, -1, sizeof(limb_t), 0, 0, FP_MODULUS);
    mpz_invert(R, A, P);
    fp_set_zero(r);
    mpz_export(r, NULL, -1, sizeof(limb_t), 0, 0, R);
    mpz_clears(A, P, R, NULL);
    fp_mul(r, r, FP_R2);
    fp_mul(r, r, FP_R2);
}

/* --- exponentiation and square roots -------------------------------------- */

void fp_exp(fp_t r, const fp_t a, const limb_t *e, int ebits)
{
    fp_t acc;
    fp_set_one(acc);
    for (int i = ebits - 1; i >= 0; i--) {
        fp_sqr(acc, acc);
        if ((e[i / 64] >> (i % 64)) & 1) fp_mul(acc, acc, a);
    }
    fp_copy(r, acc);
}

/* (p + add_one) >> shift. The three exponents the square roots need all have
 * this shape, because p is odd and congruent to 3 mod 4:
 *
 *   (p+1)/4  = (p+1) >> 2      the Fp square root
 *   (p-1)/2  = p >> 1          the quadratic character, and "lexicographically
 *                              largest" for the compressed encodings
 *   (p-3)/4  = p >> 2          the first step of the Fp2 square root
 *
 * All three are functions of the modulus alone, so this runs on public data.
 */
static void modulus_shifted(limb_t *out, int add_one, int shift)
{
    limb_t t[FP_LIMBS + 1];
    for (int i = 0; i < NLIMB; i++) t[i] = FP_MODULUS[i];
    t[NLIMB] = 0;

    if (add_one) {
        limb_t carry = 1;
        for (int i = 0; i <= NLIMB && carry; i++) {
            t[i] += carry;
            carry = (t[i] == 0);
        }
    }

    for (int s = 0; s < shift; s++) {
        for (int i = 0; i < NLIMB; i++)
            t[i] = (t[i] >> 1) | (t[i + 1] << 63);
        t[NLIMB] >>= 1;
    }
    for (int i = 0; i < NLIMB; i++) out[i] = t[i];
}

void fp_exp_constants(limb_t *sqrt_e, limb_t *half_e, limb_t *quarter_e)
{
    if (sqrt_e)    modulus_shifted(sqrt_e,    1, 2);   /* (p+1)/4 */
    if (half_e)    modulus_shifted(half_e,    0, 1);   /* (p-1)/2 */
    if (quarter_e) modulus_shifted(quarter_e, 0, 2);   /* (p-3)/4 */
}

/* All three supported primes satisfy p = 3 (mod 4). Checked rather than
 * assumed: a curve added later with p = 1 (mod 4) needs Tonelli-Shanks, and
 * silently returning wrong roots is exactly the class of defect this project
 * exists to remove. */
static int modulus_is_3_mod_4(void) { return (FP_MODULUS[0] & 3u) == 3u; }

int fp_sqrt(fp_t r, const fp_t a)
{
    if (!modulus_is_3_mod_4()) { fp_set_zero(r); return 0; }

    limb_t e[FP_LIMBS];
    fp_exp_constants(e, NULL, NULL);

    fp_t c, chk;
    fp_exp(c, a, e, FP_BITS);
    fp_sqr(chk, c);
    int ok = fp_eq(chk, a);

    /* Branch-free: a non-residue yields zero, a residue yields the root. */
    limb_t mask = (limb_t)0 - (limb_t)ok;
    fp_t zero;
    fp_set_zero(zero);
    fp_cselect(r, c, zero, mask);
    return ok;
}

int fp_is_lex_largest(const fp_t a)
{
    /* Compare the canonical residue against (p-1)/2. Larger means the borrow
     * out of half - a is set, i.e. a > (p-1)/2. */
    limb_t plain[FP_LIMBS], half[FP_LIMBS];
    fp_to_limbs(plain, a);
    fp_exp_constants(NULL, half, NULL);   /* (p-1)/2, since p is odd */

    limb_t borrow = 0;
    for (int i = 0; i < NLIMB; i++) {
        limb_t x = half[i], y = plain[i];
        limb_t d  = x - y;
        limb_t b1 = (x < y);
        limb_t d2 = d - borrow;
        limb_t b2 = (d < borrow);
        (void)d2;
        borrow = b1 | b2;
    }
    return (int)borrow;
}
