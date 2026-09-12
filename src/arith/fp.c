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
#include "ct.h"
#include <stdint.h>

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

/* Montgomery product: r = a*b*R^-1 mod p.
 *
 * The portable one. It is not a fallback in the sense of being second best
 * and untested: it is the definition of what fp_mul means, every assembly
 * backend is checked against it, and test/fp_difftest.c runs both through
 * GMP on every build. */
void fp_mul_portable(fp_t r, const fp_t a, const fp_t b)
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

/* ------------------------------------------------- backend selection ----
 *
 * fp_mul is half the time in a BLS12-381 pairing and two thirds of one on
 * BN-462, measured by running a pairing with fp_mul doing its work twice and
 * subtracting. That is what justifies assembly here and nowhere else yet.
 *
 * What the assembly has that C does not is two carry chains. adcx carries
 * through the carry flag and adox through the overflow flag, so the low and
 * high halves of the partial products accumulate at the same time. A C
 * compiler has no way to say that. Fixing the width and unrolling by hand, the
 * best plain C managed, was 1.19x; the assembly is 1.9x to 2.1x. AArch64 has
 * no second carry chain, so it wins differently: mul and umulh leave the flags
 * alone, so the multiplies for the next limbs issue inside the carry chain
 * that consumes the previous ones.
 *
 * The choice is made once, from CPUID, and stored. It depends on the machine
 * and never on an operand, so the branch it costs is outside the constant-time
 * argument. Setting ELIPS_NO_ASM in the environment forces the portable path,
 * which is how CI runs the entire suite both ways.
 */
#if FP_LIMBS == 6 || FP_LIMBS == 8
#  if defined(__x86_64__) && !defined(__ILP32__)
#    define ELIPS_FP_MUL_ASM   ELIPS_FP_MUL_X86_64
#  elif defined(__aarch64__) || defined(__arm64__)
#    define ELIPS_FP_MUL_ASM   ELIPS_FP_MUL_AARCH64
#  endif
#endif

#if defined(ELIPS_FP_MUL_ASM)
#include <stdlib.h>

/* The assembly keeps FP_LIMBS+1 words and never propagates a carry past the
 * top one. That holds exactly when the modulus is sparse. All three curves
 * are, with room to spare, but a future curve might not be. */
_Static_assert(FP_BITS < 64 * FP_LIMBS - 1,
               "the fp_mul assembly needs a sparse modulus, p < 2^(64n-1)");

#define ELIPS_FP_MUL_DECL(fn) \
    void fn(limb_t *r, const limb_t *a, const limb_t *b, \
            const limb_t *p, limb_t n0)

#if ELIPS_FP_MUL_ASM == ELIPS_FP_MUL_X86_64
#include <cpuid.h>
ELIPS_FP_MUL_DECL(elips_fp_mul_6_x86_64);
ELIPS_FP_MUL_DECL(elips_fp_mul_8_x86_64);
#  if FP_LIMBS == 6
#    define ELIPS_FP_MUL_ASM_FN elips_fp_mul_6_x86_64
#  else
#    define ELIPS_FP_MUL_ASM_FN elips_fp_mul_8_x86_64
#  endif

/* mulx needs BMI2 and adcx/adox need ADX, and neither is baseline, so the
 * backend has to be chosen from CPUID. */
static int asm_available(void)
{
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid_max(0, NULL) < 7) return 0;
    if (!__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) return 0;
    /* Leaf 7 subleaf 0: EBX bit 8 is BMI2, bit 19 is ADX. Both are needed;
     * a machine with one and not the other would fault on the other. */
    return (ebx & (1u << 8)) && (ebx & (1u << 19));
}
#else
ELIPS_FP_MUL_DECL(elips_fp_mul_6_aarch64);
ELIPS_FP_MUL_DECL(elips_fp_mul_8_aarch64);
#  if FP_LIMBS == 6
#    define ELIPS_FP_MUL_ASM_FN elips_fp_mul_6_aarch64
#  else
#    define ELIPS_FP_MUL_ASM_FN elips_fp_mul_8_aarch64
#  endif

/* mul, umulh, adcs and csel are all baseline ARMv8-A. Nothing to detect. */
static int asm_available(void) { return 1; }
#endif

static int fp_mul_use_asm;

__attribute__((constructor)) static void fp_mul_select_backend(void)
{
    if (getenv("ELIPS_NO_ASM") != NULL) return;
    fp_mul_use_asm = asm_available();
}

void fp_mul(fp_t r, const fp_t a, const fp_t b)
{
    if (fp_mul_use_asm) {
        ELIPS_FP_MUL_ASM_FN(r, a, b, FP_MODULUS, FP_MONT_N0);
        return;
    }
    fp_mul_portable(r, a, b);
}

int fp_mul_backend(void)
{
    return fp_mul_use_asm ? ELIPS_FP_MUL_ASM : ELIPS_FP_MUL_PORTABLE;
}

#else /* no assembly for this target */

void fp_mul(fp_t r, const fp_t a, const fp_t b) { fp_mul_portable(r, a, b); }

int fp_mul_backend(void) { return ELIPS_FP_MUL_PORTABLE; }

#endif

const char *fp_mul_backend_name(void)
{
    switch (fp_mul_backend()) {
    case ELIPS_FP_MUL_X86_64:  return "x86-64 mulx/adcx/adox";
    case ELIPS_FP_MUL_AARCH64: return "aarch64 mul/umulh";
    default:                   return "portable C";
    }
}

void fp_sqr(fp_t r, const fp_t a)
{
    /* A dedicated squaring saves roughly a third of the partial products, and
     * an earlier note here deferred one to the Phase 6 assembly. Phase 6 then
     * measured what it would be worth and did not write it.
     *
     * Only 4.7% of the fp_mul calls in a BLS12-381 pairing have a == b, and
     * 5.4% on BN-462, because the Fp2, Fp6 and Fp12 layers carry their own
     * squaring formulas and never reach this path with equal operands. A third
     * off 5% of the multiplies is under 1% of a pairing, for four more
     * hand-written routines to keep correct. */
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
/* ------------------------------------------- constant-time inversion ----
 *
 * Batched "divsteps", after Bernstein and Yang (ePrint 2019/266) and Pornin
 * (ePrint 2020/972). The routine this replaces was mpn_sec_invert, which is
 * correct but does one full-width pass per bit: 417 field multiplications'
 * worth on BLS12-381, against 94 ns for a multiply.
 *
 * ONE DIVSTEP, on state (delta, f, g) with f odd:
 *
 *   delta > 0 and g odd :  (1-delta,  g, (g-f)/2)
 *   otherwise           :  (1+delta,  f, (g + (g&1) f)/2)
 *
 * WHY IT BATCHES SAFELY. The branch depends only on delta and the LOW BIT of
 * g. So 62 consecutive divsteps can be run on nothing but delta and the low
 * limbs of f and g, in registers, accumulating a 2x2 integer matrix; the
 * expensive full-width work then happens once per 62 steps instead of once per
 * step. This needs no approximation of the high bits, which is the part of
 * Pornin's variant carrying a correctness obligation. Truncating to 64 bits is
 * exact for the decisions: each step consumes one bit, and 62 < 64.
 *
 * Each divstep is one half of an integer matrix,
 *
 *   case A: (1/2) [[0, 2], [-1, 1]]      case B: (1/2) [[2, 0], [c, 1]]
 *
 * so a block of K steps is (1/2^K) times an integer matrix M, and
 * (f, g) <- M (f, g) >> K exactly.
 *
 * THE INVERSE. Track d, e with f == d*a and g == e*a (mod p), starting from
 * (f,g) = (p,a) and (d,e) = (0,1). The same matrix drives them. When g reaches
 * zero f is +-1, and d is the inverse up to sign and a power of two.
 *
 * d and e get ONE Montgomery step per block, radix 2^64, because that is the
 * reduction this file already has. So they pick up 2^-64 a block while f and g
 * pick up 2^-62. FP_INV_FIX closes that gap and restores the Montgomery factor
 * in a single multiplication at the end.
 *
 * ITERATION COUNT. Fixed at 3*FP_BITS: above the bound in 2019/266 (about
 * 2.88*bits) and 38% above the worst case measured over ~6000 inputs per curve,
 * including Fibonacci pairs and all-ones patterns. tools/reference/divstep_ref.py
 * carries that measurement and says plainly which part is measured and which
 * part is recalled. The count is fixed and public, so it leaks nothing.
 *
 * CONSTANT TIME IS NOT OBVIOUS HERE, AND WAS NOT FREE. Two branches on
 * secret-derived values -- the signs of f and g in mat_shift, and the delta > 0
 * test -- made dudect read |t| = 80 on fp_inv. Both are masks now. That is why
 * test/dudect_test.c exists; reading the code had not found them.
 *
 * The whole algorithm was modelled at limb level in Python and checked against
 * exact inverses before any of this was written. test/edge_test.c checks it
 * against GMP in test/edge_test.c, where GMP is still a dependency. */

#define INV_K       62
#define INV_ITERS   (3 * FP_BITS)
#define INV_BLOCKS  ((INV_ITERS + INV_K - 1) / INV_K)
#define INV_N       (FP_LIMBS + 1)

/* 62 divsteps from delta and the low limbs alone. Branchless: both cases are
 * evaluated and selected under a mask. */
static void divsteps(int64_t *deltap, limb_t f, limb_t g,
                     int64_t *up, int64_t *vp, int64_t *qp, int64_t *rp)
{
    int64_t u = 1, v = 0, q = 0, r = 1, delta = *deltap;

    for (int i = 0; i < INV_K; i++) {
        limb_t podd = (limb_t)0 - (limb_t)(g & 1u);          /* g odd     */
        /* delta evolves from the secret's low bits, so its sign is secret too.
         * An arithmetic shift of (delta-1), not a compare, so there is nothing
         * a compiler could turn back into a branch. */
        limb_t pos  = ~(limb_t)(uint64_t)((int64_t)(delta - 1) >> 63);
        limb_t m    = podd & pos;                            /* case A    */
        int64_t im  = (int64_t)m;

        delta = 1 + ((delta ^ im) - im);                     /* 1 -+ delta */

        limb_t nf = f ^ ((f ^ g) & m);                       /* m ? g : f  */
        limb_t cf = f & podd;                                /* (g&1) * f  */
        limb_t t  = g + ((cf ^ m) - m);                      /* g -+ that  */

        /* case A: (u,v,q,r) <- (2q, 2r, q-u, r-v)
         * case B: (u,v,q,r) <- (2u, 2v, c*u+q, c*v+r)
         * With m all ones, (cu ^ m) - m is -cu = -u, so one form serves both. */
        int64_t nu = (int64_t)((limb_t)u ^ (((limb_t)u ^ (limb_t)q) & m));
        int64_t nv = (int64_t)((limb_t)v ^ (((limb_t)v ^ (limb_t)r) & m));
        limb_t  cu = (limb_t)u & podd;
        limb_t  cv = (limb_t)v & podd;
        int64_t aq = q + (int64_t)((cu ^ m) - m);
        int64_t ar = r + (int64_t)((cv ^ m) - m);

        u = (int64_t)((limb_t)nu << 1);
        v = (int64_t)((limb_t)nv << 1);
        q = aq;
        r = ar;
        f = nf;
        g = t >> 1;
    }
    *deltap = delta; *up = u; *vp = v; *qp = q; *rp = r;
}

/* out = (u*F + v*G) >> INV_K, with F and G two's complement in INV_N limbs. */
static void mat_shift(limb_t *out, int64_t u, const limb_t *F,
                      int64_t v, const limb_t *G)
{
    limb_t t[INV_N + 1];
    __int128 carry = 0;
    for (int i = 0; i < INV_N; i++) {
        carry += (__int128)u * (__int128)(unsigned long long)F[i];
        carry += (__int128)v * (__int128)(unsigned long long)G[i];
        t[i] = (limb_t)carry;
        carry >>= 64;
    }
    /* The loop treated F and G as unsigned digits. Two's complement means the
     * true value is that sum minus 2^(64*INV_N) when the top bit is set, so
     * correct the final carry rather than every limb.
     *
     * By mask, not by "if". The signs here depend on the operand, and when
     * these were two "if" statements dudect read |t| = 80 on fp_inv. */
    limb_t mF = (limb_t)0 - (F[INV_N - 1] >> 63);
    limb_t mG = (limb_t)0 - (G[INV_N - 1] >> 63);
    carry -= (__int128)(int64_t)((limb_t)u & mF);
    carry -= (__int128)(int64_t)((limb_t)v & mG);
    t[INV_N] = (limb_t)carry;

    for (int i = 0; i < INV_N; i++)
        out[i] = (t[i] >> INV_K) | (t[i + 1] << (64 - INV_K));
}

/* out = (u*d + v*e) with one Montgomery step, normalised into [0, p). */
static void redc_lin(fp_t out, int64_t u, const fp_t d, int64_t v, const fp_t e)
{
    limb_t t[FP_LIMBS + 2];
    __int128 carry = 0;
    for (int i = 0; i < FP_LIMBS; i++) {
        carry += (__int128)u * (__int128)(unsigned long long)d[i];
        carry += (__int128)v * (__int128)(unsigned long long)e[i];
        t[i] = (limb_t)carry;
        carry >>= 64;
    }
    t[FP_LIMBS]     = (limb_t)carry;
    t[FP_LIMBS + 1] = (limb_t)(carry >> 64);          /* sign extension */

    limb_t m = (limb_t)(t[0] * FP_MONT_N0);
    limb_t c = 0;
    for (int i = 0; i < FP_LIMBS; i++) {
        /* Unsigned: m and the modulus limb are both full 64-bit, so their
         * product reaches 2^128 and would overflow a SIGNED __int128. */
        unsigned __int128 s = (unsigned __int128)m
                            * (unsigned __int128)FP_MODULUS[i]
                            + (unsigned __int128)t[i]
                            + (unsigned __int128)c;
        t[i] = (limb_t)s;
        c = (limb_t)(s >> 64);
    }
    /* Built unsigned on purpose: this value is genuinely negative sometimes,
     * and shifting a negative signed __int128 left is undefined. Unsigned
     * wraparound is defined and gives the same bits. This was a real bug --
     * correct at -O2, wrong at -O3, until UBSan named the line. */
    unsigned __int128 hi = ((unsigned __int128)t[FP_LIMBS + 1] << 64)
                         | (unsigned __int128)t[FP_LIMBS];
    hi += (unsigned __int128)c;
    for (int i = 0; i < FP_LIMBS - 1; i++) t[i] = t[i + 1];
    t[FP_LIMBS - 1] = (limb_t)hi;
    t[FP_LIMBS]     = (limb_t)(hi >> 64);

    /* |value| < 2p, signed, in FP_LIMBS+1 limbs. Add 2p, then subtract p while
     * it stays at least p. Fixed trip count, conditional subtraction by mask. */
    limb_t acc[FP_LIMBS + 1];
    for (int i = 0; i <= FP_LIMBS; i++) acc[i] = t[i];
    for (int k = 0; k < 2; k++) {
        limb_t cc = 0;
        for (int i = 0; i < FP_LIMBS; i++) {
            limb_t s0 = acc[i] + FP_MODULUS[i];
            limb_t c1 = s0 < acc[i];
            limb_t s1 = s0 + cc;
            acc[i] = s1;
            cc = c1 | (s1 < s0);
        }
        acc[FP_LIMBS] += cc;
    }
    for (int k = 0; k < 4; k++) {
        limb_t tmp[FP_LIMBS + 1], borrow = 0;
        for (int i = 0; i < FP_LIMBS; i++) {
            limb_t s0 = acc[i] - FP_MODULUS[i];
            limb_t b1 = acc[i] < FP_MODULUS[i];
            limb_t s1 = s0 - borrow;
            tmp[i] = s1;
            borrow = b1 | (s0 < borrow);
        }
        tmp[FP_LIMBS] = acc[FP_LIMBS] - borrow;
        /* ct_mask, and this one is issue #32. Without it clang proves keep is
         * 0 or ~0, reads the select below as "take tmp or take acc", and emits
         * `sub` then `js` -- a branch on the borrow out of a subtraction of the
         * secret. Four of them, one per unrolled k. gcc emits cmov and does
         * not, which is why this only ever showed under clang.
         *
         * Same defect as #30, different instructions, and that is exactly why
         * the scan built for #30 did not see it. tools/verify/ct_branch_scan.py
         * now looks for this shape too. */
        limb_t keep = ct_mask((limb_t)0 - (limb_t)((tmp[FP_LIMBS] >> 63) == 0));
        for (int i = 0; i <= FP_LIMBS; i++)
            acc[i] = (tmp[i] & keep) | (acc[i] & ~keep);
    }
    for (int i = 0; i < FP_LIMBS; i++) out[i] = acc[i];
}

void fp_inv(fp_t r, const fp_t a)
{
    limb_t F[INV_N], G[INV_N];
    fp_t d, e, out, zero, nd;
    int64_t delta = 1;

    /* f = p, g = a in plain form: the divsteps are integer arithmetic, not
     * field arithmetic, so the Montgomery factor comes off here and goes back
     * on via FP_INV_FIX at the end. */
    for (int i = 0; i < FP_LIMBS; i++) F[i] = FP_MODULUS[i];
    F[INV_N - 1] = 0;
    fp_to_limbs(G, a);
    G[INV_N - 1] = 0;

    fp_set_zero(d);                 /* f == d*a (mod p), and p == 0 */
    fp_set_zero(e); e[0] = 1;       /* g == e*a (mod p), plain 1    */

    for (int b = 0; b < INV_BLOCKS; b++) {
        int64_t u, v, q, w;
        limb_t Fn[INV_N], Gn[INV_N];
        fp_t dn, en;
        divsteps(&delta, F[0], G[0], &u, &v, &q, &w);
        mat_shift(Fn, u, F, v, G);
        mat_shift(Gn, q, F, w, G);
        redc_lin(dn, u, d, v, e);
        redc_lin(en, q, d, w, e);
        for (int i = 0; i < INV_N; i++) { F[i] = Fn[i]; G[i] = Gn[i]; }
        fp_copy(d, dn); fp_copy(e, en);
    }

    /* f is now +-1; a negative one means d has the wrong sign. Select, because
     * the sign depends on the operand. */
    limb_t neg = (limb_t)0 - (F[INV_N - 1] >> 63);
    fp_neg(nd, d);
    fp_cselect(d, nd, d, neg);

    fp_mul(out, d, FP_INV_FIX);

    /* Zero has no inverse and fp_inv(0) is defined to be 0, selected rather
     * than branched for the same reason fp_mul has no early return. */
    fp_set_zero(zero);
    fp_cselect(r, zero, out, (limb_t)0 - (limb_t)fp_is_zero(a));
}

/* Variable-time inversion, for values that are already public: a point being
 * deserialised, a verification-only path, or a benchmark. Never call this on a
 * secret.
 *
 * Binary extended Euclid, after Knuth 4.5.2 exercise 39. It branches on the
 * operand at every step and that is the POINT: test/dudect_test.c uses this as
 * its negative control, the routine that proves the timing harness can detect
 * a leak at all. A "faster, constant-time" rewrite here would silently turn
 * that control into one that always passes, which is worse than having none.
 *
 * It replaced a call to GMP's mpz_invert when GMP became a test-only
 * dependency. test/edge_test.c checks it against GMP, where GMP still lives.
 */

/* x <- x/2 mod p, for x < p. p is odd, so an odd x is made even by adding p
 * first; that sum can exceed FP_LIMBS limbs, hence the carry. */
static void half_mod_p(fp_t x)
{
    limb_t carry = 0;
    if (x[0] & 1) {
        for (int i = 0; i < NLIMB; i++) {
            limb_t s  = x[i] + FP_MODULUS[i];
            limb_t c1 = (s < x[i]);
            limb_t s2 = s + carry;
            limb_t c2 = (s2 < s);
            x[i] = s2;
            carry = c1 | c2;
        }
    }
    for (int i = 0; i < NLIMB - 1; i++)
        x[i] = (x[i] >> 1) | (x[i + 1] << 63);
    x[NLIMB - 1] = (x[NLIMB - 1] >> 1) | (carry << 63);
}

static int is_one(const limb_t *x)
{
    if (x[0] != 1) return 0;
    for (int i = 1; i < NLIMB; i++) if (x[i]) return 0;
    return 1;
}

/* a >= b, as plain integers. */
static int geq(const limb_t *a, const limb_t *b)
{
    for (int i = NLIMB - 1; i >= 0; i--) {
        if (a[i] != b[i]) return a[i] > b[i];
    }
    return 1;
}

/* a <- a - b, as plain integers, where a >= b. */
static void sub_plain(limb_t *a, const limb_t *b)
{
    limb_t borrow = 0;
    for (int i = 0; i < NLIMB; i++) {
        limb_t ai = a[i], bi = b[i];
        limb_t d  = ai - bi;
        limb_t b1 = (ai < bi);
        limb_t d2 = d - borrow;
        limb_t b2 = (d < borrow);
        a[i] = d2;
        borrow = b1 | b2;
    }
}

void fp_inv_vartime(fp_t r, const fp_t a)
{
    if (fp_is_zero(a)) { fp_set_zero(r); return; }

    fp_t u, v, x1, x2;
    fp_copy(u, a);
    memcpy(v, FP_MODULUS, sizeof(fp_t));
    fp_set_zero(x1); x1[0] = 1;
    fp_set_zero(x2);

    while (!is_one(u) && !is_one(v)) {
        while ((u[0] & 1) == 0) {
            for (int i = 0; i < NLIMB - 1; i++) u[i] = (u[i] >> 1) | (u[i + 1] << 63);
            u[NLIMB - 1] >>= 1;
            half_mod_p(x1);
        }
        while ((v[0] & 1) == 0) {
            for (int i = 0; i < NLIMB - 1; i++) v[i] = (v[i] >> 1) | (v[i + 1] << 63);
            v[NLIMB - 1] >>= 1;
            half_mod_p(x2);
        }
        if (geq(u, v)) { sub_plain(u, v); fp_sub(x1, x1, x2); }
        else           { sub_plain(v, u); fp_sub(x2, x2, x1); }
    }

    /* The inverse here is of a read as a plain integer, so it is
     * (A*R)^-1 = A^-1 * R^-1. The Montgomery form wanted is A^-1 * R, so two
     * multiplications by R^2 supply the missing R^2. */
    fp_copy(r, is_one(u) ? x1 : x2);
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

int fp_sgn0(const fp_t a)
{
    limb_t plain[FP_LIMBS];
    fp_to_limbs(plain, a);
    return (int)(plain[0] & 1u);
}
