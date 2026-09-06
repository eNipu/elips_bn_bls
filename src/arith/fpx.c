/*
 * Extension field tower. See include/elips/fpx.h.
 *
 * Karatsuba throughout, matching the formulas the Phase 0 reference uses, so
 * the committed known-answer vectors validate this code directly.
 *
 * ponytail: no lazy reduction yet. Proper lazy reduction needs double-width
 * accumulators threaded through every routine, and the measured projection said
 * plain Karatsuba over the new fp_t already gets close to the gate. Add it only
 * if the measurement at the end of the phase falls short -- the ceiling is one
 * Montgomery reduction per fp_t operation instead of one per fp2/fp6 operation.
 */
#include "elips/fpx.h"

/* ============================== fp2 ==================================== */

void fp2_set_zero(fp2_t r) { fp_set_zero(r[0]); fp_set_zero(r[1]); }
void fp2_set_one(fp2_t r)  { fp_set_one(r[0]);  fp_set_zero(r[1]); }
void fp2_copy(fp2_t r, const fp2_t a) { fp_copy(r[0], a[0]); fp_copy(r[1], a[1]); }
void fp2_add(fp2_t r, const fp2_t a, const fp2_t b)
{ fp_add(r[0], a[0], b[0]); fp_add(r[1], a[1], b[1]); }
void fp2_sub(fp2_t r, const fp2_t a, const fp2_t b)
{ fp_sub(r[0], a[0], b[0]); fp_sub(r[1], a[1], b[1]); }
void fp2_neg(fp2_t r, const fp2_t a) { fp_neg(r[0], a[0]); fp_neg(r[1], a[1]); }
void fp2_conj(fp2_t r, const fp2_t a) { fp_copy(r[0], a[0]); fp_neg(r[1], a[1]); }
void fp2_mul_fp(fp2_t r, const fp2_t a, const fp_t b)
{ fp_mul(r[0], a[0], b); fp_mul(r[1], a[1], b); }
void fp2_cselect(fp2_t r, const fp2_t a, const fp2_t b, limb_t mask)
{ fp_cselect(r[0], a[0], b[0], mask); fp_cselect(r[1], a[1], b[1], mask); }
int fp2_is_zero(const fp2_t a) { return fp_is_zero(a[0]) & fp_is_zero(a[1]); }
int fp2_eq(const fp2_t a, const fp2_t b) { return fp_eq(a[0], b[0]) & fp_eq(a[1], b[1]); }

void fp2_mul(fp2_t r, const fp2_t a, const fp2_t b)
{
    /* (a0 + a1 u)(b0 + b1 u) = (a0b0 - a1b1) + ((a0+a1)(b0+b1) - a0b0 - a1b1) u,
     * using u^2 = -1. Three multiplies rather than four. */
    fp_t t0, t1, s, t;
    fp_mul(t0, a[0], b[0]);
    fp_mul(t1, a[1], b[1]);
    fp_add(s, a[0], a[1]);
    fp_add(t, b[0], b[1]);
    fp_mul(s, s, t);
    fp_sub(s, s, t0);
    fp_sub(s, s, t1);
    fp_sub(r[0], t0, t1);
    fp_copy(r[1], s);
}

void fp2_sqr(fp2_t r, const fp2_t a)
{
    /* (a0 + a1 u)^2 = (a0+a1)(a0-a1) + 2 a0 a1 u */
    fp_t s, d, m;
    fp_add(s, a[0], a[1]);
    fp_sub(d, a[0], a[1]);
    fp_mul(m, a[0], a[1]);
    fp_mul(r[0], s, d);
    fp_add(r[1], m, m);
}

void fp2_mul_xi(fp2_t r, const fp2_t a)
{
    /* multiply by xi = 1 + u: (a0 + a1 u)(1 + u) = (a0 - a1) + (a0 + a1) u */
    fp_t s, d;
    fp_sub(d, a[0], a[1]);
    fp_add(s, a[0], a[1]);
    fp_copy(r[0], d);
    fp_copy(r[1], s);
}

void fp2_inv(fp2_t r, const fp2_t a)
{
    /* conjugate over the norm; the norm a0^2 + a1^2 lives in Fp */
    fp_t n, t;
    fp_sqr(n, a[0]);
    fp_sqr(t, a[1]);
    fp_add(n, n, t);
    fp_inv(n, n);
    fp_mul(r[0], a[0], n);
    fp_mul(t, a[1], n);
    fp_neg(r[1], t);
}

/* ============================== fp6 ==================================== */

void fp6_set_zero(fp6_t r) { fp2_set_zero(r[0]); fp2_set_zero(r[1]); fp2_set_zero(r[2]); }
void fp6_set_one(fp6_t r)  { fp2_set_one(r[0]);  fp2_set_zero(r[1]); fp2_set_zero(r[2]); }
void fp6_copy(fp6_t r, const fp6_t a)
{ fp2_copy(r[0], a[0]); fp2_copy(r[1], a[1]); fp2_copy(r[2], a[2]); }
void fp6_add(fp6_t r, const fp6_t a, const fp6_t b)
{ fp2_add(r[0],a[0],b[0]); fp2_add(r[1],a[1],b[1]); fp2_add(r[2],a[2],b[2]); }
void fp6_sub(fp6_t r, const fp6_t a, const fp6_t b)
{ fp2_sub(r[0],a[0],b[0]); fp2_sub(r[1],a[1],b[1]); fp2_sub(r[2],a[2],b[2]); }
void fp6_neg(fp6_t r, const fp6_t a)
{ fp2_neg(r[0],a[0]); fp2_neg(r[1],a[1]); fp2_neg(r[2],a[2]); }
int fp6_is_zero(const fp6_t a)
{ return fp2_is_zero(a[0]) & fp2_is_zero(a[1]) & fp2_is_zero(a[2]); }
int fp6_eq(const fp6_t a, const fp6_t b)
{ return fp2_eq(a[0],b[0]) & fp2_eq(a[1],b[1]) & fp2_eq(a[2],b[2]); }

void fp6_mul_v(fp6_t r, const fp6_t a)
{
    /* v*(c0 + c1 v + c2 v^2) = c2*xi + c0 v + c1 v^2, since v^3 = xi */
    fp2_t t;
    fp2_mul_xi(t, a[2]);
    fp2_t c0, c1;
    fp2_copy(c0, a[0]);
    fp2_copy(c1, a[1]);
    fp2_copy(r[0], t);
    fp2_copy(r[1], c0);
    fp2_copy(r[2], c1);
}

void fp6_mul(fp6_t r, const fp6_t a, const fp6_t b)
{
    fp2_t t0, t1, t2, s, t, e0, e1, e2;
    fp2_mul(t0, a[0], b[0]);
    fp2_mul(t1, a[1], b[1]);
    fp2_mul(t2, a[2], b[2]);

    fp2_add(s, a[1], a[2]); fp2_add(t, b[1], b[2]); fp2_mul(s, s, t);
    fp2_sub(s, s, t1); fp2_sub(s, s, t2); fp2_mul_xi(s, s); fp2_add(e0, s, t0);

    fp2_add(s, a[0], a[1]); fp2_add(t, b[0], b[1]); fp2_mul(s, s, t);
    fp2_sub(s, s, t0); fp2_sub(s, s, t1); fp2_mul_xi(t, t2); fp2_add(e1, s, t);

    fp2_add(s, a[0], a[2]); fp2_add(t, b[0], b[2]); fp2_mul(s, s, t);
    fp2_sub(s, s, t0); fp2_sub(s, s, t2); fp2_add(e2, s, t1);

    fp2_copy(r[0], e0); fp2_copy(r[1], e1); fp2_copy(r[2], e2);
}

void fp6_sqr(fp6_t r, const fp6_t a)
{
    /* (a0 + a1 v + a2 v^2)^2 with v^3 = xi expands to
     *   c0 = a0^2 + 2*xi*a1*a2
     *   c1 = 2*a0*a1 + xi*a2^2
     *   c2 = a1^2 + 2*a0*a2
     * Each doubled cross term comes from a squaring rather than a
     * multiplication, via 2ab = (a+b)^2 - a^2 - b^2, so this costs six fp2
     * squarings instead of six fp2 multiplications. */
    fp2_t s0, s1, s2, t01, t12, t02, c0, c1, c2, t;

    fp2_sqr(s0, a[0]);
    fp2_sqr(s1, a[1]);
    fp2_sqr(s2, a[2]);

    fp2_add(t, a[0], a[1]); fp2_sqr(t01, t);
    fp2_sub(t01, t01, s0); fp2_sub(t01, t01, s1);   /* 2 a0 a1 */
    fp2_add(t, a[1], a[2]); fp2_sqr(t12, t);
    fp2_sub(t12, t12, s1); fp2_sub(t12, t12, s2);   /* 2 a1 a2 */
    fp2_add(t, a[0], a[2]); fp2_sqr(t02, t);
    fp2_sub(t02, t02, s0); fp2_sub(t02, t02, s2);   /* 2 a0 a2 */

    fp2_mul_xi(t, t12);  fp2_add(c0, s0, t);
    fp2_mul_xi(t, s2);   fp2_add(c1, t01, t);
    fp2_add(c2, s1, t02);

    fp2_copy(r[0], c0); fp2_copy(r[1], c1); fp2_copy(r[2], c2);
}

void fp6_inv(fp6_t r, const fp6_t a)
{
    /* A = c0^2 - xi c1 c2, B = xi c2^2 - c0 c1, C = c1^2 - c0 c2,
     * F = xi (c2 B + c1 C) + c0 A, inverse = (A,B,C) / F */
    fp2_t A, B, C, F, t;
    fp2_sqr(A, a[0]); fp2_mul(t, a[1], a[2]); fp2_mul_xi(t, t); fp2_sub(A, A, t);
    fp2_sqr(B, a[2]); fp2_mul_xi(B, B); fp2_mul(t, a[0], a[1]); fp2_sub(B, B, t);
    fp2_sqr(C, a[1]); fp2_mul(t, a[0], a[2]); fp2_sub(C, C, t);

    fp2_mul(F, a[2], B);
    fp2_mul(t, a[1], C);
    fp2_add(F, F, t);
    fp2_mul_xi(F, F);
    fp2_mul(t, a[0], A);
    fp2_add(F, F, t);
    fp2_inv(F, F);

    fp2_mul(r[0], A, F);
    fp2_mul(r[1], B, F);
    fp2_mul(r[2], C, F);
}

/* ============================== fp12 =================================== */

void fp12_set_zero(fp12_t r) { fp6_set_zero(r[0]); fp6_set_zero(r[1]); }
void fp12_set_one(fp12_t r)  { fp6_set_one(r[0]);  fp6_set_zero(r[1]); }
void fp12_copy(fp12_t r, const fp12_t a) { fp6_copy(r[0], a[0]); fp6_copy(r[1], a[1]); }
void fp12_add(fp12_t r, const fp12_t a, const fp12_t b)
{ fp6_add(r[0],a[0],b[0]); fp6_add(r[1],a[1],b[1]); }
void fp12_sub(fp12_t r, const fp12_t a, const fp12_t b)
{ fp6_sub(r[0],a[0],b[0]); fp6_sub(r[1],a[1],b[1]); }
void fp12_neg(fp12_t r, const fp12_t a) { fp6_neg(r[0],a[0]); fp6_neg(r[1],a[1]); }
void fp12_conj(fp12_t r, const fp12_t a) { fp6_copy(r[0],a[0]); fp6_neg(r[1],a[1]); }
int fp12_is_zero(const fp12_t a) { return fp6_is_zero(a[0]) & fp6_is_zero(a[1]); }
int fp12_eq(const fp12_t a, const fp12_t b)
{ return fp6_eq(a[0],b[0]) & fp6_eq(a[1],b[1]); }

void fp12_mul(fp12_t r, const fp12_t a, const fp12_t b)
{
    /* (d0 + d1 w)(e0 + e1 w) = (d0e0 + v d1e1)
     *                        + ((d0+d1)(e0+e1) - d0e0 - d1e1) w,  w^2 = v */
    fp6_t t0, t1, s, t;
    fp6_mul(t0, a[0], b[0]);
    fp6_mul(t1, a[1], b[1]);
    fp6_add(s, a[0], a[1]);
    fp6_add(t, b[0], b[1]);
    fp6_mul(s, s, t);
    fp6_sub(s, s, t0);
    fp6_sub(s, s, t1);
    fp6_mul_v(t, t1);
    fp6_add(r[0], t0, t);
    fp6_copy(r[1], s);
}

void fp12_sqr(fp12_t r, const fp12_t a)
{
    /* Complex squaring: two fp6 multiplications rather than three.
     *   (d0 + d1 w)^2 = (d0^2 + v d1^2) + 2 d0 d1 w
     * and (d0+d1)(d0 + v d1) - d0d1 - v d0d1 = d0^2 + v d1^2, so the whole
     * thing costs one product of sums plus one d0*d1. */
    fp6_t s, u, m, t;
    fp6_add(s, a[0], a[1]);
    fp6_mul_v(u, a[1]);
    fp6_add(u, u, a[0]);
    fp6_mul(m, a[0], a[1]);
    fp6_mul(s, s, u);
    fp6_sub(s, s, m);
    fp6_mul_v(t, m);
    fp6_sub(s, s, t);
    fp6_add(r[1], m, m);
    fp6_copy(r[0], s);
}

/* Frobenius.
 *
 * In fp12 = fp6[w] with w^2 = v and v^3 = xi, an element is
 *     c0 + c3 w + c1 w^2 + c4 w^3 + c2 w^4 + c5 w^5
 * where the storage is d0 = (c0,c1,c2) and d1 = (c3,c4,c5). Raising to the
 * p^k power conjugates each fp2 coefficient (for odd k) and multiplies the
 * coefficient of w^i by gamma^i, with gamma = w^(p^k - 1) = xi^((p^k-1)/6).
 * The gamma powers are generated into fp_params.h already in Montgomery form,
 * so this needs no setup and no global state. */
void fp12_frobenius(fp12_t r, const fp12_t a, int k)
{
    const limb_t (*g[5])[FP_LIMBS];
    int odd;
    switch (k) {
        case 1: g[0]=FROB_P1_1; g[1]=FROB_P1_2; g[2]=FROB_P1_3;
                g[3]=FROB_P1_4; g[4]=FROB_P1_5; odd = 1; break;
        case 2: g[0]=FROB_P2_1; g[1]=FROB_P2_2; g[2]=FROB_P2_3;
                g[3]=FROB_P2_4; g[4]=FROB_P2_5; odd = 0; break;
        case 3: g[0]=FROB_P3_1; g[1]=FROB_P3_2; g[2]=FROB_P3_3;
                g[3]=FROB_P3_4; g[4]=FROB_P3_5; odd = 1; break;
        case 6: fp12_conj(r, a); return;
        default: fp12_copy(r, a); return;      /* k = 0 or 12 */
    }

    /* coefficient of w^i, in storage order */
    fp12_t out;
    fp2_copy(out[0][0], a[0][0]);              /* w^0, gamma^0 = 1 */
    fp2_copy(out[1][0], a[1][0]);              /* w^1 */
    fp2_copy(out[0][1], a[0][1]);              /* w^2 */
    fp2_copy(out[1][1], a[1][1]);              /* w^3 */
    fp2_copy(out[0][2], a[0][2]);              /* w^4 */
    fp2_copy(out[1][2], a[1][2]);              /* w^5 */

    if (odd) {
        fp2_conj(out[0][0], out[0][0]); fp2_conj(out[1][0], out[1][0]);
        fp2_conj(out[0][1], out[0][1]); fp2_conj(out[1][1], out[1][1]);
        fp2_conj(out[0][2], out[0][2]); fp2_conj(out[1][2], out[1][2]);
    }
    fp2_mul(out[1][0], out[1][0], g[0]);       /* w^1 * gamma^1 */
    fp2_mul(out[0][1], out[0][1], g[1]);       /* w^2 * gamma^2 */
    fp2_mul(out[1][1], out[1][1], g[2]);       /* w^3 * gamma^3 */
    fp2_mul(out[0][2], out[0][2], g[3]);       /* w^4 * gamma^4 */
    fp2_mul(out[1][2], out[1][2], g[4]);       /* w^5 * gamma^5 */
    fp12_copy(r, out);
}

void fp12_sqr_cyc(fp12_t r, const fp12_t a)
{
    /* An element of the cyclotomic subgroup satisfies conj(a) = a^-1, so
     *     (d0 + d1 w)(d0 - d1 w) = d0^2 - v*d1^2 = 1.
     * Squaring normally needs d0^2 + v*d1^2, and the constraint rewrites that
     * as 2*d0^2 - 1, so one squaring and one multiplication suffice where the
     * generic routine needs two multiplications:
     *     a^2 = (2*d0^2 - 1) + (2*d0*d1) w
     * Only valid inside the subgroup; fp12_sqr remains for everything else. */
    fp6_t s, m, one;
    fp6_sqr(s, a[0]);
    fp6_mul(m, a[0], a[1]);
    fp6_add(s, s, s);
    fp6_set_one(one);
    fp6_sub(s, s, one);
    fp6_add(r[1], m, m);
    fp6_copy(r[0], s);
}

void fp12_exp(fp12_t r, const fp12_t a, const limb_t *e, int ebits)
{
    fp12_t acc;
    fp12_set_one(acc);
    for (int i = ebits - 1; i >= 0; i--) {
        fp12_sqr(acc, acc);
        if ((e[i / 64] >> (i % 64)) & 1) fp12_mul(acc, acc, a);
    }
    fp12_copy(r, acc);
}

void fp12_inv(fp12_t r, const fp12_t a)
{
    /* (d0 + d1 w)^-1 = (d0 - d1 w) / (d0^2 - v d1^2) */
    fp6_t f, t;
    fp6_sqr(f, a[0]);
    fp6_sqr(t, a[1]);
    fp6_mul_v(t, t);
    fp6_sub(f, f, t);
    fp6_inv(f, f);
    fp6_mul(r[0], a[0], f);
    fp6_mul(t, a[1], f);
    fp6_neg(r[1], t);
}

/* --- fp2 exponentiation and square roots ----------------------------------
 *
 * Needed by point decompression: a compressed point carries x and one bit of y,
 * and recovering y means taking a square root of x^3 + b in Fp (for G1) or in
 * Fp2 (for G2).
 */

void fp2_mul_u(fp2_t r, const fp2_t a)
{
    /* (c0 + c1 u) * u = -c1 + c0 u, since u^2 = -1. */
    fp_t t;
    fp_copy(t, a[0]);
    fp_neg(r[0], a[1]);
    fp_copy(r[1], t);
}

void fp2_exp(fp2_t r, const fp2_t a, const limb_t *e, int ebits)
{
    fp2_t acc;
    fp2_set_one(acc);
    for (int i = ebits - 1; i >= 0; i--) {
        fp2_sqr(acc, acc);
        if ((e[i / 64] >> (i % 64)) & 1) fp2_mul(acc, acc, a);
    }
    fp2_copy(r, acc);
}

/* Adj and Rodriguez-Henriquez, "Square root computation over even extension
 * fields", Algorithm 9, for q = p = 3 (mod 4):
 *
 *   a1    = a^((p-3)/4)
 *   alpha = a1^2 * a
 *   x0    = a1 * a
 *   x     = i * x0                    if alpha == -1
 *         = (1 + alpha)^((p-1)/2) * x0  otherwise
 *
 * The alpha == -1 case is genuinely value-dependent, so it is resolved with a
 * masked select over both candidates rather than a branch. Both candidates are
 * computed either way; that costs one extra exponentiation and buys a routine
 * whose timing does not reveal which case an input fell into.
 *
 * Rather than Algorithm 9's separate norm test for non-residues, the result is
 * squared and compared. Same decision, one fewer exponentiation, and it also
 * catches any arithmetic slip in the chain above it.
 */
int fp2_sqrt(fp2_t r, const fp2_t a)
{
    limb_t e_quarter[FP_LIMBS], e_half[FP_LIMBS];
    fp_exp_constants(NULL, e_half, e_quarter);

    fp2_t a1, alpha, x0, one, negone, cand_i, cand_b, b, chk, zero;

    fp2_exp(a1, a, e_quarter, FP_BITS);       /* a^((p-3)/4) */
    fp2_sqr(alpha, a1);
    fp2_mul(alpha, alpha, a);                 /* a^((p-1)/2) */
    fp2_mul(x0, a1, a);                       /* a^((p+1)/4) */

    fp2_set_one(one);
    fp2_neg(negone, one);

    fp2_mul_u(cand_i, x0);                    /* the alpha == -1 branch */

    fp2_add(b, one, alpha);
    fp2_exp(b, b, e_half, FP_BITS);
    fp2_mul(cand_b, b, x0);                   /* the ordinary branch */

    limb_t mask = (limb_t)0 - (limb_t)fp2_eq(alpha, negone);
    fp2_t root;
    fp2_cselect(root, cand_i, cand_b, mask);

    /* Everything above accumulates into locals, and r is written only once, at
     * the end. Writing r earlier and then verifying against a would be correct
     * only while the two do not alias -- and every routine in this header
     * promises they may. The earlier version did exactly that: it selected into
     * r, then compared r^2 against a, which by then was the same storage. No
     * caller in the library aliases these, which is why nothing caught it until
     * the aliasing cases in test/edge_test.c went looking. */
    fp2_sqr(chk, root);
    int ok = fp2_eq(chk, a);
    fp2_set_zero(zero);
    fp2_cselect(r, root, zero, (limb_t)0 - (limb_t)ok);
    return ok;
}

int fp2_is_lex_largest(const fp2_t a)
{
    /* Order by the imaginary part first, then the real part -- the rule the
     * BLS12-381 compressed encodings use. Both halves are evaluated and
     * combined with masks, so nothing branches on a. */
    int c1_big  = fp_is_lex_largest(a[1]);
    int c1_zero = fp_is_zero(a[1]);
    int c0_big  = fp_is_lex_largest(a[0]);
    return c1_big | (c1_zero & c0_big);
}

int fp2_sgn0(const fp2_t a)
{
    /* Both halves are evaluated and combined with masks; nothing branches. */
    int s0 = fp_sgn0(a[0]);
    int z0 = fp_is_zero(a[0]);
    int s1 = fp_sgn0(a[1]);
    return s0 | (z0 & s1);
}

int fp2_is_square(const fp2_t a)
{
    /* norm(a0 + a1 u) = a0^2 + a1^2, and a is a square in Fp2 exactly when its
     * norm is a square in Fp. fp_sqrt already reports that. */
    fp_t n, t, root;
    fp_sqr(n, a[0]);
    fp_sqr(t, a[1]);
    fp_add(n, n, t);
    int ok = fp_sqrt(root, n);
    return ok | fp_is_zero(n);
}
