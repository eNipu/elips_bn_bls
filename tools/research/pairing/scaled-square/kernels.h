/* Research only. Scaled squarings do not preserve a raw Miller API.
 * Fp4 = Fp2[s]/(s^2-xi), Fp12 = Fp4[w]/(w^3-s).
 */
#include "elips/pairing.h"

typedef fp2_t quartic[2];

static void qadd(quartic r, const quartic a, const quartic b)
{ fp2_add(r[0], a[0], b[0]); fp2_add(r[1], a[1], b[1]); }
static void qsub(quartic r, const quartic a, const quartic b)
{ fp2_sub(r[0], a[0], b[0]); fp2_sub(r[1], a[1], b[1]); }
static void qdouble(quartic r, const quartic a) { qadd(r, a, a); }
static void qfour(quartic r, const quartic a)
{ qdouble(r, a); qdouble(r, r); }
static void qmul_s(quartic r, const quartic a)
{
    fp2_t t;
    fp2_mul_xi(t, a[1]);
    fp2_copy(r[1], a[0]); fp2_copy(r[0], t);
}
static void qmul_u(quartic r, const quartic a)
{
    fp_t t;
    for (int j = 0; j < 2; j++) {
        fp_neg(t, a[j][1]);
        fp_copy(r[j][1], a[j][0]); fp_copy(r[j][0], t);
    }
}
static void qmul_minus_u(quartic r, const quartic a)
{
    fp_t t;
    for (int j = 0; j < 2; j++) {
        fp_neg(t, a[j][0]);
        fp_copy(r[j][0], a[j][1]); fp_copy(r[j][1], t);
    }
}
static void qsqr(quartic r, const quartic a)
{
    fp2_t aa, bb, t;
    fp2_sqr(aa, a[0]); fp2_sqr(bb, a[1]);
    fp2_add(t, a[0], a[1]); fp2_sqr(t, t);
    fp2_sub(t, t, aa); fp2_sub(t, t, bb);
    fp2_mul_xi(bb, bb); fp2_add(r[0], aa, bb); fp2_copy(r[1], t);
}
static void qmul(quartic r, const quartic a, const quartic b)
{
    fp2_t aa, bb, t, u;
    fp2_mul(aa, a[0], b[0]); fp2_mul(bb, a[1], b[1]);
    fp2_add(t, a[0], a[1]); fp2_add(u, b[0], b[1]); fp2_mul(t, t, u);
    fp2_sub(t, t, aa); fp2_sub(t, t, bb);
    fp2_mul_xi(bb, bb); fp2_add(r[0], aa, bb); fp2_copy(r[1], t);
}
static void qsplit(quartic a, quartic b, quartic c, const fp12_t f)
{
    fp2_copy(a[0], f[0][0]); fp2_copy(a[1], f[1][1]);
    fp2_copy(b[0], f[1][0]); fp2_copy(b[1], f[0][2]);
    fp2_copy(c[0], f[0][1]); fp2_copy(c[1], f[1][2]);
}
static void qjoin(fp12_t f, const quartic a, const quartic b, const quartic c)
{
    fp2_copy(f[0][0], a[0]); fp2_copy(f[1][1], a[1]);
    fp2_copy(f[1][0], b[0]); fp2_copy(f[0][2], b[1]);
    fp2_copy(f[0][1], c[0]); fp2_copy(f[1][2], c[1]);
}

/* 4*f^2 at a cost of 5 S4 = 15 S2 = 30 Mp, plus linear operations.
 * Evaluate at 0, 1, -1, u, infinity. No division or nontrivial constants.
 */
static void gaussian_sqr(fp12_t r, const fp12_t f)
{
    quartic a, b, c, p0, pi, pp, pm, pu, ac, S, D, T, J, t;
    quartic h0, h1, h2;
    qsplit(a, b, c, f);
    qsqr(p0, a); qsqr(pi, c);
    qadd(ac, a, c);
    qadd(t, ac, b); qsqr(pp, t);
    qsub(t, ac, b); qsqr(pm, t);
    qsub(t, a, c); qmul_u(pu, b); qadd(t, t, pu); qsqr(pu, t);
    qadd(S, pp, pm); qsub(D, pp, pm);
    qadd(T, p0, pi); qfour(T, T);
    qdouble(t, pu); qadd(t, t, S); qsub(t, t, T); qmul_minus_u(J, t);
    qsub(t, D, J); qmul_s(t, t); qfour(h0, p0); qadd(h0, h0, t);
    qfour(t, pi); qmul_s(t, t); qadd(h1, D, J); qadd(h1, h1, t);
    qdouble(h2, S); qsub(h2, h2, T);
    qjoin(r, h0, h1, h2);
}

/* 4*f^2, also 5 S4, using 1,-1,u,-u,infinity. The four finite points
 * form a length-four Fourier transform; interpolating needs fewer adds.
 */
static void fft4_sqr(fp12_t r, const fp12_t f)
{
    quartic a, b, c, ac, amc, ub, pp, pm, pu, mu, c4, S, T, U, V, t;
    quartic h0, h1, h2;
    qsplit(a, b, c, f);
    qadd(ac, a, c); qsub(amc, a, c); qmul_u(ub, b);
    qadd(t, ac, b); qsqr(pp, t);
    qsub(t, ac, b); qsqr(pm, t);
    qadd(t, amc, ub); qsqr(pu, t);
    qsub(t, amc, ub); qsqr(mu, t);
    qsqr(c4, c); qfour(c4, c4);
    qadd(S, pp, pm); qadd(T, pu, mu); qsub(U, pp, pm);
    qsub(t, pu, mu); qmul_minus_u(V, t);
    qadd(h0, S, T); qsub(h0, h0, c4);
    qsub(t, U, V); qmul_s(t, t); qadd(h0, h0, t);
    qadd(h1, U, V); qmul_s(t, c4); qadd(h1, h1, t);
    qsub(h2, S, T);
    qjoin(r, h0, h1, h2);
}

/* Stronger control than the repository alone: scaled Chung-Hasan-style
 * cubic squaring, 2*f^2 using 4 S4 + 1 M4 = 33 Mp.
 * This is a standard identity with its halving deferred, not new mathematics.
 */
static void cubic_sqr(fp12_t r, const fp12_t f)
{
    quartic a, b, c, p0, pi, pp, pm, bc4, t, ac, h0, h1, h2;
    qsplit(a, b, c, f);
    qsqr(p0, a); qsqr(pi, c);
    qadd(ac, a, c);
    qadd(t, ac, b); qsqr(pp, t);
    qsub(t, ac, b); qsqr(pm, t);
    qmul(bc4, b, c); qfour(bc4, bc4);
    qdouble(h0, p0); qmul_s(t, bc4); qadd(h0, h0, t);
    qsub(h1, pp, pm); qsub(h1, h1, bc4);
    qdouble(t, pi); qmul_s(t, t); qadd(h1, h1, t);
    qadd(h2, pp, pm); qadd(t, p0, pi); qdouble(t, t); qsub(h2, h2, t);
    qjoin(r, h0, h1, h2);
}

/* Existing normalized-line research control, 2*a*(B+C*v) in 4 M2. */
static void linear4(fp6_t r, const fp6_t a, const fp2_t B,
                    const fp2_t C, const fp2_t plus, const fp2_t minus)
{
    fp2_t t0, t3, tp, tm, s, u, q0, q1, q2;
    fp2_mul(t0, a[0], B); fp2_mul(t3, a[2], C);
    fp2_add(s, a[0], a[2]);
    fp2_add(u, s, a[1]); fp2_mul(tp, u, plus);
    fp2_sub(u, s, a[1]); fp2_mul(tm, u, minus);
    fp2_mul_xi(u, t3); fp2_add(q0, t0, u); fp2_add(q0, q0, q0);
    fp2_add(u, t3, t3); fp2_sub(q1, tp, tm); fp2_sub(q1, q1, u);
    fp2_add(u, t0, t0); fp2_add(q2, tp, tm); fp2_sub(q2, q2, u);
    fp2_copy(r[0], q0); fp2_copy(r[1], q1); fp2_copy(r[2], q2);
}
static void normalized8(fp12_t f, const fp2_t B, const fp2_t C)
{
    fp2_t plus, minus;
    fp6_t hA, hD, A, D;
    fp2_add(plus, B, C); fp2_sub(minus, B, C);
    linear4(hA, f[0], B, C, plus, minus);
    linear4(hD, f[1], B, C, plus, minus);
    fp6_add(A, f[0], f[0]); fp6_add(D, f[1], f[1]);
    fp6_mul_v(hD, hD); fp6_mul_v(hD, hD); fp6_mul_v(hA, hA);
    fp6_add(f[0], A, hD); fp6_add(f[1], D, hA);
}
