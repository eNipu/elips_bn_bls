/*
 * Optimal ate Miller loop over Jacobian coordinates.
 *
 * Line functions
 * --------------
 * With the D-type twist the library uses (x = x'/w^2, y = y'/w^3, and
 * xi = w^6), the line through the untwisted T evaluated at P = (xP, yP) is
 *
 *     l(P) = yP - lambda*xP*w^-1 + (lambda*x_T - y_T)*w^-3
 *
 * and since w^-1 = w^5/xi and w^-3 = w^3/xi, that is
 *
 *     l(P) = yP + [(lambda*x_T - y_T)/xi]*w^3 + [-lambda*xP/xi]*w^5
 *
 * so only the w^0, w^3 and w^5 coefficients are non-zero. This is the same
 * (0,3,5) sparsity the legacy code produced, which is a useful cross-check that
 * the derivation matches the tower's conventions.
 *
 * Rather than divide by xi, the whole line is scaled by xi and by the common
 * denominator. Both factors live in Fp2, and any Fp2 factor is killed by the
 * final exponentiation because (p^2-1) divides (p^12-1)/r, so the pairing value
 * is unchanged. That removes every division from the loop.
 *
 * For a Jacobian T = (X, Y, Z), doubling gives lambda = 3X^2/(2YZ) and
 * x_T = X/Z^2, y_T = Y/Z^3. Scaling by xi*2YZ^3 leaves
 *
 *     c0 = xi * yP * (2YZ^3)
 *     c3 = 3X^3 - 2Y^2
 *     c5 = -3X^2 * Z^2 * xP
 *
 * and for mixed addition with affine Q, writing H = xQ*Z^2 - X and
 * R = yQ*Z^3 - Y, scaling by xi*H*Z^3 leaves
 *
 *     c0 = xi * yP * (H*Z^3)
 *     c3 = R*X - Y*H
 *     c5 = -R * Z^2 * xP
 */
#include "elips/pairing.h"
#include <string.h>

/* f *= (c0 + c3 w^3 + c5 w^5).
 *
 * In storage terms the sparse element is L0 = (c0,0,0), L1 = (0,c3,c5), so
 * Karatsuba over fp6 costs 3 + 6 + 6 = 15 fp2 multiplications against 18 for a
 * dense fp12 multiply.
 *
 * ponytail: the legacy code reached 2 non-trivial coefficients by rescaling P
 * so that yP became 1. Worth copying if the final measurement asks for it; the
 * ceiling here is those 15 multiplications. */
/* Takes the two fp6 halves separately rather than the fp12.
 *
 * Not a style choice: given an fp12_t parameter, GCC narrows what it believes
 * f[1] to be the moment the body indexes into it, and then reports every later
 * whole-fp6 access to f[1] as a 288-byte write into a 96-byte object. The
 * object really is 288 bytes and the code was correct, but a false warning that
 * cannot be distinguished from a true one is worth the signature change.
 * Passing the halves gives GCC types it reasons about correctly, and costs
 * nothing at runtime. */
static void fp12_mul_sparse035(fp6_t f0, fp6_t f1,
                               const fp2_t c0, const fp2_t c3, const fp2_t c5)
{
    fp6_t t0, t1, s, u;

    /* t0 = a0 * (c0,0,0) */
    fp2_mul(t0[0], f0[0], c0);
    fp2_mul(t0[1], f0[1], c0);
    fp2_mul(t0[2], f0[2], c0);

    /* t1 = a1 * (0,c3,c5), using v^3 = xi:
     *   r0 = (a1*c5 + a2*c3)*xi
     *   r1 =  a0*c3 + a2*c5*xi
     *   r2 =  a0*c5 + a1*c3   */
    {
        const fp_t *a0 = f1[0], *a1 = f1[1], *a2 = f1[2];
        fp2_t m1, m2, m3, m4, m5, m6;
        fp2_mul(m1, a1, c5); fp2_mul(m2, a2, c3);
        fp2_add(t1[0], m1, m2); fp2_mul_xi(t1[0], t1[0]);
        fp2_mul(m3, a0, c3); fp2_mul(m4, a2, c5); fp2_mul_xi(m4, m4);
        fp2_add(t1[1], m3, m4);
        fp2_mul(m5, a0, c5); fp2_mul(m6, a1, c3);
        fp2_add(t1[2], m5, m6);
    }

    /* s = (a0 + a1) * (c0, c3, c5) */
    fp6_add(s, f0, f1);
    {
        fp6_t L; fp2_copy(L[0], c0); fp2_copy(L[1], c3); fp2_copy(L[2], c5);
        fp6_mul(s, s, L);
    }
    fp6_sub(s, s, t0);
    fp6_sub(s, s, t1);

    fp6_mul_v(u, t1);
    fp6_add(f0, t0, u);
    fp6_copy(f1, s);
}

/* T <- 2T, and f *= the tangent line at T evaluated at P.
 *
 * In homogeneous coordinates x_T = X/Z and y_T = Y/Z, so the tangent slope is
 * 3X^2/(2YZ). Scaling the line by xi*2YZ^2 clears every denominator and leaves
 *
 *     c0 = xi * yP * (2 Y Z^2)
 *     c3 = 3X^3 - 2 Y^2 Z
 *     c5 = -3 X^2 Z * xP
 *
 * which is a term simpler than the Jacobian version this replaced. The xi and
 * the common denominator both live in Fp2 and die in the final exponentiation,
 * so the pairing value is unchanged. */
static void dbl_step(fp12_t f, ep2_t *T, const fp_t px, const fp_t py)
{
    /* The line and the doubled point are computed together so that X^2, Y^2,
     * Z^2, XY and YZ are each formed once. Splitting them into a line
     * evaluation followed by a call to ep2_dbl cost 11% of the Miller loop,
     * measured -- the shared subexpressions are most of the work. */
    fp2_t XX, YY, ZZ, XY, YZ, t, u, c0, c3, c5, b3;
    fp2_t z3, x3, y3;

    fp2_sqr(XX, T->x);
    fp2_sqr(YY, T->y);
    fp2_sqr(ZZ, T->z);
    fp2_mul(XY, T->x, T->y);
    fp2_mul(YZ, T->y, T->z);

    /* ---- line: c3 = 3X^3 - 2Y^2 Z, c5 = -3X^2 Z xP, c0 = xi yP 2YZ^2 ---- */
    fp2_mul(t, XX, T->x);                       /* X^3        */
    fp2_add(u, t, t); fp2_add(c3, u, t);        /* 3X^3       */
    fp2_mul(t, YY, T->z); fp2_add(t, t, t);     /* 2 Y^2 Z    */
    fp2_sub(c3, c3, t);

    fp2_mul(t, XX, T->z);                       /* X^2 Z      */
    fp2_add(u, t, t); fp2_add(u, u, t);         /* 3 X^2 Z    */
    fp2_mul_fp(c5, u, px);
    fp2_neg(c5, c5);

    fp2_mul(t, YZ, T->z);                       /* Y Z^2      */
    fp2_add(t, t, t);
    fp2_mul_fp(c0, t, py);
    fp2_mul_xi(c0, c0);

    /* ---- RCB doubling, reusing the same squares ---- */
    ep2_curve_b(b3);
    fp2_add(t, b3, b3); fp2_add(b3, t, b3);     /* 3b */

    fp2_add(z3, YY, YY); fp2_add(z3, z3, z3); fp2_add(z3, z3, z3);  /* 8Y^2 */
    fp2_mul(t, b3, ZZ);                         /* b3 Z^2 */
    fp2_mul(x3, t, z3);
    fp2_add(y3, YY, t);
    fp2_mul(z3, YZ, z3);
    fp2_add(u, t, t); fp2_add(u, u, t);         /* 3 b3 Z^2 */
    fp2_sub(u, YY, u);                          /* Y^2 - 3b3Z^2 */
    fp2_mul(y3, u, y3);
    fp2_add(y3, x3, y3);
    fp2_mul(x3, u, XY);
    fp2_add(x3, x3, x3);

    fp2_copy(T->x, x3); fp2_copy(T->y, y3); fp2_copy(T->z, z3);

    fp12_mul_sparse035(f[0], f[1], c0, c3, c5);
}

/* T <- T + Q (Q affine), and f *= the chord line evaluated at P.
 *
 * With H = xQ*Z - X and R = yQ*Z - Y the slope is R/H, and scaling by xi*H*Z
 * leaves
 *
 *     c0 = xi * yP * H * Z
 *     c3 = R*X - Y*H
 *     c5 = -R * Z * xP
 */
static void add_step(fp12_t f, ep2_t *T, const fp2_t qx, const fp2_t qy,
                     const fp_t px, const fp_t py)
{
    fp2_t H, R, t, c0, c3, c5;
    ep2_t Qp;

    fp2_mul(H, qx, T->z); fp2_sub(H, H, T->x);
    fp2_mul(R, qy, T->z); fp2_sub(R, R, T->y);

    fp2_mul(c3, R, T->x);
    fp2_mul(t, T->y, H);
    fp2_sub(c3, c3, t);

    fp2_mul(c5, R, T->z);
    fp2_mul_fp(c5, c5, px);
    fp2_neg(c5, c5);

    fp2_mul(c0, H, T->z);
    fp2_mul_fp(c0, c0, py);
    fp2_mul_xi(c0, c0);

    ep2_from_affine(&Qp, qx, qy);
    ep2_add(T, T, &Qp);
    fp12_mul_sparse035(f[0], f[1], c0, c3, c5);
}

void pairing_miller(fp12_t f, const fp2_t qx, const fp2_t qy,
                    const fp_t px, const fp_t py)
{
    ep2_t T;
    fp2_t nqy;
    fp2_neg(nqy, qy);

    /* Seed the accumulator with the most significant digit, exactly as the
     * reference does: T = Q when it is +1, T = -Q when it is -1. */
    if (ELIPS_LOOP[ELIPS_LOOP_TOP] > 0) ep2_from_affine(&T, qx, qy);
    else                                ep2_from_affine(&T, qx, nqy);
    fp12_set_one(f);

    for (int i = ELIPS_LOOP_TOP - 1; i >= 0; i--) {
        fp12_sqr(f, f);
        dbl_step(f, &T, px, py);
        if (ELIPS_LOOP[i] > 0)      add_step(f, &T, qx, qy,  px, py);
        else if (ELIPS_LOOP[i] < 0) add_step(f, &T, qx, nqy, px, py);
    }
#ifdef ELIPS_FAMILY_BN
    /* BN's optimal ate needs two correction lines after the loop:
     *
     *     f *= l_{T, psi(Q)}(P),        T += psi(Q)
     *     f *= l_{T, -psi^2(Q)}(P),     T -= psi^2(Q)
     *
     * psi is the Frobenius carried through the untwisting map, so on affine
     * twist coordinates it is
     *     psi(x, y) = (conj(x) * gamma^-2, conj(y) * gamma^-3)
     * and psi^2 multiplies by the Fp2 norms of those, which are in Fp -- the
     * generator asserts their imaginary parts are zero. */
    {
        fp2_t q1x, q1y, q2x, q2y;
        fp2_conj(q1x, qx); fp2_mul(q1x, q1x, PSI_X);
        fp2_conj(q1y, qy); fp2_mul(q1y, q1y, PSI_Y);
        fp2_mul(q2x, qx, PSI2_X);
        fp2_mul(q2y, qy, PSI2_Y);
        fp2_neg(q2y, q2y);                 /* -psi^2(Q) */
        add_step(f, &T, q1x, q1y, px, py);
        add_step(f, &T, q2x, q2y, px, py);
    }
#endif
}

/* ------------------------------------------------------- multi-pairing ----
 *
 * A product of pairings needs only one Miller loop and one final
 * exponentiation. Every pair shares the running accumulator f, so the squaring
 * at the top of each iteration happens once rather than n times, and the whole
 * product is exponentiated once at the end instead of n times.
 *
 * The final exponentiation is the larger of the two savings and the one that
 * does not degrade: it is a fixed exponentiation, and x -> x^E is a
 * homomorphism, so
 *
 *     final_exp(f_1 * ... * f_n) == final_exp(f_1) * ... * final_exp(f_n)
 *
 * exactly. That identity is what test/pairing_test.c checks the result against,
 * and it holds bit for bit -- including the stray cube on BLS12, since both
 * sides carry it.
 *
 * The per-pair state is one ep2_t, so it is held in a fixed-size buffer and
 * long inputs are processed in chunks whose Miller values are multiplied
 * together. Chunking costs only the sharing of squarings across a chunk
 * boundary; the single final exponentiation is preserved for any n. At
 * ELIPS_MULTI_CHUNK = 8 a chunk already captures 7/8 of the available squaring
 * saving, and the common cases -- a two-term signature verification, a small
 * aggregate -- never reach a boundary at all. */

#define ELIPS_MULTI_CHUNK 8

/* One shared Miller loop over m <= ELIPS_MULTI_CHUNK pairs. */
static void miller_chunk(fp12_t f, const fp2_t *qx, const fp2_t *qy,
                         const fp_t *px, const fp_t *py, size_t m)
{
    ep2_t T[ELIPS_MULTI_CHUNK];
    fp2_t nqy[ELIPS_MULTI_CHUNK];

    for (size_t j = 0; j < m; j++) {
        fp2_neg(nqy[j], qy[j]);
        if (ELIPS_LOOP[ELIPS_LOOP_TOP] > 0) ep2_from_affine(&T[j], qx[j], qy[j]);
        else                                ep2_from_affine(&T[j], qx[j], nqy[j]);
    }
    fp12_set_one(f);

    for (int i = ELIPS_LOOP_TOP - 1; i >= 0; i--) {
        fp12_sqr(f, f);                     /* once, not once per pair */
        for (size_t j = 0; j < m; j++) {
            dbl_step(f, &T[j], px[j], py[j]);
            if (ELIPS_LOOP[i] > 0)
                add_step(f, &T[j], qx[j], qy[j],  px[j], py[j]);
            else if (ELIPS_LOOP[i] < 0)
                add_step(f, &T[j], qx[j], nqy[j], px[j], py[j]);
        }
    }
#ifdef ELIPS_FAMILY_BN
    /* BN's two correction lines are per pair, exactly as in pairing_miller. */
    for (size_t j = 0; j < m; j++) {
        fp2_t q1x, q1y, q2x, q2y;
        fp2_conj(q1x, qx[j]); fp2_mul(q1x, q1x, PSI_X);
        fp2_conj(q1y, qy[j]); fp2_mul(q1y, q1y, PSI_Y);
        fp2_mul(q2x, qx[j], PSI2_X);
        fp2_mul(q2y, qy[j], PSI2_Y);
        fp2_neg(q2y, q2y);
        add_step(f, &T[j], q1x, q1y, px[j], py[j]);
        add_step(f, &T[j], q2x, q2y, px[j], py[j]);
    }
#endif
}

void pairing_miller_multi(fp12_t f, const fp2_t *qx, const fp2_t *qy,
                          const fp_t *px, const fp_t *py, size_t n)
{
    fp12_set_one(f);                        /* n == 0 is the empty product */
    for (size_t base = 0; base < n; base += ELIPS_MULTI_CHUNK) {
        size_t m = n - base;
        if (m > ELIPS_MULTI_CHUNK) m = ELIPS_MULTI_CHUNK;
        fp12_t g;
        miller_chunk(g, qx + base, qy + base, px + base, py + base, m);
        fp12_mul(f, f, g);
    }
}

void pairing_final_exp_plain(fp12_t r, const fp12_t f)
{
    /* easy part: f^(p^6-1) then ^(p^2+1) */
    fp12_t t0, t1;
    fp12_conj(t0, f);
    fp12_inv(t1, f);
    fp12_mul(t0, t0, t1);
    fp12_frobenius(t1, t0, 2);
    fp12_mul(t0, t1, t0);

    /* hard part: raise to (p^4 - p^2 + 1)/r, computed directly.
     * Slow on purpose. This is the definition, and it is what the fast chain
     * gets validated against. */
    /* t0 is cyclotomic after the easy part, so the cheap squaring applies. */
    fp12_t acc;
    fp12_set_one(acc);
    for (int i = ELIPS_HARD_BITS - 1; i >= 0; i--) {
        fp12_sqr_cyc(acc, acc);
        if ((ELIPS_HARD_EXP[i / 64] >> (i % 64)) & 1) fp12_mul(acc, acc, t0);
    }
    fp12_copy(r, acc);
}

/* f^x over the signed digits of the MOTHER parameter x.
 *
 * Deliberately ELIPS_PARAM and not ELIPS_LOOP. They coincide on BLS12, but the
 * BN Miller loop runs over 6x+2 while its final exponentiation needs x, and
 * using the loop constant here silently computed the wrong exponent until the
 * test caught it.
 *
 * Only valid for f in the cyclotomic subgroup, where the conjugate is the
 * inverse -- which is what makes the negative digits free. Every element the
 * final exponentiation touches is cyclotomic, because the easy part put it
 * there. The digit pattern is a public curve constant, so branching on it
 * leaks nothing. */
void fp12_exp_param(fp12_t r, const fp12_t f)
{
    fp12_t acc, fi;
    fp12_conj(fi, f);
    if (ELIPS_PARAM[ELIPS_PARAM_TOP] > 0) fp12_copy(acc, f);
    else                                  fp12_copy(acc, fi);
    for (int i = ELIPS_PARAM_TOP - 1; i >= 0; i--) {
        fp12_sqr_cyc(acc, acc);       /* acc stays cyclotomic throughout */
        if (ELIPS_PARAM[i] > 0)      fp12_mul(acc, acc, f);
        else if (ELIPS_PARAM[i] < 0) fp12_mul(acc, acc, fi);
    }
    fp12_copy(r, acc);
}

#ifdef ELIPS_FAMILY_BLS12
void pairing_final_exp_fast(fp12_t r, const fp12_t f)
{
    /* easy part: f^(p^6-1)(p^2+1), which lands in the cyclotomic subgroup */
    fp12_t m, t0, t1;
    fp12_conj(t0, f);
    fp12_inv(t1, f);
    fp12_mul(t0, t0, t1);          /* f^(p^6-1) */
    fp12_frobenius(t1, t0, 2);
    fp12_mul(m, t1, t0);           /* ^(p^2+1) */

    /* hard part, computing 3*lambda via
     *     3*lambda = (x-1)^2 (x+p) (x^2+p^2-1) + 3
     * From here on every inverse is a conjugation, since m is cyclotomic. */
    fp12_t a, b, c, d, e, mi;
    fp12_conj(mi, m);

    fp12_exp_param(a, m);          /* m^x            */
    fp12_mul(a, a, mi);            /* m^(x-1)        */

    fp12_exp_param(b, a);          /* a^x            */
    fp12_conj(t0, a);
    fp12_mul(b, b, t0);            /* a^(x-1) = m^((x-1)^2) */

    fp12_exp_param(c, b);          /* b^x            */
    fp12_frobenius(t0, b, 1);
    fp12_mul(c, c, t0);            /* b^(x+p)        */

    fp12_exp_param(d, c);          /* c^x            */
    fp12_exp_param(d, d);          /* c^(x^2)        */
    fp12_frobenius(t0, c, 2);
    fp12_mul(d, d, t0);            /* * c^(p^2)      */
    fp12_conj(t0, c);
    fp12_mul(d, d, t0);            /* * c^-1  => c^(x^2+p^2-1) */

    fp12_sqr_cyc(e, m);
    fp12_mul(e, e, m);             /* m^3            */
    fp12_mul(r, d, e);
}
#else
/* Fast final exponentiation for BN.
 *
 * Unlike BLS12, the BN hard part decomposes exactly, with no stray factor.
 * Derived symbolically and reproduced here:
 *
 *     lambda = d0 + d1*p + d2*p^2 + d3*p^3
 *     d0 = -36x^3 - 30x^2 - 18x - 2
 *     d1 = -36x^3 - 18x^2 - 12x + 1
 *     d2 =            6x^2      + 1
 *     d3 =                        1
 *
 * Regrouping by which power of x each term needs, with u1 = f^x, u2 = f^(x^2)
 * and u3 = f^(x^3):
 *
 *     f^lambda = u3^(-36 - 36p)
 *              * u2^(-30 - 18p + 6p^2)
 *              * u1^(-18 - 12p)
 *              * f ^( -2 +   p +  p^2 + p^3)
 *
 * so three parameter exponentiations plus a handful of small powers. Every
 * negative exponent is a conjugation, because the easy part already put f in
 * the cyclotomic subgroup.
 *
 * This yields e exactly. The legacy BN chain raised to roughly 12*X^3 times
 * this (issue #16), so the new path is not just faster but right. */

/* a^6, then the multiples the chain needs, sharing intermediates. */
static void small_powers(fp12_t a6, fp12_t a12, fp12_t a18,
                         fp12_t a30, fp12_t a36, const fp12_t a)
{
    fp12_t t2, t4;
    fp12_sqr_cyc(t2, a);        /* a^2  ; a is cyclotomic, so all of these are */
    fp12_sqr_cyc(t4, t2);       /* a^4  */
    fp12_mul(a6, t4, t2);       /* a^6  */
    fp12_sqr_cyc(a12, a6);      /* a^12 */
    fp12_mul(a18, a12, a6);     /* a^18 */
    fp12_mul(a30, a18, a12);    /* a^30 */
    fp12_sqr_cyc(a36, a18);     /* a^36 */
}

void pairing_final_exp_fast(fp12_t r, const fp12_t f)
{
    /* easy part */
    fp12_t m, t0, t1;
    fp12_conj(t0, f);
    fp12_inv(t1, f);
    fp12_mul(t0, t0, t1);
    fp12_frobenius(t1, t0, 2);
    fp12_mul(m, t1, t0);

    fp12_t u1, u2, u3;
    fp12_exp_param(u1, m);      /* m^x     */
    fp12_exp_param(u2, u1);     /* m^(x^2) */
    fp12_exp_param(u3, u2);     /* m^(x^3) */

    fp12_t a6, a12, a18, a30, a36, acc, tmp;

    /* u3^(-36 - 36p) */
    small_powers(a6, a12, a18, a30, a36, u3);
    fp12_conj(acc, a36);                       /* u3^-36        */
    fp12_frobenius(tmp, a36, 1);
    fp12_conj(tmp, tmp);
    fp12_mul(acc, acc, tmp);                   /* * u3^(-36p)   */

    /* u2^(-30 - 18p + 6p^2) */
    small_powers(a6, a12, a18, a30, a36, u2);
    fp12_conj(tmp, a30);         fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, a18, 1); fp12_conj(tmp, tmp); fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, a6, 2);  fp12_mul(acc, acc, tmp);

    /* u1^(-18 - 12p) */
    small_powers(a6, a12, a18, a30, a36, u1);
    fp12_conj(tmp, a18);         fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, a12, 1); fp12_conj(tmp, tmp); fp12_mul(acc, acc, tmp);

    /* m^(-2 + p + p^2 + p^3) */
    fp12_sqr_cyc(tmp, m); fp12_conj(tmp, tmp); /* m^-2 */
    fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, m, 1); fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, m, 2); fp12_mul(acc, acc, tmp);
    fp12_frobenius(tmp, m, 3); fp12_mul(acc, acc, tmp);

    fp12_copy(r, acc);
}

/* Note for anyone comparing families: the identity the BLS12 chain rests on,
 *     3*lambda = (x-1)^2 (x+p) (x^2+p^2-1) + 3,
 * is specific to BLS12. Applying it to BN gives a wrong exponent, which the
 * test suite caught the moment BN was wired in. BN's own decomposition above
 * is exact, so BN returns e while BLS12 returns e^3. */
#endif

/* ---------------------------------------------------------- public API ---- */

void ep_generator(ep_t *g)
{
    fp_copy(g->x, EP_GEN_X);
    fp_copy(g->y, EP_GEN_Y);
    fp_set_one(g->z);
}

void ep2_generator(ep2_t *g)
{
    fp2_copy(g->x, EP2_GEN_X);
    fp2_copy(g->y, EP2_GEN_Y);
    fp2_set_one(g->z);
}

/* ------------------------------------------------- subgroup membership ----
 *
 * These used to be [r]P, a full-length ladder, and they dominated the cost of
 * elips_pairing. They are now endomorphism tests, following Scott (ePrint
 * 2021/1130) as corrected by El Housni, Guillevic and Piellard (ePrint
 * 2022/352) for G2, and Dai, Lin, Zhao and Zhou (ePrint 2022/348) for the
 * generalisation BN needs.
 *
 * A test that is too permissive does not fail loudly -- it accepts
 * attacker-chosen points of small order, which is the hole the test exists to
 * close. So neither the relation nor its exactness is recalled here.
 * tools/reference/subgroup_ref.py derives both per curve, and gen_params.py
 * re-asserts the exactness conditions every time it regenerates fp_params.h.
 *
 * The argument, in short. Each test asks whether E(P) == [m]P for an
 * endomorphism E satisfying a known quadratic. A point of prime order l
 * satisfying it forces m to be a root of that quadratic mod l:
 *
 *     G1 (BLS12)   phi^2 + phi + 1 = 0    needs  m^2 + m + 1 == 0 mod l
 *     G2 (both)    psi^2 - t psi + p = 0  needs  m^2 - t m + p == 0 mod l
 *
 * so the test is exact when no prime dividing the cofactor divides the
 * corresponding integer. For G1 with m = -x^2 that integer is x^4 - x^2 + 1,
 * which is r itself, so BLS12 G1 is exact for every seed. G2 has no such
 * identity and is checked per curve by a gcd. */

int ep_in_subgroup(const ep_t *p)
{
#ifdef ELIPS_FAMILY_BN
    /* #E(Fp) = r exactly on a BN curve, so the cofactor is 1 and being on the
     * curve is already the whole condition. The generator asserts h1 == 1. */
    return ep_is_infinity(p) || ep_on_curve(p);
#else
    if (ep_is_infinity(p)) return 1;
    if (!ep_on_curve(p))   return 0;

    /* phi(P) == [-x^2]P. Squaring the parameter makes the sign irrelevant, so
     * this is two ladders over |x| rather than one over r: 128 bits against
     * 255 on BLS12-381. */
    ep_t lhs, rhs;
    ep_phi(&lhs, p);
    ep_mul(&rhs, p,    ELIPS_ABSX, ELIPS_ABSX_BITS);
    ep_mul(&rhs, &rhs, ELIPS_ABSX, ELIPS_ABSX_BITS);
    ep_neg(&rhs, &rhs);
    return ep_eq(&lhs, &rhs);
#endif
}

int ep2_in_subgroup(const ep2_t *q)
{
    if (ep2_is_infinity(q)) return 1;
    if (!ep2_on_curve(q))   return 0;

    /* psi(Q) == [p mod r]Q. That multiplier is x on BLS12 and 6x^2 on BN. */
    ep2_t lhs, rhs;
    ep2_psi(&lhs, q);
#ifdef ELIPS_FAMILY_BLS12
    ep2_mul(&rhs, q, ELIPS_ABSX, ELIPS_ABSX_BITS);
#if ELIPS_X_NEGATIVE
    ep2_neg(&rhs, &rhs);
#endif
#else
    ep2_mul(&rhs, q, ELIPS_6XSQ, ELIPS_6XSQ_BITS);
#endif
    return ep2_eq(&lhs, &rhs);
}

int elips_pairing(fp12_t out, const ep_t *P, const ep2_t *Q)
{
    fp12_set_one(out);
    if (ep_is_infinity(P) || ep2_is_infinity(Q)) return 0;
    if (!ep_in_subgroup(P) || !ep2_in_subgroup(Q)) return 0;

    fp_t px, py;
    fp2_t qx, qy;
    if (!ep_to_affine(px, py, P))  return 0;
    if (!ep2_to_affine(qx, qy, Q)) return 0;

    fp12_t f;
    pairing_miller(f, qx, qy, px, py);
    pairing_final_exp_fast(out, f);
    return 1;
}

int elips_pairing_multi(fp12_t out, const ep_t *P, const ep2_t *Q, size_t n)
{
    fp12_set_one(out);

    /* Validate everything before computing anything. A caller that ignores the
     * return value then gets one, not a product over the pairs that happened to
     * be checked before the bad one. The contract matches elips_pairing: the
     * identity is rejected rather than treated as a factor of one. */
    for (size_t j = 0; j < n; j++) {
        if (ep_is_infinity(&P[j]) || ep2_is_infinity(&Q[j])) return 0;
        if (!ep_in_subgroup(&P[j]) || !ep2_in_subgroup(&Q[j])) return 0;
    }

    fp_t  px[ELIPS_MULTI_CHUNK], py[ELIPS_MULTI_CHUNK];
    fp2_t qx[ELIPS_MULTI_CHUNK], qy[ELIPS_MULTI_CHUNK];
    fp12_t f, g;
    fp12_set_one(f);

    size_t m = 0;
    for (size_t j = 0; j < n; j++) {
        ep_to_affine(px[m], py[m], &P[j]);
        ep2_to_affine(qx[m], qy[m], &Q[j]);
        if (++m == ELIPS_MULTI_CHUNK) {
            miller_chunk(g, qx, qy, px, py, m);
            fp12_mul(f, f, g);
            m = 0;
        }
    }
    if (m) {
        miller_chunk(g, qx, qy, px, py, m);
        fp12_mul(f, f, g);
    }

    pairing_final_exp_fast(out, f);         /* once, whatever n is */
    return 1;
}
