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
static void fp12_mul_sparse035(fp12_t f, const fp2_t c0, const fp2_t c3, const fp2_t c5)
{
    fp6_t t0, t1, s, u;

    /* t0 = a0 * (c0,0,0) */
    fp2_mul(t0[0], f[0][0], c0);
    fp2_mul(t0[1], f[0][1], c0);
    fp2_mul(t0[2], f[0][2], c0);

    /* t1 = a1 * (0,c3,c5), using v^3 = xi:
     *   r0 = (a1*c5 + a2*c3)*xi
     *   r1 =  a0*c3 + a2*c5*xi
     *   r2 =  a0*c5 + a1*c3   */
    {
        const fp_t *a0 = f[1][0], *a1 = f[1][1], *a2 = f[1][2];
        fp2_t m1, m2, m3, m4, m5, m6;
        fp2_mul(m1, a1, c5); fp2_mul(m2, a2, c3);
        fp2_add(t1[0], m1, m2); fp2_mul_xi(t1[0], t1[0]);
        fp2_mul(m3, a0, c3); fp2_mul(m4, a2, c5); fp2_mul_xi(m4, m4);
        fp2_add(t1[1], m3, m4);
        fp2_mul(m5, a0, c5); fp2_mul(m6, a1, c3);
        fp2_add(t1[2], m5, m6);
    }

    /* s = (a0 + a1) * (c0, c3, c5) */
    fp6_add(s, f[0], f[1]);
    {
        fp6_t L; fp2_copy(L[0], c0); fp2_copy(L[1], c3); fp2_copy(L[2], c5);
        fp6_mul(s, s, L);
    }
    fp6_sub(s, s, t0);
    fp6_sub(s, s, t1);

    fp6_mul_v(u, t1);
    fp6_add(f[0], t0, u);
    fp6_copy(f[1], s);
}

/* T <- 2T, and f *= the tangent line at T evaluated at P. */
static void dbl_step(fp12_t f, ep2_t *T, const fp_t px, const fp_t py)
{
    fp2_t A, B, C, D, E, FF, t, ZZ, c0, c3, c5;

    fp2_sqr(A, T->x);                       /* X^2   */
    fp2_sqr(B, T->y);                       /* Y^2   */
    fp2_sqr(C, B);                          /* Y^4   */
    fp2_sqr(ZZ, T->z);                      /* Z^2   */

    /* c3 = 3X^3 - 2Y^2 */
    fp2_mul(c3, A, T->x);                   /* X^3 */
    fp2_add(t, c3, c3); fp2_add(c3, t, c3); /* 3X^3 */
    fp2_add(t, B, B);                       /* 2Y^2 */
    fp2_sub(c3, c3, t);

    /* c5 = -3 X^2 Z^2 xP */
    fp2_add(t, A, A); fp2_add(t, t, A);     /* 3X^2 */
    fp2_mul(c5, t, ZZ);
    fp2_mul_fp(c5, c5, px);
    fp2_neg(c5, c5);

    /* c0 = xi * yP * (2 Y Z^3) ; note Z3 = 2YZ so 2YZ^3 = Z3 * Z^2 */
    fp2_mul(t, T->y, T->z);
    fp2_add(t, t, t);                       /* Z3 = 2YZ */
    fp2_t Z3; fp2_copy(Z3, t);
    fp2_mul(c0, t, ZZ);                     /* 2 Y Z^3 */
    fp2_mul_fp(c0, c0, py);
    fp2_mul_xi(c0, c0);

    /* the doubling itself (dbl-2009-l, a = 0) */
    fp2_add(D, T->x, B); fp2_sqr(D, D); fp2_sub(D, D, A); fp2_sub(D, D, C);
    fp2_add(D, D, D);
    fp2_add(E, A, A); fp2_add(E, E, A);
    fp2_sqr(FF, E);
    fp2_add(t, D, D); fp2_sub(t, FF, t);
    fp2_copy(T->z, Z3);
    fp2_sub(D, D, t); fp2_mul(D, E, D);
    fp2_add(C, C, C); fp2_add(C, C, C); fp2_add(C, C, C);
    fp2_sub(T->y, D, C);
    fp2_copy(T->x, t);

    fp12_mul_sparse035(f, c0, c3, c5);
}

/* T <- T + Q (Q affine), and f *= the chord line evaluated at P. */
static void add_step(fp12_t f, ep2_t *T, const fp2_t qx, const fp2_t qy,
                     const fp_t px, const fp_t py)
{
    fp2_t ZZ, U, S, H, R, HH, I, J, V, t, c0, c3, c5;

    fp2_sqr(ZZ, T->z);
    fp2_mul(U, qx, ZZ);                     /* xQ Z^2      */
    fp2_mul(S, qy, T->z); fp2_mul(S, S, ZZ);/* yQ Z^3      */
    fp2_sub(H, U, T->x);                    /* H = xQZ^2-X */
    fp2_sub(R, S, T->y);                    /* R = yQZ^3-Y */

    /* c3 = R*X - Y*H */
    fp2_mul(c3, R, T->x);
    fp2_mul(t, T->y, H);
    fp2_sub(c3, c3, t);

    /* c5 = -R Z^2 xP */
    fp2_mul(c5, R, ZZ);
    fp2_mul_fp(c5, c5, px);
    fp2_neg(c5, c5);

    /* c0 = xi * yP * H * Z^3 */
    fp2_mul(c0, H, ZZ); fp2_mul(c0, c0, T->z);
    fp2_mul_fp(c0, c0, py);
    fp2_mul_xi(c0, c0);

    /* mixed addition (madd-2007-bl) */
    fp2_sqr(HH, H);
    fp2_add(I, HH, HH); fp2_add(I, I, I);   /* 4 H^2 */
    fp2_mul(J, H, I);
    fp2_add(R, R, R);                       /* r = 2R */
    fp2_mul(V, T->x, I);
    fp2_sqr(t, R); fp2_sub(t, t, J);
    fp2_t tv; fp2_add(tv, V, V); fp2_sub(t, t, tv);
    fp2_t X3; fp2_copy(X3, t);
    fp2_sub(t, V, X3); fp2_mul(t, R, t);
    fp2_mul(tv, T->y, J); fp2_add(tv, tv, tv);
    fp2_sub(T->y, t, tv);
    fp2_add(t, T->z, H); fp2_sqr(t, t); fp2_sub(t, t, ZZ); fp2_sub(t, t, HH);
    fp2_copy(T->z, t);
    fp2_copy(T->x, X3);

    fp12_mul_sparse035(f, c0, c3, c5);
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
#  error "BN needs the two post-loop correction lines; not implemented yet."
#endif
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
    fp12_exp(r, t0, ELIPS_HARD_EXP, ELIPS_HARD_BITS);
}
