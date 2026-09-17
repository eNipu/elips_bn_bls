/*
 * Group law in JACOBIAN coordinates: x = X/Z^2, y = Y/Z^3, and the identity is
 * Z = 0. Included from ec_tmpl.h, and instantiated for both E(Fp) and E'(Fp2).
 *
 * WHY, and what it costs. The doubling is 2M + 5S here against 4M + 5S for the
 * dedicated homogeneous form, and the doubling is 89% of ep2_in_subgroup and
 * most of every scalar multiplication. The addition goes the other way: a
 * complete Jacobian addition has to be the unified add-or-double below at
 * 13M + 5S, against 14M for the homogeneous RCB formulas this replaced. The
 * trade is worth taking wherever doublings outnumber additions, and it is a
 * trade and not a free win. Issues #50 (the twist) and #51 (G1).
 *
 * The trade is sharper over Fp than over Fp2, because fp_sqr is literally
 * fp_mul(a, a) -- S = M exactly -- while Fp2 has a real squaring at S/M = 0.62.
 * So on G1 the doubling saves 2M and the addition costs 4M, and which way a
 * caller comes out depends only on its ratio of the two. Measured on
 * BLS12-381, against the homogeneous formulas with their curve constant
 * hoisted out of the inner loop:
 *
 *      ep_in_subgroup   128 dbl :   6 add     -11.5%
 *      ep_mul             4 dbl :   1 add      -3.9%
 *      ep_mul_glv         2 dbl :   1 add      +4.4%   (slower)
 *
 * The GLV ladder is the one caller that loses, and it is the signing path.
 * That is the cost of having one group law rather than two; it was taken
 * knowingly and issue #51 records the numbers. A wider ladder window would
 * shift the ratio back, but the arithmetic in ec.c refutes it: at width three
 * the 64-entry table costs more additions to build than the narrower digits
 * save.
 *
 * NOT COMPLETE, unlike the RCB formulas that came before. Both routines below
 * carry an explicit exception argument instead, and test/ec_group_test.c
 * checks both of them, for both groups, against an affine oracle on random
 * off-subgroup points. An argument is not a check.
 */

/* Unified addition: 13M + 5S, correct for every input pair.
 *
 * Jacobian addition is not complete the way the RCB formulas were, and this
 * routine is reached from ep_in_subgroup and ep2_in_subgroup with points an
 * attacker chose, so
 * every degenerate pair has to be handled rather than argued away. It is
 * handled the way blst's POINTonE1_dadd and POINTonE2_dadd handle it: compute
 * the addition and
 * the doubling helpers side by side, then select between them with a mask.
 *
 *   P == Q  gives H = 0 and R = 0. The doubling helpers are selected, and
 *           substituting them reproduces dbl-2009-l exactly: Z3 = 2 Y1 Z1,
 *           X3 = 9X1^4 - 8 X1 Y1^2, Y3 = 3X1^2 (4 X1 Y1^2 - X3) - 8 Y1^4.
 *   P == -Q gives H = 0 but R != 0, so the addition path runs and
 *           Z3 = H Z1 Z2 = 0, which is the identity. Correct.
 *   Either input at infinity is handled by the two masked selects at the end,
 *           which return the other input.
 *
 * Every choice is a mask over both candidates. Nothing branches on a
 * coordinate, which tools/verify/ct_branch_scan.py checks. */
void PT(add)(PTT *r, const PTT *p, const PTT *q)
{
    EC_FT z1z1, z2z2, u1, u2, s1, s2, h, rr, sx, zz;
    EC_FT dh, dr, dsx, t0, t1, x3, y3, z3;
    limb_t is_dbl, p1inf, p2inf;

    /* the doubling helpers, from p alone */
    F(add)(dh, p->y, p->y);                     /* 2 Y1   */
    F(sqr)(dr, p->x);
    F(add)(t0, dr, dr); F(add)(dr, t0, dr);     /* 3 X1^2 */
    F(add)(dsx, p->x, p->x);                    /* 2 X1   */

    F(sqr)(z1z1, p->z);
    F(sqr)(z2z2, q->z);
    F(mul)(zz, p->z, q->z);                     /* Z1 Z2  */

    F(mul)(u1, p->x, z2z2);                     /* U1 = X1 Z2^2 */
    F(mul)(u2, q->x, z1z1);                     /* U2 = X2 Z1^2 */
    F(mul)(s1, p->y, q->z); F(mul)(s1, s1, z2z2);   /* S1 = Y1 Z2^3 */
    F(mul)(s2, q->y, p->z); F(mul)(s2, s2, z1z1);   /* S2 = Y2 Z1^3 */

    F(sub)(h,  u2, u1);
    F(sub)(rr, s2, s1);
    F(add)(sx, u1, u2);

    is_dbl = (limb_t)0 - (limb_t)(F(is_zero)(h) & F(is_zero)(rr));
    F(cselect)(h,   dh,    h,   is_dbl);
    F(cselect)(rr,  dr,    rr,  is_dbl);
    F(cselect)(sx,  dsx,   sx,  is_dbl);
    F(cselect)(u1,  p->x,  u1,  is_dbl);
    F(cselect)(s1,  p->y,  s1,  is_dbl);
    F(cselect)(zz,  p->z,  zz,  is_dbl);

    F(mul)(z3, zz, h);                          /* Z3 = H Z1 Z2 */
    F(sqr)(t0, h);                              /* HH           */
    F(mul)(t1, t0, h);                          /* HHH          */
    F(mul)(t1, t1, s1);                         /* HHH S1       */
    F(mul)(u1, t0, u1);                         /* HH U1        */
    F(mul)(t0, t0, sx);                         /* HH sx        */
    F(sqr)(x3, rr);
    F(sub)(x3, x3, t0);                         /* X3 = R^2 - HH sx */
    F(sub)(y3, u1, x3);
    F(mul)(y3, y3, rr);
    F(sub)(y3, y3, t1);                         /* Y3 = R(HH U1 - X3) - HHH S1 */

    /* Either input at infinity returns the other. Computed from the inputs,
     * and r is written only at the end, so r may alias p or q. */
    p2inf = (limb_t)0 - (limb_t)F(is_zero)(q->z);
    F(cselect)(x3, p->x, x3, p2inf);
    F(cselect)(y3, p->y, y3, p2inf);
    F(cselect)(z3, p->z, z3, p2inf);
    p1inf = (limb_t)0 - (limb_t)F(is_zero)(p->z);
    F(cselect)(x3, q->x, x3, p1inf);
    F(cselect)(y3, q->y, y3, p1inf);
    F(cselect)(z3, q->z, z3, p1inf);

    F(copy)(r->x, x3); F(copy)(r->y, y3); F(copy)(r->z, z3);
}

/* dbl-2009-l, the dedicated doubling for a = 0: 2M + 5S.
 *
 *   A = X^2,  B = Y^2,  C = B^2,  D = 2((X+B)^2 - A - C),  E = 3A,  F = E^2
 *   X3 = F - 2D,   Y3 = E(D - X3) - 8C,   Z3 = 2 Y Z
 *
 * Exception-free, which matters because this runs on untrusted points:
 *
 *   Z = 0. Z3 = 2 Y * 0 = 0, so the identity stays the identity.
 *   Y = 0, that is 2-torsion. B = C = 0, so D = 2((X+0)^2 - X^2 - 0) = 0, and
 *          Z3 = 2 * 0 * Z = 0. The identity, which is what doubling a
 *          2-torsion point should give.
 *
 * Neither case can produce a point that is not the identity, and the identity
 * is Z = 0 regardless of X and Y, so there is no degenerate output to guard
 * against the way there is in homogeneous coordinates. */
void PT(dbl)(PTT *r, const PTT *p)
{
    EC_FT a, b, c, d, e, f, t, x3, y3, z3;

    F(sqr)(a, p->x);                            /* A = X^2 */
    F(sqr)(b, p->y);                            /* B = Y^2 */
    F(sqr)(c, b);                               /* C = B^2 */

    F(add)(t, p->x, b);
    F(sqr)(t, t);                               /* (X+B)^2 */
    F(sub)(t, t, a);
    F(sub)(t, t, c);
    F(add)(d, t, t);                            /* D = 2((X+B)^2 - A - C) */

    F(add)(e, a, a); F(add)(e, e, a);           /* E = 3A  */
    F(sqr)(f, e);                               /* F = E^2 */

    F(sub)(x3, f, d);
    F(sub)(x3, x3, d);                          /* X3 = F - 2D */

    F(mul)(z3, p->y, p->z);
    F(add)(z3, z3, z3);                         /* Z3 = 2 Y Z  */

    F(add)(c, c, c); F(add)(c, c, c); F(add)(c, c, c);   /* 8C */
    F(sub)(y3, d, x3);
    F(mul)(y3, y3, e);
    F(sub)(y3, y3, c);                          /* Y3 = E(D - X3) - 8C */

    F(copy)(r->x, x3); F(copy)(r->y, y3); F(copy)(r->z, z3);
}

int PT(to_affine)(EC_FT x, EC_FT y, const PTT *p)
{
    if (PT(is_infinity)(p)) { F(set_zero)(x); F(set_zero)(y); return 0; }
    EC_FT zi, zi2, zi3;
    F(inv)(zi, p->z);
    F(sqr)(zi2, zi);                            /* 1/Z^2 */
    F(mul)(zi3, zi2, zi);                       /* 1/Z^3 */
    F(mul)(x, p->x, zi2);
    F(mul)(y, p->y, zi3);
    return 1;
}

/* Y^2 = X^3 + b Z^6. */
int PT(on_curve)(const PTT *p)
{
    if (PT(is_infinity)(p)) return 1;
    EC_FT lhs, rhs, z2, z6, b;
    PT(curve_b)(b);
    F(sqr)(lhs, p->y);
    F(sqr)(rhs, p->x);  F(mul)(rhs, rhs, p->x);
    F(sqr)(z2, p->z);
    F(sqr)(z6, z2);     F(mul)(z6, z6, z2);     /* Z^6 */
    F(mul)(z6, z6, b);
    F(add)(rhs, rhs, z6);
    return F(eq)(lhs, rhs);
}

/* X1 Z2^2 == X2 Z1^2 and Y1 Z2^3 == Y2 Z1^3, cross-multiplied rather than
 * normalised, so there is no inversion and no branch on a coordinate. Two
 * points at infinity have Z = 0 on both sides and compare equal, which is what
 * they should do. */
int PT(eq)(const PTT *a, const PTT *b)
{
    EC_FT az2, bz2, t0, t1;
    F(sqr)(az2, a->z);
    F(sqr)(bz2, b->z);
    F(mul)(t0, a->x, bz2); F(mul)(t1, b->x, az2);
    if (!F(eq)(t0, t1)) return 0;
    F(mul)(az2, az2, a->z);                     /* Z1^3 */
    F(mul)(bz2, bz2, b->z);                     /* Z2^3 */
    F(mul)(t0, a->y, bz2); F(mul)(t1, b->y, az2);
    return F(eq)(t0, t1);
}
