/*
 * Group law in HOMOGENEOUS PROJECTIVE coordinates: x = X/Z, y = Y/Z, and the
 * identity is Z = 0. Included from ec_tmpl.h when EC_JACOBIAN is not set.
 *
 * This is what G1 uses, and since the twist went Jacobian it is all that uses
 * it. Over Fp a squaring IS a multiplication -- fp_sqr is fp_mul(a, a) -- so
 * the Jacobian trade of two doubling multiplies for four addition ones comes
 * out roughly level, and the completeness of the RCB formulas is worth more
 * than level. See issue #50.
 */

/* RCB Algorithm 7: complete addition for a = 0.
 * Correct for every input pair, including equal, opposite and identity. */
void PT(add)(PTT *r, const PTT *p, const PTT *q)
{
    EC_FT t0, t1, t2, t3, t4, x3, y3, z3, b3;
    PT(curve_b3)(b3);

    F(mul)(t0, p->x, q->x);
    F(mul)(t1, p->y, q->y);
    F(mul)(t2, p->z, q->z);

    F(add)(t3, p->x, p->y);
    F(add)(t4, q->x, q->y);
    F(mul)(t3, t3, t4);
    F(add)(t4, t0, t1);
    F(sub)(t3, t3, t4);

    F(add)(t4, p->y, p->z);
    F(add)(x3, q->y, q->z);
    F(mul)(t4, t4, x3);
    F(add)(x3, t1, t2);
    F(sub)(t4, t4, x3);

    F(add)(x3, p->x, p->z);
    F(add)(y3, q->x, q->z);
    F(mul)(x3, x3, y3);
    F(add)(y3, t0, t2);
    F(sub)(y3, x3, y3);

    F(add)(x3, t0, t0);
    F(add)(t0, x3, t0);          /* t0 = 3 X1X2 */
    F(mul)(t2, b3, t2);          /* t2 = b3 Z1Z2 */

    F(add)(z3, t1, t2);
    F(sub)(t1, t1, t2);
    F(mul)(y3, b3, y3);

    F(mul)(x3, t4, y3);
    F(mul)(t2, t3, t1);
    F(sub)(x3, t2, x3);

    F(mul)(y3, y3, t0);
    F(mul)(t1, t1, z3);
    F(add)(y3, t1, y3);

    F(mul)(t0, t0, t3);
    F(mul)(z3, z3, t4);
    F(add)(z3, z3, t0);

    F(copy)(r->x, x3); F(copy)(r->y, y3); F(copy)(r->z, z3);
}

/* RCB Algorithm 9: exception-free doubling for a = 0.
 *
 * Y^2 and Z^2 are written F(mul)(t, a, a) rather than F(sqr): over Fp, which is
 * the only field that reaches this file now, fp_sqr IS fp_mul(a, a), so routing
 * through it buys nothing and adds a call, measured at +1.3% on ep_mul at 100%
 * agreement. The dedicated 4M + 5S homogeneous doubling that used to sit beside
 * this one was for the twist and is gone with it; see issue #50. */
void PT(dbl)(PTT *r, const PTT *p)
{
    EC_FT t0, t1, t2, x3, y3, z3, b3;
    PT(curve_b3)(b3);

    F(mul)(t0, p->y, p->y);
    F(add)(z3, t0, t0);
    F(add)(z3, z3, z3);
    F(add)(z3, z3, z3);          /* z3 = 8 Y^2 */

    F(mul)(t1, p->y, p->z);
    F(mul)(t2, p->z, p->z);
    F(mul)(t2, b3, t2);          /* t2 = b3 Z^2 */

    F(mul)(x3, t2, z3);
    F(add)(y3, t0, t2);
    F(mul)(z3, t1, z3);

    F(add)(t1, t2, t2);
    F(add)(t2, t1, t2);          /* t2 = 3 b3 Z^2 */
    F(sub)(t0, t0, t2);          /* t0 = Y^2 - 3 b3 Z^2 */

    F(mul)(y3, t0, y3);
    F(add)(y3, x3, y3);

    F(mul)(t1, p->x, p->y);
    F(mul)(x3, t0, t1);
    F(add)(x3, x3, x3);

    F(copy)(r->x, x3); F(copy)(r->y, y3); F(copy)(r->z, z3);
}

int PT(to_affine)(EC_FT x, EC_FT y, const PTT *p)
{
    if (PT(is_infinity)(p)) { F(set_zero)(x); F(set_zero)(y); return 0; }
    EC_FT zi;
    F(inv)(zi, p->z);
    F(mul)(x, p->x, zi);
    F(mul)(y, p->y, zi);
    return 1;
}

int PT(on_curve)(const PTT *p)
{
    if (PT(is_infinity)(p)) return 1;
    EC_FT lhs, rhs, z3, b;
    PT(curve_b)(b);
    F(mul)(lhs, p->y, p->y); F(mul)(lhs, lhs, p->z);
    F(mul)(rhs, p->x, p->x); F(mul)(rhs, rhs, p->x);
    F(mul)(z3, p->z, p->z);  F(mul)(z3, z3, p->z);
    F(mul)(z3, z3, b);
    F(add)(rhs, rhs, z3);
    return F(eq)(lhs, rhs);
}

/* (X1:Y1:Z1) == (X2:Y2:Z2) iff X1*Z2 == X2*Z1 and Y1*Z2 == Y2*Z1.
 *
 * Two points at infinity have Z = 0 and satisfy both, as they should. An
 * infinity never compares equal to an affine point: with (0:Y1:0) the X test
 * is 0 == 0 and passes, but the Y test needs Y1*Z2 == 0 with Y1 and Z2 both
 * nonzero, so it fails. No inversion and no branch on coordinate values. */
int PT(eq)(const PTT *a, const PTT *b)
{
    EC_FT t0, t1;
    F(mul)(t0, a->x, b->z); F(mul)(t1, b->x, a->z);
    if (!F(eq)(t0, t1)) return 0;
    F(mul)(t0, a->y, b->z); F(mul)(t1, b->y, a->z);
    return F(eq)(t0, t1);
}

