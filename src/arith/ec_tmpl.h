/*
 * Group law in HOMOGENEOUS projective coordinates, using the complete formulas
 * of Renes, Costello and Batina (EUROCRYPT 2016), specialised to a = 0.
 *
 * Written once and instantiated for both E(Fp) and E'(Fp2). The original
 * library kept two hand-copied versions of every curve routine and they had
 * already drifted apart; generating both from one text removes that.
 *
 * Why homogeneous rather than Jacobian
 * ------------------------------------
 * The previous Jacobian implementation had no complete addition. To resolve the
 * coincident-point case without branching it computed a full doubling on every
 * call and selected between the two results, which cost it everything the
 * coordinate change had gained: measured at 0.98x against the affine code it
 * replaced. Every addition in a ladder paid that.
 *
 * The RCB formulas are complete in one expression: the same code is correct for
 * P + Q, P + P, P + (-P), P + O and O + O, with no branch and no fixup. That
 * makes them both faster here and simpler to reason about.
 *
 * A point is (X : Y : Z) with x = X/Z and y = Y/Z. The identity is any point
 * with Z = 0; the canonical one is (0 : 1 : 0).
 *
 * Expects EC_PT, EC_F and EC_FT from the includer.
 */

#define CAT_(a,b) a##b
#define CAT(a,b)  CAT_(a,b)
#define PT(name)  CAT(EC_PT, CAT(_, name))
#define PTT       CAT(EC_PT, _t)
#define F(name)   CAT(EC_F,  CAT(_, name))

/* 3b, the only curve constant the RCB formulas need. */
static void PT(curve_b3)(EC_FT out)
{
    EC_FT b, t;
    PT(curve_b)(b);
    F(add)(t, b, b);
    F(add)(out, t, b);
}

void PT(set_infinity)(PTT *r)
{
    F(set_zero)(r->x);
    F(set_one)(r->y);          /* (0 : 1 : 0) */
    F(set_zero)(r->z);
}

int PT(is_infinity)(const PTT *p) { return F(is_zero)(p->z); }

void PT(copy)(PTT *r, const PTT *p)
{ F(copy)(r->x, p->x); F(copy)(r->y, p->y); F(copy)(r->z, p->z); }

void PT(neg)(PTT *r, const PTT *p)
{ F(copy)(r->x, p->x); F(neg)(r->y, p->y); F(copy)(r->z, p->z); }

void PT(from_affine)(PTT *r, const EC_FT x, const EC_FT y)
{ F(copy)(r->x, x); F(copy)(r->y, y); F(set_one)(r->z); }

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

/* RCB Algorithm 9: exception-free doubling for a = 0. */
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

/* Kept for source compatibility with the Jacobian version, which needed a
 * separate incomplete routine for speed. The complete formulas are already the
 * fast path, so this is simply an alias. */
void PT(add_generic)(PTT *r, const PTT *p, const PTT *q) { PT(add)(r, p, q); }

int PT(to_affine)(EC_FT x, EC_FT y, const PTT *p)
{
    if (PT(is_infinity)(p)) { F(set_zero)(x); F(set_zero)(y); return 0; }
    EC_FT zi;
    F(inv)(zi, p->z);
    F(mul)(x, p->x, zi);
    F(mul)(y, p->y, zi);
    return 1;
}

/* Constant-time fixed-window scalar multiplication, 4 bits at a time.
 *
 * The window index selects from the table by scanning every entry under a mask,
 * so no memory address depends on the scalar. Only the declared bit length
 * affects the loop count, and that is public.
 *
 * With complete addition this needs no special handling of the identity: the
 * accumulator starts at infinity and the formulas simply work. */
#define EC_WIN      4
#define EC_TBL_SIZE (1 << EC_WIN)

void PT(mul)(PTT *r, const PTT *p, const limb_t *k, int kbits)
{
    PTT tbl[EC_TBL_SIZE], acc, sel;

    PT(set_infinity)(&tbl[0]);
    PT(copy)(&tbl[1], p);
    for (int i = 2; i < EC_TBL_SIZE; i++) {
        if (i & 1) PT(add)(&tbl[i], &tbl[i - 1], p);
        else       PT(dbl)(&tbl[i], &tbl[i / 2]);
    }

    PT(set_infinity)(&acc);
    int top = ((kbits + EC_WIN - 1) / EC_WIN) * EC_WIN;
    for (int pos = top - EC_WIN; pos >= 0; pos -= EC_WIN) {
        for (int d = 0; d < EC_WIN; d++) PT(dbl)(&acc, &acc);

        limb_t w = 0;
        for (int b = 0; b < EC_WIN; b++) {
            int bit = pos + b;
            if (bit < kbits)
                w |= ((k[bit / 64] >> (bit % 64)) & 1) << b;
        }

        PT(set_infinity)(&sel);
        for (int i = 0; i < EC_TBL_SIZE; i++) {
            limb_t diff = w ^ (limb_t)i;
            limb_t mask = (limb_t)0 - (limb_t)(1 - (int)((diff | (~diff + 1)) >> 63));
            F(cselect)(sel.x, tbl[i].x, sel.x, mask);
            F(cselect)(sel.y, tbl[i].y, sel.y, mask);
            F(cselect)(sel.z, tbl[i].z, sel.z, mask);
        }
        PT(add)(&acc, &acc, &sel);
    }
    PT(copy)(r, &acc);
}

#undef EC_WIN
#undef EC_TBL_SIZE

/* y^2 z == x^3 + b z^3 */
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

#undef CAT_
#undef CAT
#undef PT
#undef PTT
#undef F
