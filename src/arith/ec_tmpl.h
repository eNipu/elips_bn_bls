/*
 * Jacobian group law, written once and instantiated for both E(Fp) and
 * E'(Fp2). The original library kept two hand-copied versions of every curve
 * routine and they had already drifted apart; generating both from one text
 * removes that failure mode entirely.
 *
 * Expects EC_PT, EC_F (field prefix) and EC_FT (field element type) to be
 * defined by the includer.
 */

#define CAT_(a,b) a##b
#define CAT(a,b)  CAT_(a,b)
#define PT(name)  CAT(EC_PT, CAT(_, name))
#define PTT       CAT(EC_PT, _t)
#define F(name)   CAT(EC_F,  CAT(_, name))

void PT(set_infinity)(PTT *r)
{
    F(set_zero)(r->x); F(set_zero)(r->y); F(set_zero)(r->z);
    /* X = Y = Z = 0 marks infinity; only Z is tested. */
}

int PT(is_infinity)(const PTT *p) { return F(is_zero)(p->z); }

void PT(copy)(PTT *r, const PTT *p)
{ F(copy)(r->x, p->x); F(copy)(r->y, p->y); F(copy)(r->z, p->z); }

void PT(neg)(PTT *r, const PTT *p)
{ F(copy)(r->x, p->x); F(neg)(r->y, p->y); F(copy)(r->z, p->z); }

void PT(from_affine)(PTT *r, const EC_FT x, const EC_FT y)
{ F(copy)(r->x, x); F(copy)(r->y, y); F(set_one)(r->z); }

/* dbl-2009-l for a = 0: 2M + 5S. */
void PT(dbl)(PTT *r, const PTT *p)
{
    EC_FT A, B, C, D, E, FF, t;

    F(sqr)(A, p->x);                     /* A = X^2            */
    F(sqr)(B, p->y);                     /* B = Y^2            */
    F(sqr)(C, B);                        /* C = B^2            */

    F(add)(D, p->x, B);                  /* D = 2((X+B)^2-A-C) */
    F(sqr)(D, D);
    F(sub)(D, D, A);
    F(sub)(D, D, C);
    F(add)(D, D, D);

    F(add)(E, A, A);                     /* E = 3A             */
    F(add)(E, E, A);
    F(sqr)(FF, E);                       /* F = E^2            */

    F(add)(t, D, D);                     /* X3 = F - 2D        */
    F(sub)(t, FF, t);

    F(mul)(r->z, p->y, p->z);            /* Z3 = 2 Y Z         */
    F(add)(r->z, r->z, r->z);

    F(sub)(D, D, t);                     /* Y3 = E(D-X3) - 8C  */
    F(mul)(D, E, D);
    F(add)(C, C, C); F(add)(C, C, C); F(add)(C, C, C);
    F(sub)(r->y, D, C);
    F(copy)(r->x, t);

    /* Doubling the point at infinity, or a point of order two, must give
     * infinity. Z3 = 2YZ already gives that, so nothing extra is needed. */
}

/* add-2007-bl, with the degenerate cases resolved by masked select rather than
 * by branching, so the timing does not reveal which case occurred. */
void PT(add)(PTT *r, const PTT *p, const PTT *q)
{
    EC_FT Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, rr, V, t;
    PTT out, dbl_res;

    F(sqr)(Z1Z1, p->z);
    F(sqr)(Z2Z2, q->z);
    F(mul)(U1, p->x, Z2Z2);
    F(mul)(U2, q->x, Z1Z1);
    F(mul)(S1, p->y, q->z); F(mul)(S1, S1, Z2Z2);
    F(mul)(S2, q->y, p->z); F(mul)(S2, S2, Z1Z1);

    F(sub)(H, U2, U1);
    F(sub)(rr, S2, S1);

    int h_zero = F(is_zero)(H);
    int r_zero = F(is_zero)(rr);

    F(add)(I, H, H); F(sqr)(I, I);
    F(mul)(J, H, I);
    F(add)(rr, rr, rr);
    F(mul)(V, U1, I);

    F(sqr)(out.x, rr);
    F(sub)(out.x, out.x, J);
    F(add)(t, V, V);
    F(sub)(out.x, out.x, t);

    F(sub)(t, V, out.x);
    F(mul)(t, rr, t);
    F(mul)(S1, S1, J);
    F(add)(S1, S1, S1);
    F(sub)(out.y, t, S1);

    F(add)(t, p->z, q->z);
    F(sqr)(t, t);
    F(sub)(t, t, Z1Z1);
    F(sub)(t, t, Z2Z2);
    F(mul)(out.z, t, H);

    /* Same x and same y: this is a doubling, and the formula above produced
     * zeros. Same x, different y: the sum is infinity. */
    PT(dbl)(&dbl_res, p);
    limb_t m_dbl = (limb_t)0 - (limb_t)(h_zero & r_zero);
    F(cselect)(out.x, dbl_res.x, out.x, m_dbl);
    F(cselect)(out.y, dbl_res.y, out.y, m_dbl);
    F(cselect)(out.z, dbl_res.z, out.z, m_dbl);

    PTT inf; PT(set_infinity)(&inf);
    limb_t m_inf = (limb_t)0 - (limb_t)(h_zero & (1 - r_zero));
    F(cselect)(out.x, inf.x, out.x, m_inf);
    F(cselect)(out.y, inf.y, out.y, m_inf);
    F(cselect)(out.z, inf.z, out.z, m_inf);

    /* Either operand at infinity: return the other one. */
    limb_t m_p_inf = (limb_t)0 - (limb_t)PT(is_infinity)(p);
    limb_t m_q_inf = (limb_t)0 - (limb_t)PT(is_infinity)(q);
    F(cselect)(out.x, p->x, out.x, m_q_inf);
    F(cselect)(out.y, p->y, out.y, m_q_inf);
    F(cselect)(out.z, p->z, out.z, m_q_inf);
    F(cselect)(out.x, q->x, out.x, m_p_inf);
    F(cselect)(out.y, q->y, out.y, m_p_inf);
    F(cselect)(out.z, q->z, out.z, m_p_inf);

    PT(copy)(r, &out);
}

/* Addition without the degenerate-case fixup.
 *
 * Still branch-free and therefore still constant time, but INCORRECT when
 * p == q, p == -q, or either is infinity. The complete PT(add) above pays for
 * those cases by computing a full doubling on every call and selecting, which
 * roughly doubles its cost -- measured, the complete version was no faster than
 * the affine code it replaced, while this one is.
 *
 * Only call it where the inputs are structurally distinct. The Miller loop
 * qualifies: T is a running multiple of Q and never coincides with it. */
void PT(add_generic)(PTT *r, const PTT *p, const PTT *q)
{
    EC_FT Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, rr, V, t;
    PTT out;

    F(sqr)(Z1Z1, p->z);
    F(sqr)(Z2Z2, q->z);
    F(mul)(U1, p->x, Z2Z2);
    F(mul)(U2, q->x, Z1Z1);
    F(mul)(S1, p->y, q->z); F(mul)(S1, S1, Z2Z2);
    F(mul)(S2, q->y, p->z); F(mul)(S2, S2, Z1Z1);

    F(sub)(H, U2, U1);
    F(sub)(rr, S2, S1);
    F(add)(I, H, H); F(sqr)(I, I);
    F(mul)(J, H, I);
    F(add)(rr, rr, rr);
    F(mul)(V, U1, I);

    F(sqr)(out.x, rr);
    F(sub)(out.x, out.x, J);
    F(add)(t, V, V);
    F(sub)(out.x, out.x, t);

    F(sub)(t, V, out.x);
    F(mul)(t, rr, t);
    F(mul)(S1, S1, J);
    F(add)(S1, S1, S1);
    F(sub)(out.y, t, S1);

    F(add)(t, p->z, q->z);
    F(sqr)(t, t);
    F(sub)(t, t, Z1Z1);
    F(sub)(t, t, Z2Z2);
    F(mul)(out.z, t, H);

    PT(copy)(r, &out);
}

int PT(to_affine)(EC_FT x, EC_FT y, const PTT *p)
{
    if (PT(is_infinity)(p)) { F(set_zero)(x); F(set_zero)(y); return 0; }
    EC_FT zi, zi2, zi3;
    F(inv)(zi, p->z);
    F(sqr)(zi2, zi);
    F(mul)(zi3, zi2, zi);
    F(mul)(x, p->x, zi2);
    F(mul)(y, p->y, zi3);
    return 1;
}

/* Double-and-add-always: the addition happens on every bit and its result is
 * discarded by a masked select when the bit is zero, so the sequence of field
 * operations does not depend on the scalar. Costs twice a plain ladder; window
 * methods are Phase 4's job. */
void PT(mul)(PTT *r, const PTT *p, const limb_t *k, int kbits)
{
    PTT acc, sum;
    PT(set_infinity)(&acc);
    for (int i = kbits - 1; i >= 0; i--) {
        PT(dbl)(&acc, &acc);
        PT(add)(&sum, &acc, p);
        limb_t bit  = (k[i / 64] >> (i % 64)) & 1;
        limb_t mask = (limb_t)0 - bit;
        F(cselect)(acc.x, sum.x, acc.x, mask);
        F(cselect)(acc.y, sum.y, acc.y, mask);
        F(cselect)(acc.z, sum.z, acc.z, mask);
    }
    PT(copy)(r, &acc);
}

/* y^2 == x^3 + b, in Jacobian form: Y^2 == X^3 + b Z^6 */
int PT(on_curve)(const PTT *p)
{
    if (PT(is_infinity)(p)) return 1;
    EC_FT lhs, rhs, z2, z6, b;
    PT(curve_b)(b);
    F(sqr)(lhs, p->y);
    F(sqr)(rhs, p->x); F(mul)(rhs, rhs, p->x);
    F(sqr)(z2, p->z); F(sqr)(z6, z2); F(mul)(z6, z6, z2);
    F(mul)(z6, z6, b);
    F(add)(rhs, rhs, z6);
    return F(eq)(lhs, rhs);
}

#undef CAT_
#undef CAT
#undef PT
#undef PTT
#undef F
