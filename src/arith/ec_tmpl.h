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

/* Constant-time fixed-window scalar multiplication, 4 bits at a time.
 *
 * The window index selects from a precomputed table by scanning every entry and
 * combining with a mask, so no memory address depends on the scalar. Only the
 * declared bit length affects the loop count, and that is public.
 *
 * This replaces double-and-add-always, which performed an addition on every bit
 * and threw half of them away. For a 256-bit scalar the window does 256
 * doublings and 64 additions where the old routine did 256 of each.
 *
 * Window 4 is a deliberate middle: the table is 16 points, and the masked scan
 * that keeps the lookup constant-time costs 16 selects per window, so a wider
 * window buys fewer additions but pays more scanning.
 * ponytail: 4 bits until a measurement says otherwise. */
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
    int top = ((kbits + EC_WIN - 1) / EC_WIN) * EC_WIN;   /* round up */
    for (int pos = top - EC_WIN; pos >= 0; pos -= EC_WIN) {
        for (int d = 0; d < EC_WIN; d++) PT(dbl)(&acc, &acc);

        /* extract the window without branching on its value */
        limb_t w = 0;
        for (int b = 0; b < EC_WIN; b++) {
            int bit = pos + b;
            if (bit < kbits)
                w |= ((k[bit / 64] >> (bit % 64)) & 1) << b;
        }

        /* masked linear scan: touch every entry, keep one */
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
