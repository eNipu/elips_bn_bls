/*
 * The scalar multiplications and the representation-agnostic point routines,
 * written once and instantiated for both E(Fp) and E'(Fp2). The original
 * library kept two hand-copied versions of every curve routine and they had
 * already drifted apart; generating both from one text removes that.
 *
 * The group law itself is in ec_jacobian.h. Nothing in THIS file knows which
 * coordinates are in use: set_infinity, is_infinity, copy, neg and from_affine
 * mean the same thing in any of them, and both ladders only ever call add, dbl
 * and cselect.
 *
 * Both groups were homogeneous projective (x = X/Z) until issue #50 moved the
 * twist to Jacobian and issue #51 moved G1 after it, so the homogeneous RCB
 * formulas that used to sit beside them are gone. What that cost is argued at
 * the top of ec_jacobian.h: those formulas were complete in one expression and
 * these are not, which is why test/ec_group_test.c now checks both groups
 * against an affine oracle rather than one.
 *
 * Expects EC_PT, EC_F and EC_FT from the includer.
 */

#define CAT_(a,b) a##b
#define CAT(a,b)  CAT_(a,b)
#define PT(name)  CAT(EC_PT, CAT(_, name))
#define PTT       CAT(EC_PT, _t)
#define F(name)   CAT(EC_F,  CAT(_, name))

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

/* The group law: add, dbl, to_affine, on_curve, eq. Its own file rather than
 * an inlined block, because it is a formula set carrying its own exception
 * argument and it is the code that decides whether an attacker-supplied point
 * is accepted, so it wants to be readable as a unit. */
#include "arith/ec_jacobian.h"

/* [k]P where k is a PUBLIC CONSTANT of the curve, not a caller's scalar.
 *
 * The two subgroup tests multiply by |x| and by 6x^2, which are curve
 * parameters compiled into fp_params.h. Nothing about them is secret, so this
 * expands them into non-adjacent form and skips the zero digits.
 *
 * NAF rather than plain binary, and the difference is not small. Binary was
 * tried first and made two of the three curves SLOWER, because the seeds are
 * sparse in signed-digit form and dense in bits:
 *
 *      constant                  bits   popcount   NAF
 *      BLS12-381  |x|              64          6     6
 *      BLS12-461  |x|              77         43     3
 *      BN-462     |x|             115        101     4
 *      BN-462     6x^2            231        106    18
 *
 * BLS12-461 pays 42 additions in binary and 2 in NAF. These seeds are chosen
 * near powers of two, which is exactly the shape binary represents worst.
 *
 * Still constant time in P, which is the property that matters here. The
 * control flow depends on k alone, so the sequence of operations is identical
 * for every point, and P is what an attacker supplies. Do NOT call this with a
 * scalar derived from a secret: use PT(mul), which is constant time in both.
 *
 * The unified addition again, and here its cost is not an optimisation to give
 * up. The caller is validating a point it does not trust, so every exceptional
 * case has to work rather than be argued away. PT(add) is not complete in the
 * RCB sense -- it selects between an addition and a doubling under a mask --
 * but it is correct for every input pair, which is the property this needs. */
#define PUBCONST_LIMBS (((ELIPS_ORDER_BITS + 63) / 64) + 1)
#define PUBCONST_BITS  (64 * PUBCONST_LIMBS)

void PT(mul_pubconst)(PTT *r, const PTT *p, const limb_t *k, int kbits)
{
    signed char naf[PUBCONST_BITS + 1];
    limb_t t[PUBCONST_LIMBS];
    PTT acc, np;
    int n = 0, i, nz;

    /* Nothing here is sized for a caller's scalar, and nothing should reach
     * it. Fall back rather than overrun if one ever does. */
    if (kbits <= 0 || kbits > PUBCONST_BITS) { PT(mul)(r, p, k, kbits > 0 ? kbits : 0); return; }

    for (i = 0; i < PUBCONST_LIMBS; i++)
        t[i] = (i < (kbits + 63) / 64) ? k[i] : 0;
    if (kbits % 64) t[(kbits - 1) / 64] &= (limb_t)-1 >> (64 - kbits % 64);

    for (;;) {
        for (nz = 0, i = 0; i < PUBCONST_LIMBS; i++) nz |= (t[i] != 0);
        if (!nz) break;
        if (t[0] & 1) {
            /* digit is +1 when the next bit up is 0, else -1 */
            int z = 2 - (int)(t[0] & 3);
            naf[n] = (signed char)z;
            if (z == 1) {                       /* t -= 1 */
                for (i = 0; i < PUBCONST_LIMBS && t[i]-- == 0; i++) { }
            } else {                            /* t += 1 */
                for (i = 0; i < PUBCONST_LIMBS && ++t[i] == 0; i++) { }
            }
        } else {
            naf[n] = 0;
        }
        for (i = 0; i + 1 < PUBCONST_LIMBS; i++)
            t[i] = (t[i] >> 1) | (t[i + 1] << 63);
        t[PUBCONST_LIMBS - 1] >>= 1;
        n++;
    }

    if (n == 0) { PT(set_infinity)(r); return; }

    PT(neg)(&np, p);
    PT(copy)(&acc, naf[n - 1] > 0 ? p : &np);
    for (i = n - 2; i >= 0; i--) {
        PT(dbl)(&acc, &acc);
        if      (naf[i] > 0) PT(add)(&acc, &acc, p);
        else if (naf[i] < 0) PT(add)(&acc, &acc, &np);
    }
    PT(copy)(r, &acc);
}

/* Constant-time fixed-window scalar multiplication, 4 bits at a time.
 *
 * The window index selects from the table by scanning every entry under a mask,
 * so no memory address depends on the scalar. Only the declared bit length
 * affects the loop count, and that is public.
 *
 * The accumulator starts at infinity and needs no special handling, because
 * PT(add) returns the other operand when either one is the identity. */
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

#undef CAT_
#undef CAT
#undef PT
#undef PTT
#undef F
