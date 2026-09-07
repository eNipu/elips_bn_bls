/*
 * Instantiates the Jacobian group law for both curves from one template.
 * See src/arith/ec_tmpl.h for the formulas.
 */
#include "elips/ec.h"

/* --- curve constants ------------------------------------------------------
 * Both curves are y^2 = x^3 + b with b = +-4. The sextic twist that carries G2
 * is y^2 = x^3 + b*xi, which follows from the untwisting map the legacy code
 * uses: x = x'/w^2, y = y'/w^3 with xi = w^6. */

void ep_curve_b(fp_t b)
{
    limb_t v[FP_LIMBS];
    memset(v, 0, sizeof v);
    v[0] = (limb_t)ELIPS_CURVE_B;
    fp_from_limbs(b, v);
#if ELIPS_CURVE_B_SIGN < 0
    fp_neg(b, b);
#endif
}

void ep2_curve_b(fp2_t b)
{
    fp_t bb;
    ep_curve_b(bb);
    fp_copy(b[0], bb);
    fp_set_zero(b[1]);
    fp2_mul_xi(b, b);
}

#define EC_PT ep
#define EC_F  fp
#define EC_FT fp_t
#include "arith/ec_tmpl.h"
#undef EC_PT
#undef EC_F
#undef EC_FT

#ifdef ELIPS_FAMILY_BLS12
void ep_phi(ep_t *r, const ep_t *p)
{
    /* Safe in projective form: x = X/Z, and scaling X by beta scales x by beta
     * while leaving Z alone. Infinity has X = 0 and stays infinity. */
    fp_mul(r->x, p->x, EP_BETA);
    fp_copy(r->y, p->y);
    fp_copy(r->z, p->z);
}
#endif

void ep2_psi(ep2_t *r, const ep2_t *p)
{
    /* psi(x, y) = (conj(x) * gamma^-2, conj(y) * gamma^-3) on affine
     * coordinates.
     *
     * In the homogeneous projective form used here, x = X/Z and y = Y/Z, and
     * conjugation is a field homomorphism, so conjugating all three coordinates
     * already gives (conj(x), conj(y)); scaling X and Y by the two constants
     * finishes it. No special case and no normalisation.
     *
     * (An earlier version of this comment described Jacobian coordinates with
     * x = X/Z^2. The code was and is right for the homogeneous form the curve
     * layer moved to in Phase 4; the comment was not.) */
    fp2_conj(r->x, p->x); fp2_mul(r->x, r->x, PSI_X);
    fp2_conj(r->y, p->y); fp2_mul(r->y, r->y, PSI_Y);
    fp2_conj(r->z, p->z);
}

#define EC_PT ep2
#define EC_F  fp2
#define EC_FT fp2_t
#include "arith/ec_tmpl.h"
#undef EC_PT
#undef EC_F
#undef EC_FT

/* ---- constant-time division, for the GLV decomposition -------------------
 *
 * Restoring division one bit at a time, with the conditional subtraction done
 * by mask rather than by branch. The loop count comes from a public bit length,
 * never from the value, so the timing reveals nothing about the scalar.
 *
 * This replaces an earlier version that used GMP's mpz_tdiv_qr. That was
 * correct but data dependent, which made the whole GLV path unusable on secret
 * scalars and forced it to be opt-in. It is now the default. */

#define GLV_W ((int)FP_LIMBS)

static limb_t glv_sub(limb_t *r, const limb_t *a, const limb_t *b)
{
    limb_t borrow = 0;
    for (int i = 0; i < GLV_W; i++) {
        limb_t ai = a[i], bi = b[i];
        limb_t d  = ai - bi;
        limb_t b1 = (ai < bi);
        limb_t d2 = d - borrow;
        limb_t b2 = (d < borrow);
        r[i] = d2;
        borrow = b1 | b2;
    }
    return borrow;
}

static void glv_shl1(limb_t *a, limb_t in)
{
    limb_t carry = in;
    for (int i = 0; i < GLV_W; i++) {
        limb_t next = a[i] >> 63;
        a[i] = (a[i] << 1) | carry;
        carry = next;
    }
}

/* q = k / d, rem = k mod d. kbits is public. */
static void glv_divrem(limb_t *q, limb_t *rem, const limb_t *k, int kbits,
                       const limb_t *d)
{
    limb_t t[GLV_W];
    memset(q, 0, sizeof(limb_t) * GLV_W);
    memset(rem, 0, sizeof(limb_t) * GLV_W);

    for (int i = kbits - 1; i >= 0; i--) {
        glv_shl1(rem, (k[i / 64] >> (i % 64)) & 1);
        limb_t borrow = glv_sub(t, rem, d);
        limb_t mask   = (limb_t)0 - (1 - borrow);      /* all ones if rem >= d */
        for (int j = 0; j < GLV_W; j++)
            rem[j] = (t[j] & mask) | (rem[j] & ~mask);
        q[i / 64] |= (mask & 1) << (i % 64);
    }
}

#ifdef ELIPS_FAMILY_BLS12

/* Two-dimensional GLV on G1.
 *
 * phi(x, y) = (beta x, y) acts on G1 as multiplication by -x^2, and |x^2| is
 * sqrt(r) to within a bit. So writing the scalar in base x^2,
 *
 *     k = e0 + e1 x^2,   [k]P = [e0]P + [e1][x^2]P = [e0]P + [e1](-phi(P))
 *
 * gives two digits of half length and halves the ladder. No lattice reduction
 * and no Babai rounding: one constant-time division by x^2 does the whole
 * decomposition, reusing glv_divrem.
 *
 * BN is excluded because its lambda is 348 bits against sqrt(r) = 231, so no
 * base-B split exists there; ep_mul stays the general routine for it.
 *
 * PRECONDITION, same as ep2_mul_glv: p must be in G1. phi acts as [-x^2] only
 * on that eigenspace. */
void ep_mul_glv(ep_t *r, const ep_t *p, const limb_t *k, int kbits)
{
    limb_t kbuf[GLV_W], ord[GLV_W], u2[GLV_W];
    limb_t e[2][GLV_W], quo[GLV_W], tmp[GLV_W];

    memset(kbuf, 0, sizeof kbuf);
    memset(ord,  0, sizeof ord);
    memset(u2,   0, sizeof u2);
    for (int i = 0; i < GLV_W && i < (kbits + 63) / 64; i++) kbuf[i] = k[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ORDER_BITS + 63) / 64; i++) ord[i] = ELIPS_ORDER[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ABSX2_BITS + 63) / 64; i++) u2[i] = ELIPS_ABSX2[i];

    int inbits = kbits > ELIPS_ORDER_BITS ? kbits : ELIPS_ORDER_BITS;
    glv_divrem(quo, tmp, kbuf, inbits, ord);              /* tmp = k mod r    */
    glv_divrem(e[1], e[0], tmp, ELIPS_ORDER_BITS, u2);    /* tmp = e0 + e1 x^2 */

    /* P[0] = P, P[1] = [x^2]P = -phi(P). x^2 is positive whatever the sign of
     * x, so there is no sign case here as there is on G2. */
    ep_t P[2];
    ep_copy(&P[0], p);
    ep_phi(&P[1], p);
    ep_neg(&P[1], &P[1]);

    ep_t tbl[4];
    ep_set_infinity(&tbl[0]);
    ep_copy(&tbl[1], &P[0]);
    ep_copy(&tbl[2], &P[1]);
    ep_add(&tbl[3], &P[0], &P[1]);

    ep_t acc, sel;
    ep_set_infinity(&acc);
    for (int i = ELIPS_ABSX2_BITS; i >= 0; i--) {
        ep_dbl(&acc, &acc);
        limb_t w = 0;
        for (int j = 0; j < 2; j++)
            w |= ((e[j][i / 64] >> (i % 64)) & 1) << j;

        ep_set_infinity(&sel);
        for (int j = 0; j < 4; j++) {
            limb_t diff = w ^ (limb_t)j;
            limb_t mask = (limb_t)0 - (limb_t)(1 - (int)((diff | (~diff + 1)) >> 63));
            fp_cselect(sel.x, tbl[j].x, sel.x, mask);
            fp_cselect(sel.y, tbl[j].y, sel.y, mask);
            fp_cselect(sel.z, tbl[j].z, sel.z, mask);
        }
        ep_add(&acc, &acc, &sel);
    }
    ep_copy(r, &acc);
}

/* Four-dimensional GLV on G2.
 *
 * psi acts on G2 as multiplication by p mod r, which on BLS12 reduces to the
 * mother parameter x -- only 64 to 77 bits against a 255 to 308 bit order. So
 * the scalar in base |x| has four short digits and the ladder is a quarter as
 * long.
 *
 * Constant time throughout: the decomposition above, and a table lookup that
 * scans every entry under a mask. */
void ep2_mul_glv(ep2_t *r, const ep2_t *q, const limb_t *k, int kbits)
{
    limb_t kbuf[GLV_W], ord[GLV_W], u[GLV_W];
    limb_t d[4][GLV_W], quo[GLV_W], tmp[GLV_W];

    memset(kbuf, 0, sizeof kbuf);
    memset(ord,  0, sizeof ord);
    memset(u,    0, sizeof u);
    for (int i = 0; i < GLV_W && i < (kbits + 63) / 64; i++) kbuf[i] = k[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ORDER_BITS + 63) / 64; i++) ord[i] = ELIPS_ORDER[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ABSX_BITS + 63) / 64; i++) u[i] = ELIPS_ABSX[i];

    /* reduce mod r, then peel off four base-|x| digits */
    int inbits = kbits > ELIPS_ORDER_BITS ? kbits : ELIPS_ORDER_BITS;
    glv_divrem(quo, tmp, kbuf, inbits, ord);          /* tmp = k mod r */

    int bits = ELIPS_ORDER_BITS;
    for (int i = 0; i < 3; i++) {
        glv_divrem(quo, d[i], tmp, bits, u);
        memcpy(tmp, quo, sizeof tmp);
        bits -= ELIPS_ABSX_BITS - 1;
        if (bits < ELIPS_ABSX_BITS) bits = ELIPS_ABSX_BITS;
    }
    memcpy(d[3], tmp, sizeof d[3]);

    /* base points: psi = [x], and x is negative here, so [|x|]Q = -psi(Q) */
    ep2_t P[4], t1;
    ep2_copy(&P[0], q);
    ep2_psi(&t1, q);
    ep2_copy(&P[1], &t1);
    ep2_psi(&P[2], &t1);
    ep2_psi(&P[3], &P[2]);
#if ELIPS_X_NEGATIVE
    ep2_neg(&P[1], &P[1]);
    ep2_neg(&P[3], &P[3]);
#endif

    ep2_t tbl[16];
    ep2_set_infinity(&tbl[0]);
    for (int j = 1; j < 16; j++) {
        int low = j & (-j), idx = 0;
        while ((1 << idx) != low) idx++;
        ep2_add(&tbl[j], &tbl[j ^ low], &P[idx]);
    }

    ep2_t acc, sel;
    ep2_set_infinity(&acc);
    for (int i = ELIPS_ABSX_BITS; i >= 0; i--) {
        ep2_dbl(&acc, &acc);
        limb_t w = 0;
        for (int j = 0; j < 4; j++)
            w |= ((d[j][i / 64] >> (i % 64)) & 1) << j;

        ep2_set_infinity(&sel);
        for (int j = 0; j < 16; j++) {
            limb_t diff = w ^ (limb_t)j;
            limb_t mask = (limb_t)0 - (limb_t)(1 - (int)((diff | (~diff + 1)) >> 63));
            fp2_cselect(sel.x, tbl[j].x, sel.x, mask);
            fp2_cselect(sel.y, tbl[j].y, sel.y, mask);
            fp2_cselect(sel.z, tbl[j].z, sel.z, mask);
        }
        ep2_add(&acc, &acc, &sel);
    }
    ep2_copy(r, &acc);
}

#endif  /* ELIPS_FAMILY_BLS12 */

#ifdef ELIPS_FAMILY_BN
/* Two-dimensional GLV on G2, for BN.
 *
 * psi acts on G2 as multiplication by p mod r, which on BN is 6x^2 -- 231 bits
 * against a 462-bit order, so exactly sqrt(r). Writing the scalar in base 6x^2,
 *
 *     k = e0 + e1 (6x^2),   [k]Q = [e0]Q + [e1] psi(Q)
 *
 * gives two half-length digits and halves the ladder. As on BLS12 G1 this needs
 * no lattice reduction: one constant-time division by 6x^2 does it, and
 * ELIPS_6XSQ already exists because the fast subgroup test uses the same
 * multiplier.
 *
 * Two dimensions, not four: BLS12 gets four because psi acts as [x] there and
 * x is a quarter of r, but on BN the smallest available multiplier is 6x^2.
 *
 * PRECONDITION: q must be in G2. psi acts as [6x^2] only on that eigenspace. */
void ep2_mul_glv(ep2_t *r, const ep2_t *q, const limb_t *k, int kbits)
{
    limb_t kbuf[GLV_W], ord[GLV_W], b[GLV_W];
    limb_t e[2][GLV_W], quo[GLV_W], tmp[GLV_W];

    memset(kbuf, 0, sizeof kbuf);
    memset(ord,  0, sizeof ord);
    memset(b,    0, sizeof b);
    for (int i = 0; i < GLV_W && i < (kbits + 63) / 64; i++) kbuf[i] = k[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ORDER_BITS + 63) / 64; i++) ord[i] = ELIPS_ORDER[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_6XSQ_BITS + 63) / 64; i++) b[i] = ELIPS_6XSQ[i];

    int inbits = kbits > ELIPS_ORDER_BITS ? kbits : ELIPS_ORDER_BITS;
    glv_divrem(quo, tmp, kbuf, inbits, ord);              /* tmp = k mod r      */
    glv_divrem(e[1], e[0], tmp, ELIPS_ORDER_BITS, b);     /* tmp = e0 + e1 6x^2 */

    ep2_t P[2];
    ep2_copy(&P[0], q);
    ep2_psi(&P[1], q);

    ep2_t tbl[4];
    ep2_set_infinity(&tbl[0]);
    ep2_copy(&tbl[1], &P[0]);
    ep2_copy(&tbl[2], &P[1]);
    ep2_add(&tbl[3], &P[0], &P[1]);

    ep2_t acc, sel;
    ep2_set_infinity(&acc);
    for (int i = ELIPS_6XSQ_BITS; i >= 0; i--) {
        ep2_dbl(&acc, &acc);
        limb_t w = 0;
        for (int j = 0; j < 2; j++)
            w |= ((e[j][i / 64] >> (i % 64)) & 1) << j;

        ep2_set_infinity(&sel);
        for (int j = 0; j < 4; j++) {
            limb_t diff = w ^ (limb_t)j;
            limb_t mask = (limb_t)0 - (limb_t)(1 - (int)((diff | (~diff + 1)) >> 63));
            fp2_cselect(sel.x, tbl[j].x, sel.x, mask);
            fp2_cselect(sel.y, tbl[j].y, sel.y, mask);
            fp2_cselect(sel.z, tbl[j].z, sel.z, mask);
        }
        ep2_add(&acc, &acc, &sel);
    }
    ep2_copy(r, &acc);
}
#endif  /* ELIPS_FAMILY_BN */

#undef GLV_W
