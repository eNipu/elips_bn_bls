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

void ep_phi(ep_t *r, const ep_t *p)
{
    /* Safe in projective form: x = X/Z, and scaling X by beta scales x by beta
     * while leaving Z alone. Infinity has X = 0 and stays infinity.
     *
     * beta is the cube root of unity in Fp that makes phi act as [lambda] on
     * G1, and which lambda that is differs by family: -x^2 on BLS12, where phi
     * also drives the subgroup test, and 36x^3 + 18x^2 + 6x + 1 on BN, where it
     * only drives GLV. tools/reference/glv_ref.py picks it and gen_params.py
     * re-asserts it on a real point. */
    fp_mul(r->x, p->x, EP_BETA);
    fp_copy(r->y, p->y);
    fp_copy(r->z, p->z);
}

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

#include "arith/glv_scalar.h"

/* ---- the two-digit ladder ------------------------------------------------
 *
 * Every two-dimensional split here ends the same way: accumulate
 * [d0]P0 + [d1]P1 with both digits the same public length. Written once and
 * instantiated for E(Fp) and E'(Fp2).
 *
 * Two bits of each digit per step, not one. A bitwise interleaved ladder costs
 * one doubling AND one addition per bit, so halving the digit length halves the
 * doublings but leaves the additions where ep_mul's width-4 window already had
 * them -- and additions are the expensive half. On BN-462, ep_add is 3.49 us
 * against ep_dbl at 2.23 us, so a 231-step bitwise ladder is 231(D + A) = 1323
 * us against ep_mul's 462 D + 115 A = 1433 us. That is a 1.07x "speedup" for a
 * decomposition that halves the work, which is most of the gain thrown away.
 *
 * Taking two bits of each digit at a time indexes a 16-entry table by
 * (d1 << 2) | d0 and drops the additions to one per two bits: 231 D + 115 A =
 * 918 us by the same model. The table costs 14 point operations to build, paid
 * back several times over by the 115 additions it removes.
 *
 * Three is not better. A width-3 window needs 64 entries, and both the table
 * build and the masked scan grow as 2^(2w) while the additions saved grow only
 * as 1/w.
 *
 * Constant time: the loop bound is a public digit length, the digits are
 * extracted with fixed shifts, and every table entry is read under a mask on
 * every step. Two-bit digits never straddle a limb, since the position is
 * always even. */
#define GLV_DEF_LADDER2(SUF, PTT, PFX, FPFX)                                   \
static void glv_ladder2_##SUF(PTT *r, const PTT *P0, const PTT *P1,            \
                              const limb_t *d0, const limb_t *d1, int top)     \
{                                                                              \
    PTT tbl[16], acc, sel;                                                     \
                                                                               \
    /* tbl[j] = [j & 3] P0 + [j >> 2] P1. The addition law is complete, so the \
     * entries built on the identity need no special case. */                  \
    PFX##_set_infinity(&tbl[0]);                                               \
    PFX##_copy(&tbl[1], P0);                                                   \
    PFX##_dbl(&tbl[2], P0);                                                    \
    PFX##_add(&tbl[3], &tbl[2], P0);                                           \
    for (int j = 4; j < 16; j++) PFX##_add(&tbl[j], &tbl[j - 4], P1);          \
                                                                               \
    PFX##_set_infinity(&acc);                                                  \
    for (int pos = top - 2; pos >= 0; pos -= 2) {                              \
        PFX##_dbl(&acc, &acc);                                                 \
        PFX##_dbl(&acc, &acc);                                                 \
                                                                               \
        limb_t w = ((d0[pos / 64] >> (pos % 64)) & 3)                          \
                 | (((d1[pos / 64] >> (pos % 64)) & 3) << 2);                  \
                                                                               \
        PFX##_set_infinity(&sel);                                              \
        for (int j = 0; j < 16; j++) {                                         \
            limb_t diff = w ^ (limb_t)j;                                       \
            limb_t mask = (limb_t)0 -                                          \
                (limb_t)(1 - (int)((diff | (~diff + 1)) >> 63));               \
            FPFX##_cselect(sel.x, tbl[j].x, sel.x, mask);                      \
            FPFX##_cselect(sel.y, tbl[j].y, sel.y, mask);                      \
            FPFX##_cselect(sel.z, tbl[j].z, sel.z, mask);                      \
        }                                                                      \
        PFX##_add(&acc, &acc, &sel);                                           \
    }                                                                          \
    PFX##_copy(r, &acc);                                                       \
}

GLV_DEF_LADDER2(ep,  ep_t,  ep,  fp)
#ifdef ELIPS_FAMILY_BN
/* Only BN has a two-dimensional G2 split. On BLS12, psi acts as [x] and the G2
 * scalar has four digits, which a 2-bit window would index with a 256-entry
 * table -- more point operations to build than the whole ladder saves. */
GLV_DEF_LADDER2(ep2, ep2_t, ep2, fp2)
#endif

/* Digit length rounded up to a whole number of two-bit windows. */
#define GLV_TOP(bits) ((((bits) + 1 + 1) / 2) * 2)

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
    glv_divrem(quo, tmp, kbuf, inbits, ord, GLV_W);              /* tmp = k mod r    */
    glv_divrem(e[1], e[0], tmp, ELIPS_ORDER_BITS, u2, GLV_W);    /* tmp = e0 + e1 x^2 */

    /* P[0] = P, P[1] = [x^2]P = -phi(P). x^2 is positive whatever the sign of
     * x, so there is no sign case here as there is on G2. */
    ep_t P[2];
    ep_copy(&P[0], p);
    ep_phi(&P[1], p);
    ep_neg(&P[1], &P[1]);

    glv_ladder2_ep(r, &P[0], &P[1], e[0], e[1], GLV_TOP(ELIPS_ABSX2_BITS));
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
    glv_divrem(quo, tmp, kbuf, inbits, ord, GLV_W);   /* tmp = k mod r */

    int bits = ELIPS_ORDER_BITS;
    for (int i = 0; i < 3; i++) {
        glv_divrem(quo, d[i], tmp, bits, u, GLV_W);
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
    glv_divrem(quo, tmp, kbuf, inbits, ord, GLV_W);              /* tmp = k mod r      */
    glv_divrem(e[1], e[0], tmp, ELIPS_ORDER_BITS, b, GLV_W);     /* tmp = e0 + e1 6x^2 */

    ep2_t P[2];
    ep2_copy(&P[0], q);
    ep2_psi(&P[1], q);

    glv_ladder2_ep2(r, &P[0], &P[1], e[0], e[1], GLV_TOP(ELIPS_6XSQ_BITS));
}

/* Two-dimensional GLV on G1, for BN.
 *
 * phi(x, y) = (beta x, y) acts on G1 as [lambda] with
 * lambda = 36x^3 + 18x^2 + 6x + 1. That is 348 bits against a 462-bit order,
 * so unlike every other split in this file it is NOT short enough to use as a
 * base: there is no B with k = e0 + e1 B and both digits near sqrt(r). BN G1
 * needs the general GLV machinery, a reduced lattice basis and Babai rounding.
 *
 * tools/reference/glv_ref.py reduces the lattice; gen_params.py re-runs the
 * reduction and checks the closed forms against it before emitting them. With
 * a = 2x + 1, c = 6x^2 + 2x and d = a + c, the reduced basis is (-a, c) and
 * (-d, -a), its determinant is exactly r, and rounding (k, 0) onto the lattice
 * gives
 *
 *     q1 = round(k a / r)      q2 = round(k c / r)
 *     k1 = k - q1 a - q2 d     k2 = q1 c - q2 a
 *
 * with k = k1 + k2 lambda mod r. Both digits are bounded by half a basis
 * vector, ELIPS_GLV_BITS = 230 on BN-462, so the joint ladder is 231 steps
 * against 462.
 *
 * Unlike the other three routines, the digits here are SIGNED. They are kept
 * in two's complement, their signs are extracted as masks, and the two base
 * points are conditionally negated by masked select, so no branch sees them.
 *
 * The rounded division is floor((2N + r) / 2r), which is exactly round(N/r)
 * for the non-negative N here and needs no separate comparison.
 *
 * PRECONDITION: p must be in G1. On BN that is nearly free -- #E(Fp) = r is
 * prime, so every point on the curve except infinity is already in G1, which
 * is why ep_in_subgroup is a curve equation there and nothing more. */
void ep_mul_glv(ep_t *r, const ep_t *p, const limb_t *k, int kbits)
{
    limb_t kbuf[GLV_W], ord_n[GLV_W], quo_n[GLV_W], red_n[GLV_W];
    limb_t ord[GLV_WW], two_r[GLV_WW], kred[GLV_WW];
    limb_t A[GLV_WW], C[GLV_WW], D[GLV_WW];
    limb_t num[GLV_WW], tmp[GLV_WW], q1[GLV_WW], q2[GLV_WW];
    limb_t e1[GLV_WW], e2[GLV_WW];

    memset(kbuf, 0, sizeof kbuf);
    memset(ord_n, 0, sizeof ord_n);
    for (int i = 0; i < GLV_W && i < (kbits + 63) / 64; i++) kbuf[i] = k[i];
    for (int i = 0; i < GLV_W && i < (ELIPS_ORDER_BITS + 63) / 64; i++)
        ord_n[i] = ELIPS_ORDER[i];

    /* k mod r, at the narrow width: k never exceeds it and the divider costs
     * one pass per bit per limb, so there is nothing to gain from widening. */
    int inbits = kbits > ELIPS_ORDER_BITS ? kbits : ELIPS_ORDER_BITS;
    glv_divrem(quo_n, red_n, kbuf, inbits, ord_n, GLV_W);

    memset(ord,  0, sizeof ord);
    memset(kred, 0, sizeof kred);
    memset(A,    0, sizeof A);
    memset(C,    0, sizeof C);
    memset(D,    0, sizeof D);
    for (int i = 0; i < GLV_W; i++) { ord[i] = ord_n[i]; kred[i] = red_n[i]; }
    for (int i = 0; i < GLV_WW && i < (ELIPS_GLV_A_BITS + 63) / 64; i++) A[i] = ELIPS_GLV_A[i];
    for (int i = 0; i < GLV_WW && i < (ELIPS_GLV_C_BITS + 63) / 64; i++) C[i] = ELIPS_GLV_C[i];
    for (int i = 0; i < GLV_WW && i < (ELIPS_GLV_D_BITS + 63) / 64; i++) D[i] = ELIPS_GLV_D[i];
    (void)glv_add(two_r, ord, ord, GLV_WW);

    /* q1 = round(k a / r), q2 = round(k c / r). The numerator bounds are public
     * -- k < r and the basis entries are constants -- so the loop counts are
     * too. */
    glv_mul(num, kred, A, GLV_WW);
    (void)glv_add(num, num, num, GLV_WW);
    (void)glv_add(num, num, ord, GLV_WW);
    glv_divrem(q1, tmp, num, ELIPS_ORDER_BITS + ELIPS_GLV_A_BITS + 1, two_r, GLV_WW);

    glv_mul(num, kred, C, GLV_WW);
    (void)glv_add(num, num, num, GLV_WW);
    (void)glv_add(num, num, ord, GLV_WW);
    glv_divrem(q2, tmp, num, ELIPS_ORDER_BITS + ELIPS_GLV_C_BITS + 1, two_r, GLV_WW);

    /* e1 = k - q1 a - q2 d,  e2 = q1 c - q2 a, both signed. */
    glv_mul(tmp, q1, A, GLV_WW);
    (void)glv_sub(e1, kred, tmp, GLV_WW);
    glv_mul(tmp, q2, D, GLV_WW);
    (void)glv_sub(e1, e1, tmp, GLV_WW);

    glv_mul(e2, q1, C, GLV_WW);
    glv_mul(tmp, q2, A, GLV_WW);
    (void)glv_sub(e2, e2, tmp, GLV_WW);

    limb_t s1 = (limb_t)0 - (e1[GLV_WW - 1] >> 63);
    limb_t s2 = (limb_t)0 - (e2[GLV_WW - 1] >> 63);
    glv_cneg(e1, e1, s1, GLV_WW);
    glv_cneg(e2, e2, s2, GLV_WW);

    /* P[0] = [sign(e1)] P, P[1] = [sign(e2)] phi(P). ep_neg touches only y, so
     * one masked select per point applies the sign without a branch. */
    ep_t P[2];
    fp_t ny;
    ep_copy(&P[0], p);
    fp_neg(ny, P[0].y);
    fp_cselect(P[0].y, ny, P[0].y, s1);
    ep_phi(&P[1], p);
    fp_neg(ny, P[1].y);
    fp_cselect(P[1].y, ny, P[1].y, s2);

    glv_ladder2_ep(r, &P[0], &P[1], e1, e2, GLV_TOP(ELIPS_GLV_BITS));
}

#endif  /* ELIPS_FAMILY_BN */

#undef GLV_TOP
#undef GLV_DEF_LADDER2
