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

void ep2_psi(ep2_t *r, const ep2_t *p)
{
    /* psi(x, y) = (conj(x) * gamma^-2, conj(y) * gamma^-3) on affine
     * coordinates. In Jacobian form the z coordinate conjugates too, because
     * the map is applied coordinate-wise to (X : Y : Z) with x = X/Z^2. */
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

#ifdef ELIPS_FAMILY_BLS12
#include <gmp.h>

/* Four-dimensional GLV on G2. See the note in ec.h about the decomposition. */
void ep2_mul_glv(ep2_t *r, const ep2_t *q, const limb_t *k, int kbits)
{
    /* --- decompose k in base |x| -------------------------------------- */
    mpz_t K, U, R_, dig[4];
    mpz_inits(K, U, R_, dig[0], dig[1], dig[2], dig[3], NULL);
    mpz_import(K, (size_t)((kbits + 63) / 64), -1, sizeof(limb_t), 0, 0, k);
    mpz_import(U, (size_t)((ELIPS_ABSX_BITS + 63) / 64), -1, sizeof(limb_t), 0, 0, ELIPS_ABSX);
    mpz_import(R_, (size_t)((ELIPS_ORDER_BITS + 63) / 64), -1, sizeof(limb_t), 0, 0, ELIPS_ORDER);
    mpz_mod(K, K, R_);                     /* GLV needs k reduced */
    for (int i = 0; i < 4; i++) {
        mpz_tdiv_qr(K, dig[i], K, U);      /* dig[i] = K mod |x|, K /= |x| */
    }

    limb_t d[4][FP_LIMBS];
    for (int i = 0; i < 4; i++) {
        memset(d[i], 0, sizeof d[i]);
        mpz_export(d[i], NULL, -1, sizeof(limb_t), 0, 0, dig[i]);
    }
    mpz_clears(K, U, R_, dig[0], dig[1], dig[2], dig[3], NULL);

    /* --- the four base points ------------------------------------------
     * psi = [x]. With x negative, [|x|]Q = -psi(Q), so
     *   P0 = Q, P1 = [|x|]Q = -psi(Q), P2 = [|x|^2]Q = psi^2(Q),
     *   P3 = [|x|^3]Q = -psi^3(Q).
     * With x positive the signs simply do not alternate. */
    ep2_t P[4], t1;
    ep2_copy(&P[0], q);
    ep2_psi(&t1, q);                       /* psi(Q)   */
    ep2_copy(&P[1], &t1);
    ep2_psi(&P[2], &t1);                    /* psi^2(Q) */
    ep2_psi(&P[3], &P[2]);                  /* psi^3(Q) */
#if ELIPS_X_NEGATIVE
    ep2_neg(&P[1], &P[1]);
    ep2_neg(&P[3], &P[3]);
#endif

    /* --- joint table: every subset sum of the four base points --------- */
    ep2_t tbl[16];
    ep2_set_infinity(&tbl[0]);
    for (int j = 1; j < 16; j++) {
        int low = j & (-j);                 /* lowest set bit */
        int idx = 0; while ((1 << idx) != low) idx++;
        ep2_add(&tbl[j], &tbl[j ^ low], &P[idx]);
    }

    /* --- interleaved ladder, one bit of each digit per step ------------ */
    ep2_t acc, sel;
    ep2_set_infinity(&acc);
    for (int i = ELIPS_ABSX_BITS - 1; i >= 0; i--) {
        ep2_dbl(&acc, &acc);
        limb_t w = 0;
        for (int j = 0; j < 4; j++)
            w |= ((d[j][i / 64] >> (i % 64)) & 1) << j;

        ep2_set_infinity(&sel);
        for (int j = 0; j < 16; j++) {      /* masked scan, no secret index */
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
#endif
