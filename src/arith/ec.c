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

#define EC_PT ep2
#define EC_F  fp2
#define EC_FT fp2_t
#include "arith/ec_tmpl.h"
#undef EC_PT
#undef EC_F
#undef EC_FT
