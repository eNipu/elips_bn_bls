/*
 * Does the projective group law agree with affine arithmetic, everywhere?
 *
 * Both groups are in Jacobian coordinates (x = X/Z^2, y = Y/Z^3), where the
 * doubling is 2M + 5S and the addition is the unified add-or-double form.
 * Neither is complete the way the RCB formulas were, so both carry an
 * exception argument on their definition. This is the check on those
 * arguments, because ep_dbl/ep_add and ep2_dbl/ep2_add are what the subgroup
 * tests run on points an attacker chose, and an argument is not a check.
 *
 * THE ORACLE IS AFFINE, with a real field inversion and explicit special
 * cases. That is deliberately not a second projective formula: it shares no
 * algebra with the code under test, it is slow and obviously correct, and it
 * cannot be wrong in the same direction.
 *
 * Inputs: random on-curve points NOT restricted to the subgroup, each in a
 * random projective representative, and every degenerate pair -- P+P, P+(-P),
 * P+O, O+P, O+O, dbl(O) -- plus the aliasing cases.
 *
 * Both groups, not just the twist. G1 moved to these formulas in issue #51 and
 * ep_read_compressed hands ep_in_subgroup whatever bytes arrived, so G1 needs
 * the coverage for the same reason G2 does.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elips/pairing.h"
#include "elips/ec.h"

static int fails, checks;
static void ok(int c, const char *what)
{ checks++; if (!c) { fails++; printf("  [FAIL] %s\n", what); } }

/* ---- random field elements ---- */
static unsigned long long st = 0x853c49e6748fea9bull;
static unsigned long long xs(void)
{ st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }

static void rand_fp(fp_t r)
{
    limb_t k[FP_LIMBS];
    for (int i = 0; i < FP_LIMBS; i++) k[i] = xs();
    k[FP_LIMBS - 1] >>= 8;              /* comfortably below p */
    fp_from_limbs(r, k);
}
static void rand_fp2(fp2_t r) { rand_fp(r[0]); rand_fp(r[1]); }

/* G1. On BN the cofactor is 1 -- #E(Fp) = r -- so every curve point is in G1
 * and there is no off-subgroup case to reach. On BLS12 the cofactor is
 * (x-1)^2/3, around 2^126 on BLS12-381, so essentially no random point is. */
#define EC_PT ep
#define EC_F  fp
#define EC_FT fp_t
#define EC_RANDF rand_fp
#ifdef ELIPS_FAMILY_BN
#define EC_EXPECT_OFFSUB 0
#else
#define EC_EXPECT_OFFSUB 1
#endif
#include "ec_group_tmpl.h"
#undef EC_PT
#undef EC_F
#undef EC_FT
#undef EC_RANDF
#undef EC_EXPECT_OFFSUB

/* G2. The cofactor is enormous on both families. */
#define EC_PT ep2
#define EC_F  fp2
#define EC_FT fp2_t
#define EC_RANDF rand_fp2
#define EC_EXPECT_OFFSUB 1
#include "ec_group_tmpl.h"
#undef EC_PT
#undef EC_F
#undef EC_FT
#undef EC_RANDF
#undef EC_EXPECT_OFFSUB

int main(int argc, char **argv)
{
    long want = (argc > 1) ? atol(argv[1]) : 20000;
    int rc1, rc2;

    printf("projective group law against affine [%s]\n", ELIPS_CURVE_NAME);

    rc1 = ep_grp_run(want, "G1");
    rc2 = ep2_grp_run(want, "G2");

    if (rc1 == 2 || rc2 == 2) return 2;
    return (rc1 || rc2) ? 1 : 0;
}
