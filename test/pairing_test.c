/* Does the new pairing equal the value Phase 0 established as correct?
 *
 * The legacy finalexp_plain path was proved in Phase 0 to equal
 * f^((p^12-1)/r) exactly (issue #16: it is the *optimal* variant that is off by
 * a factor). So legacy_miller followed by legacy_finalexp_plain is the
 * reference here, and it shares no code with the new layer. */
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "elips/pairing.h"
#include <ELiPS_bn_bls/curve_settings.h>
#ifdef ELIPS_FAMILY_BN
#  include <ELiPS_bn_bls/bn_inits.h>
#  include <ELiPS_bn_bls/bn_generate_points.h>
#  include <ELiPS_bn_bls/bn_twist.h>
#  include <ELiPS_bn_bls/bn_miller_optate.h>
#  include <ELiPS_bn_bls/bn_final_exp.h>
#  define CURVE_INIT()        init_bn()
#  define GEN_G1(P)           bn12_generate_G1_point(P)
#  define GEN_G2(Q)           bn12_generate_G2_point(Q)
#  define TO_EFP(o, P)        do { Fp_set(&(o).x,&(P).x.x0.x0.x0); \
                                   Fp_set(&(o).y,&(P).y.x0.x0.x0); } while (0)
#  define TO_EFP2(o, Q)       EFp12_to_EFp2(&(o), &(Q))
#  define REF_MILLER(f,P,Q)   Miller_algo_for_opt_ate(f,P,Q)
#  define REF_FINALEXP(f)     bn_final_exp_plain(f,f)
#else
#  include <ELiPS_bn_bls/bls12_inits.h>
#  include <ELiPS_bn_bls/bls12_generate_points.h>
#  include <ELiPS_bn_bls/bls12_twist.h>
#  include <ELiPS_bn_bls/bls12_miller_optate.h>
#  include <ELiPS_bn_bls/bls12_finalexp.h>
#  define CURVE_INIT()        bls12_inits()
#  define GEN_G1(P)           bls12_generate_G1_point(P)
#  define GEN_G2(Q)           bls12_generate_G2_point(Q)
#  define TO_EFP(o, P)        bls12_EFp12_to_EFp(&(o), &(P))
#  define TO_EFP2(o, Q)       bls12_EFp12_to_EFp2(&(o), &(Q))
#  define REF_MILLER(f,P,Q)   bls12_Miller_algo_for_opt_ate(f,P,Q)
#  define REF_FINALEXP(f)     bls12_finalexp_plain(f,f)
#endif

static void ld(fp_t r, const mpz_t v)
{ limb_t l[FP_LIMBS]; memset(l,0,sizeof l);
  mpz_export(l,NULL,-1,sizeof(limb_t),0,0,v); fp_from_limbs(r,l); }
static void st(mpz_t o, const fp_t a)
{ limb_t l[FP_LIMBS]; fp_to_limbs(l,a); mpz_import(o,FP_LIMBS,-1,sizeof(limb_t),0,0,l); }

int main(void)
{
    CURVE_INIT();
    printf("pairing test [%s]\n", ELIPS_CURVE_NAME);
    EFp12 P12, Q12; EFp12_init(&P12); EFp12_init(&Q12);
    GEN_G1(&P12);
    GEN_G2(&Q12);

    EFp mp; EFp2 mq; EFp_init(&mp); EFp2_init(&mq);
    TO_EFP(mp, P12);
    TO_EFP2(mq, Q12);

    /* new layer */
    fp_t px, py; fp2_t qx, qy;
    ld(px, mp.x.x0); ld(py, mp.y.x0);
    ld(qx[0], mq.x.x0.x0); ld(qx[1], mq.x.x1.x0);
    ld(qy[0], mq.y.x0.x0); ld(qy[1], mq.y.x1.x0);
    fp12_t nf, nz;
    pairing_miller(nf, qx, qy, px, py);
    pairing_final_exp_plain(nz, nf);

    /* legacy reference (Miller + the plain, provably correct final exp) */
    Fp12 lf; Fp12_init(&lf);
    REF_MILLER(&lf, &P12, &Q12);
    REF_FINALEXP(&lf);

    mpz_t a, b; mpz_inits(a, b, NULL);
    int diff = 0;
    const limb_t *nc[12] = { nz[0][0][0], nz[0][0][1], nz[0][1][0], nz[0][1][1],
                           nz[0][2][0], nz[0][2][1], nz[1][0][0], nz[1][0][1],
                           nz[1][1][0], nz[1][1][1], nz[1][2][0], nz[1][2][1] };
    mpz_t *lc[12];
    Fp6 *six[2] = { &lf.x0, &lf.x1 };
    int k = 0;
    for (int s2 = 0; s2 < 2; s2++) {
        Fp2 *cc[3] = { &six[s2]->x0, &six[s2]->x1, &six[s2]->x2 };
        for (int j = 0; j < 3; j++) { lc[k++] = &cc[j]->x0.x0; lc[k++] = &cc[j]->x1.x0; }
    }
    for (int i = 0; i < 12; i++) {
        st(a, nc[i]);
        mpz_mod(b, *lc[i], curve_parameters.prime);
        if (mpz_cmp(a, b) != 0) {
            diff++;
            if (diff <= 3) gmp_printf("  coord %2d  new=%Zx\n            ref=%Zx\n", i, a, b);
        }
    }
    printf("  [%s] new pairing == Phase 0 reference value\n", diff ? "FAIL" : "PASS");

    /* Independent properties, so a shared misunderstanding of the conventions
     * cannot make both sides agree on a wrong answer. */
    int fails = diff ? 1 : 0;
    fp12_t one, chk;
    fp12_set_one(one);
    if (fp12_eq(nz, one)) { printf("  [FAIL] pairing is degenerate\n"); fails++; }
    else printf("  [PASS] e(P,Q) != 1\n");

    {   /* e(P,Q)^r == 1 */
        limb_t rr[FP_LIMBS]; memset(rr, 0, sizeof rr);
        mpz_export(rr, NULL, -1, sizeof(limb_t), 0, 0, curve_parameters.order);
        int rbits = (int)mpz_sizeinbase(curve_parameters.order, 2);
        fp12_exp(chk, nz, rr, rbits);
        if (!fp12_eq(chk, one)) { printf("  [FAIL] e(P,Q)^r != 1\n"); fails++; }
        else printf("  [PASS] e(P,Q)^r == 1\n");
    }

#ifdef ELIPS_FAMILY_BLS12
    {   /* The fast chain must equal the exact value cubed. It computes e^3 by
         * design, because 3*lambda = (x-1)^2 (x+p) (x^2+p^2-1) + 3 has a short
         * evaluation and lambda alone does not. See issue #16. */
        fp12_t fast, cube;
        pairing_final_exp_fast(fast, nf);
        fp12_sqr(cube, nz); fp12_mul(cube, cube, nz);
        if (!fp12_eq(fast, cube)) { printf("  [FAIL] fast chain != exact^3\n"); fails++; }
        else printf("  [PASS] fast final exponentiation == exact^3\n");
    }
#else
    {   /* BN's hard part decomposes exactly, so its fast chain must equal the
         * exact value with no factor -- unlike BLS12. */
        fp12_t fast;
        pairing_final_exp_fast(fast, nf);
        if (!fp12_eq(fast, nz)) { printf("  [FAIL] BN fast chain != exact e\n"); fails++; }
        else printf("  [PASS] fast final exponentiation == exact e\n");
    }
#endif

    printf("%s\n", fails ? "FAILED" : "OK");
    mpz_clears(a, b, NULL);
    return fails ? 1 : 0;
}
