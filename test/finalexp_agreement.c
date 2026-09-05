/*
 * The library ships two final exponentiations per curve family. They raise to
 * exponents that are algebraically identical:
 *
 *     (p^6-1)(p^2+1)(p^4-p^2+1)/r  ==  (p^12-1)/r
 *
 * so they must produce the same value. They do not. See issue #16.
 *
 * Registered in CTest with WILL_FAIL, so the defect stays visible without
 * blocking the build. Drop WILL_FAIL when #16 is fixed and this becomes an
 * ordinary regression test.
 *
 * Both routines write through their input pointer (defect A3), so each call
 * gets a private copy.
 */
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_generate_points.h>
#include <ELiPS_bn_bls/bn_miller_optate.h>
#include <ELiPS_bn_bls/bn_final_exp.h>
#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_miller_optate.h>
#include <ELiPS_bn_bls/bls12_finalexp.h>
#include <ELiPS_bn_bls/curve_settings.h>

int main(int argc, char **argv)
{
    int bls = (argc > 1 && strcmp(argv[1], "bls12") == 0);
    if (bls) bls12_inits(); else init_bn();

    EFp12 P, Q; EFp12_init(&P); EFp12_init(&Q);
    if (bls) { bls12_generate_G1_point(&P); bls12_generate_G2_point(&Q); }
    else     { bn12_generate_G1_point(&P);  bn12_generate_G2_point(&Q);  }

    Fp12 raw, a, b;
    Fp12_init(&raw); Fp12_init(&a); Fp12_init(&b);
    if (bls) bls12_Miller_algo_for_opt_ate(&raw, &P, &Q);
    else     Miller_algo_for_opt_ate(&raw, &P, &Q);

    Fp12_set(&a, &raw);
    if (bls) bls12_finalexp_plain(&a, &a);   else bn_final_exp_plain(&a, &a);
    Fp12_set(&b, &raw);
    if (bls) bls12_finalexp_optimal(&b, &b); else bn_final_exp_optimal(&b, &b);

    int agree = (Fp12_cmp(&a, &b) == 0);
    printf("%s: final_exp_plain == final_exp_optimal : %s\n",
           bls ? "BLS12-461" : "BN-462", agree ? "YES" : "NO (issue #16)");
    return agree ? 0 : 1;
}
