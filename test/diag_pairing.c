/*
 * Isolate which pairing stage disagrees with the reference.
 *
 * The library ships two final exponentiations that must be mathematically
 * identical: the "plain" one raises to (p^12-1)/r by the direct route
 * (easy part, then (p^4-p^2+1)/r), and the "optimal" one uses an addition
 * chain. Since (p^6-1)(p^2+1)(p^4-p^2+1)/r == (p^12-1)/r, both must agree.
 * Comparing them needs no external oracle at all.
 *
 * Both final exponentiations write through their input pointer (defect A3),
 * so every call gets its own copy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_pairings.h>
#include <ELiPS_bn_bls/bn_generate_points.h>
#include <ELiPS_bn_bls/bn_miller_optate.h>
#include <ELiPS_bn_bls/bn_final_exp.h>
#include <ELiPS_bn_bls/bn_twist.h>
#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_miller_optate.h>
#include <ELiPS_bn_bls/bls12_finalexp.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/curve_settings.h>

static void flat(Fp12 *f, mpz_t *c)
{
    Fp6 *six[2] = { &f->x0, &f->x1 };
    for (int s = 0; s < 2; s++) {
        Fp2 *cc[3] = { &six[s]->x0, &six[s]->x1, &six[s]->x2 };
        for (int j = 0; j < 3; j++) {
            mpz_mod(c[s*6+j*2  ], cc[j]->x0.x0, curve_parameters.prime);
            mpz_mod(c[s*6+j*2+1], cc[j]->x1.x0, curve_parameters.prime);
        }
    }
}
static void jarr(FILE *o, const char *k, mpz_t *c, int comma)
{
    fprintf(o, "  \"%s\": [", k);
    for (int i = 0; i < 12; i++) {
        char *s = mpz_get_str(NULL, 16, c[i]);
        fprintf(o, "%s\"%s\"", i ? ", " : "", s); free(s);
    }
    fprintf(o, "]%s\n", comma ? "," : "");
}
static int same(mpz_t *a, mpz_t *b)
{ for (int i = 0; i < 12; i++) if (mpz_cmp(a[i], b[i])) return 0; return 1; }

int main(int argc, char **argv)
{
    int bls = (argc > 1 && strcmp(argv[1], "bls12") == 0);
    if (bls) bls12_inits(); else init_bn();

    EFp12 P, Q; EFp12_init(&P); EFp12_init(&Q);
    if (bls) { bls12_generate_G1_point(&P); bls12_generate_G2_point(&Q); }
    else     { bn12_generate_G1_point(&P);  bn12_generate_G2_point(&Q);  }

    EFp mp; EFp2 mq; EFp_init(&mp); EFp2_init(&mq);
    if (bls) { bls12_EFp12_to_EFp(&mp, &P); bls12_EFp12_to_EFp2(&mq, &Q); }
    else     { Fp_set(&mp.x, &P.x.x0.x0.x0); Fp_set(&mp.y, &P.y.x0.x0.x0);
               EFp12_to_EFp2(&mq, &Q); }

    Fp12 raw, a, b; Fp12_init(&raw); Fp12_init(&a); Fp12_init(&b);
    if (bls) bls12_Miller_algo_for_opt_ate(&raw, &P, &Q);
    else     Miller_algo_for_opt_ate(&raw, &P, &Q);

    mpz_t craw[12], cpl[12], copt[12];
    for (int i = 0; i < 12; i++) { mpz_init(craw[i]); mpz_init(cpl[i]); mpz_init(copt[i]); }
    flat(&raw, craw);

    Fp12_set(&a, &raw);
    if (bls) bls12_finalexp_plain(&a, &a); else bn_final_exp_plain(&a, &a);
    flat(&a, cpl);

    Fp12_set(&b, &raw);
    if (bls) bls12_finalexp_optimal(&b, &b); else bn_final_exp_optimal(&b, &b);
    flat(&b, copt);

    printf("curve: %s\n", bls ? "BLS12-461" : "BN-462");
    printf("final_exp_plain == final_exp_optimal : %s\n",
           same(cpl, copt) ? "YES" : "NO  <-- the two must be mathematically identical");

    FILE *o = fopen(argc > 2 ? argv[2] : "diag.json", "w");
    fprintf(o, "{\n  \"curve\": \"%s\",\n", bls ? "BLS12-461" : "BN-462");
    { char *s = mpz_get_str(NULL,16,curve_parameters.prime);
      fprintf(o, "  \"p\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,curve_parameters.order);
      fprintf(o, "  \"r\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mp.x.x0); fprintf(o,"  \"px\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mp.y.x0); fprintf(o,"  \"py\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mq.x.x0.x0); fprintf(o,"  \"qx0\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mq.x.x1.x0); fprintf(o,"  \"qx1\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mq.y.x0.x0); fprintf(o,"  \"qy0\": \"%s\",\n", s); free(s); }
    { char *s = mpz_get_str(NULL,16,mq.y.x1.x0); fprintf(o,"  \"qy1\": \"%s\",\n", s); free(s); }
    jarr(o, "miller_raw", craw, 1);
    jarr(o, "fe_plain",   cpl,  1);
    jarr(o, "fe_optimal", copt, 0);
    fprintf(o, "}\n");
    fclose(o);
    return 0;
}
