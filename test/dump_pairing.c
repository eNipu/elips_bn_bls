/*
 * Export the library's own pairing inputs and output as JSON, so the reference
 * implementation can recompute e(P,Q) from the same points and compare.
 *
 * Export-then-verify sidesteps having to reproduce the library's random point
 * generation in Python, and it tests exactly the value the library produces.
 *
 * The 8-sparse mapping mutates P and Q in place, so the coordinates are dumped
 * BEFORE the pairing is called.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_clears.h>
#include <ELiPS_bn_bls/bn_pairings.h>
#include <ELiPS_bn_bls/bn_generate_points.h>
#include <ELiPS_bn_bls/bn_twist.h>
#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/curve_settings.h>

static void jz(FILE *f, const char *k, mpz_t v, int comma)
{
    char *s = mpz_get_str(NULL, 16, v);
    fprintf(f, "  \"%s\": \"%s\"%s\n", k, s, comma ? "," : "");
    free(s);
}

int main(int argc, char **argv)
{
    int bls = (argc > 1 && strcmp(argv[1], "bls12") == 0);
    const char *curve = bls ? "BLS12-461" : "BN-462";

    EFp12 P, Q;
    Fp12 f;
    EFp  mp;  EFp2 mq;

    if (bls) bls12_inits(); else init_bn();
    EFp12_init(&P); EFp12_init(&Q); Fp12_init(&f);
    EFp_init(&mp);  EFp2_init(&mq);

    if (bls) { bls12_generate_G1_point(&P); bls12_generate_G2_point(&Q); }
    else     { bn12_generate_G1_point(&P);  bn12_generate_G2_point(&Q);  }

    /* Pull out the affine G1 point and the twist representative of G2. */
    if (bls) { bls12_EFp12_to_EFp(&mp, &P); bls12_EFp12_to_EFp2(&mq, &Q); }
    else     { Fp_set(&mp.x, &P.x.x0.x0.x0); Fp_set(&mp.y, &P.y.x0.x0.x0);
               EFp12_to_EFp2(&mq, &Q); }

    FILE *out = fopen(argc > 2 ? argv[2] : "pairing_dump.json", "w");
    fprintf(out, "{\n  \"curve\": \"%s\",\n", curve);
    jz(out, "p",   curve_parameters.prime, 1);
    jz(out, "r",   curve_parameters.order, 1);
    jz(out, "px",  mp.x.x0, 1);
    jz(out, "py",  mp.y.x0, 1);
    jz(out, "qx0", mq.x.x0.x0, 1);
    jz(out, "qx1", mq.x.x1.x0, 1);
    jz(out, "qy0", mq.y.x0.x0, 1);
    jz(out, "qy1", mq.y.x1.x0, 1);

    if (bls) bls12_opt_ate(&f, &P, &Q); else bn12_opt_ate(&f, &P, &Q);

    mpz_t c[12];
    for (int i = 0; i < 12; i++) mpz_init(c[i]);
    Fp6 *six[2] = { &f.x0, &f.x1 };
    for (int s = 0; s < 2; s++) {
        Fp2 *cc[3] = { &six[s]->x0, &six[s]->x1, &six[s]->x2 };
        for (int j = 0; j < 3; j++) {
            mpz_set(c[s*6 + j*2    ], cc[j]->x0.x0);
            mpz_set(c[s*6 + j*2 + 1], cc[j]->x1.x0);
        }
    }
    fprintf(out, "  \"f\": [");
    for (int i = 0; i < 12; i++) {
        char *s = mpz_get_str(NULL, 16, c[i]);
        fprintf(out, "%s\"%s\"", i ? ", " : "", s);
        free(s);
    }
    fprintf(out, "]\n}\n");
    fclose(out);
    printf("dumped %s\n", curve);
    return 0;
}
