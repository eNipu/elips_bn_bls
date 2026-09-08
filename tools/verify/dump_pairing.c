/*
 * Prints the generators and the EXACT pairing e(G1, G2), for
 * tools/verify/crosscheck_pyecc.py to compare against an implementation that
 * shares nothing with this one.
 *
 * The exact pairing, from pairing_final_exp_plain, not elips_pairing: the
 * latter returns e^3 on BLS12 by design, and comparing that against another
 * library would fail for a reason that is not a bug.
 *
 * The generators are printed too. A cross-check that agrees on the output
 * while disagreeing on the input has proved nothing.
 */
#include <stdio.h>
#include <gmp.h>

#include "elips/ec.h"
#include "elips/fpx.h"
#include "elips/pairing.h"

static void emit(const char *label, const fp_t a)
{
    limb_t l[FP_LIMBS];
    mpz_t v;
    mpz_init(v);
    fp_to_limbs(l, a);
    mpz_import(v, FP_LIMBS, -1, sizeof(limb_t), 0, 0, l);
    gmp_printf("%s %Zd\n", label, v);
    mpz_clear(v);
}

int main(void)
{
    ep_t P;
    ep2_t Q;
    fp2_t qx, qy;
    fp_t px, py;
    fp12_t f, e;

    ep_generator(&P);
    ep2_generator(&Q);
    if (!ep_to_affine(px, py, &P) || !ep2_to_affine(qx, qy, &Q)) {
        fprintf(stderr, "dump_pairing: a generator is the identity\n");
        return 1;
    }

    pairing_miller(f, qx, qy, px, py);
    pairing_final_exp_plain(e, f);

    emit("P.x", px);
    emit("P.y", py);
    emit("Q.x0", qx[0]); emit("Q.x1", qx[1]);
    emit("Q.y0", qy[0]); emit("Q.y1", qy[1]);

    /* a[0] = (g0, g2, g4) and a[1] = (g1, g3, g5), by powers of w. */
    for (int d = 0; d < 2; d++)
        for (int c = 0; c < 3; c++)
            for (int b = 0; b < 2; b++) {
                char lbl[32];
                snprintf(lbl, sizeof lbl, "e %d %d %d", d, c, b);
                emit(lbl, e[d][c][b]);
            }
    return 0;
}
