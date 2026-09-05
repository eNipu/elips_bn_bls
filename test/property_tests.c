/*
 * Property tests (issue #14).
 *
 * Fixed vectors pin down values; these pin down the algebraic laws that must
 * survive any refactor. Several of these checks did not exist anywhere in the
 * repository before, most importantly the subgroup checks: nothing verified
 * that a generated point actually has order r.
 *
 * Every failure sets the exit status. The pre-existing drivers printed
 * "failed" and returned 0, so CI could not see them.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_pairings.h>
#include <ELiPS_bn_bls/bn_generate_points.h>
#include <ELiPS_bn_bls/bn_twist.h>
#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/bls12_scm.h>
#include <ELiPS_bn_bls/bls12_G3_exp.h>
#include <ELiPS_bn_bls/curve_settings.h>

static int failures;

/* EFp12 has no comparison function in the library; points are equal when both
 * are at infinity or both coordinates match. */
static int EFp12_cmp_for_test(EFp12 *A, EFp12 *B)
{
    if (A->infinity || B->infinity) return !(A->infinity && B->infinity);
    return (Fp12_cmp(&A->x, &B->x) != 0) || (Fp12_cmp(&A->y, &B->y) != 0);
}

static void ok(int cond, const char *what)
{
    printf("  [%s] %s\n", cond ? "PASS" : "FAIL", what);
    if (!cond) failures++;
}

/* Map a G1 point between its EFp12 embedding and plain affine EFp. */
static void g1_to_efp(EFp *out, EFp12 *P, int bls)
{
    if (bls) { bls12_EFp12_to_EFp(out, P); return; }
    Fp_set(&out->x, &P->x.x0.x0.x0);
    Fp_set(&out->y, &P->y.x0.x0.x0);
    out->infinity = P->infinity;
}
static void efp_to_g1(EFp12 *out, EFp *P, int bls)
{
    if (bls) { bls12_EFp_to_EFp12(out, P); return; }
    Fp12_set_ui(&out->x, 0); Fp12_set_ui(&out->y, 0);
    Fp_set(&out->x.x0.x0.x0, &P->x);
    Fp_set(&out->y.x0.x0.x0, &P->y);
    out->infinity = P->infinity;
}
static void g2_to_efp2(EFp2 *out, EFp12 *Q, int bls)
{ if (bls) bls12_EFp12_to_EFp2(out, Q); else EFp12_to_EFp2(out, Q); }

int main(int argc, char **argv)
{
    int bls = (argc > 1 && strcmp(argv[1], "bls12") == 0);
    const char *name = bls ? "BLS12-461" : "BN-462";
    printf("property tests: %s\n", name);

    if (bls) bls12_inits(); else init_bn();

    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, 20260906UL);      /* fixed seed: reproducible failures */

    EFp12 P, Q;
    EFp12_init(&P); EFp12_init(&Q);
    if (bls) { bls12_generate_G1_point(&P); bls12_generate_G2_point(&Q); }
    else     { bn12_generate_G1_point(&P);  bn12_generate_G2_point(&Q);  }

    /* ---- 1. the generated G1 point satisfies the curve equation ---- */
    {
        EFp p1; EFp_init(&p1); g1_to_efp(&p1, &P, bls);
        Fp lhs, rhs; Fp_init(&lhs); Fp_init(&rhs);
        Fp_mul(&lhs, &p1.y, &p1.y);                 /* y^2       */
        Fp_mul(&rhs, &p1.x, &p1.x);
        Fp_mul(&rhs, &rhs, &p1.x);                  /* x^3       */
        if (bls) Fp_add_mpz(&rhs, &rhs, curve_parameters.curve_b);   /* +b */
        else     Fp_sub_mpz(&rhs, &rhs, curve_parameters.curve_b);   /* -b */
        ok(Fp_cmp(&lhs, &rhs) == 0, "G1 point satisfies the curve equation");
        Fp_clear(&lhs); Fp_clear(&rhs); EFp_clear(&p1);
    }

    /* ---- 2. subgroup order: [r]P == O.  Absent from the repo until now. ---- */
    {
        EFp p1, rp; EFp_init(&p1); EFp_init(&rp);
        g1_to_efp(&p1, &P, bls);
        EFp_SCM(&rp, &p1, curve_parameters.order);
        ok(rp.infinity == 1, "[r]P == O for the G1 generator");
        EFp_clear(&p1); EFp_clear(&rp);
    }
    {
        EFp2 q2, rq; EFp2_init(&q2); EFp2_init(&rq);
        g2_to_efp2(&q2, &Q, bls);
        EFp2_SCM(&rq, &q2, curve_parameters.order);
        ok(rq.infinity == 1, "[r]Q == O for the G2 generator");
        EFp2_clear(&q2); EFp2_clear(&rq);
    }

    /* ---- 3. the pairing is non-degenerate and lands in mu_r ---- */
    Fp12 z; Fp12_init(&z);
    if (bls) bls12_opt_ate(&z, &P, &Q); else bn12_opt_ate(&z, &P, &Q);
    ok(Fp12_cmp_one(&z) != 0, "e(P,Q) != 1  (non-degenerate)");
    ok(Fp12_cmp_zero(&z) != 0, "e(P,Q) != 0");
    {
        Fp12 zr; Fp12_init(&zr);
        Fp12_pow(&zr, &z, curve_parameters.order);
        ok(Fp12_cmp_one(&zr) == 0, "e(P,Q)^r == 1  (lies in mu_r)");
        Fp12_clear(&zr);
    }

    /* ---- 4. bilinearity in both arguments simultaneously ---- */
    {
        mpz_t a, b, ab;
        mpz_inits(a, b, ab, NULL);
        mpz_urandomm(a, st, curve_parameters.order);
        mpz_urandomm(b, st, curve_parameters.order);
        mpz_mul(ab, a, b); mpz_mod(ab, ab, curve_parameters.order);

        EFp p1, ap; EFp_init(&p1); EFp_init(&ap);
        g1_to_efp(&p1, &P, bls); EFp_SCM(&ap, &p1, a);
        EFp12 aP; EFp12_init(&aP); efp_to_g1(&aP, &ap, bls);

        EFp2 q2, bq; EFp2_init(&q2); EFp2_init(&bq);
        g2_to_efp2(&q2, &Q, bls); EFp2_SCM(&bq, &q2, b);
        EFp12 bQ; EFp12_init(&bQ);
        if (bls) bls12_EFp2_to_EFp12(&bQ, &bq); else EFp2_to_EFp12(&bQ, &bq);

        Fp12 lhs, rhs; Fp12_init(&lhs); Fp12_init(&rhs);
        if (bls) bls12_opt_ate(&lhs, &aP, &bQ); else bn12_opt_ate(&lhs, &aP, &bQ);
        Fp12_pow(&rhs, &z, ab);
        ok(Fp12_cmp(&lhs, &rhs) == 0, "e([a]P,[b]Q) == e(P,Q)^(ab)");

        Fp12_clear(&lhs); Fp12_clear(&rhs);
        EFp_clear(&p1); EFp_clear(&ap); EFp12_clear(&aP);
        EFp2_clear(&q2); EFp2_clear(&bq); EFp12_clear(&bQ);
        mpz_clears(a, b, ab, NULL);
    }

    /* ---- 5. additivity in the first argument.
     * Stronger than the scalar bilinearity the old suite tested, because it
     * exercises the group law rather than just repeated scaling. ---- */
    {
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        mpz_urandomm(a, st, curve_parameters.order);
        mpz_urandomm(b, st, curve_parameters.order);

        EFp p1, p_a, p_b, p_sum;
        EFp_init(&p1); EFp_init(&p_a); EFp_init(&p_b); EFp_init(&p_sum);
        g1_to_efp(&p1, &P, bls);
        EFp_SCM(&p_a, &p1, a);
        EFp_SCM(&p_b, &p1, b);
        EFp_ECA(&p_sum, &p_a, &p_b);

        EFp12 A, B, S; EFp12_init(&A); EFp12_init(&B); EFp12_init(&S);
        efp_to_g1(&A, &p_a, bls); efp_to_g1(&B, &p_b, bls); efp_to_g1(&S, &p_sum, bls);

        Fp12 ea, eb, es, prod;
        Fp12_init(&ea); Fp12_init(&eb); Fp12_init(&es); Fp12_init(&prod);
        if (bls) { bls12_opt_ate(&ea,&A,&Q); bls12_opt_ate(&eb,&B,&Q); bls12_opt_ate(&es,&S,&Q); }
        else     { bn12_opt_ate(&ea,&A,&Q);  bn12_opt_ate(&eb,&B,&Q);  bn12_opt_ate(&es,&S,&Q);  }
        Fp12_mul(&prod, &ea, &eb);
        ok(Fp12_cmp(&es, &prod) == 0, "e(P1+P2,Q) == e(P1,Q)*e(P2,Q)");

        Fp12_clear(&ea); Fp12_clear(&eb); Fp12_clear(&es); Fp12_clear(&prod);
        EFp12_clear(&A); EFp12_clear(&B); EFp12_clear(&S);
        EFp_clear(&p1); EFp_clear(&p_a); EFp_clear(&p_b); EFp_clear(&p_sum);
        mpz_clears(a, b, NULL);
    }

    /* ---- 6. split-scalar routines agree with the plain ones.
     * Guards the M8 fix: those routines right-align each sub-scalar's binary
     * expansion into a fixed-width row, and the memmove that did it copied
     * sizeof(buffer) rather than the digit count, running off the end of the
     * row for any sub-scalar shorter than the widest one -- about two thirds
     * of random scalars. BLS12 only; BN has no split variants. ---- */
    if (bls) {
        int trials = 6, bad = 0;
        for (int t = 0; t < trials; t++) {
            mpz_t k; mpz_init(k);
            mpz_urandomm(k, st, curve_parameters.order);

            EFp12 r_plain, r_2, r_4;
            EFp12_init(&r_plain); EFp12_init(&r_2); EFp12_init(&r_4);

            bls12_plain_G1_scm(&r_plain, &P, k);
            bls12_2split_G1_scm(&r_2, &P, k);
            if (EFp12_cmp_for_test(&r_plain, &r_2)) bad++;

            bls12_plain_G2_scm(&r_plain, &Q, k);
            bls12_2split_G2_scm(&r_2, &Q, k);
            bls12_4split_G2_scm(&r_4, &Q, k);
            if (EFp12_cmp_for_test(&r_plain, &r_2)) bad++;
            if (EFp12_cmp_for_test(&r_plain, &r_4)) bad++;

            Fp12 e_plain, e_2, e_4;
            Fp12_init(&e_plain); Fp12_init(&e_2); Fp12_init(&e_4);
            bls12_plain_G3_exp(&e_plain, &z, k);
            bls12_2split_G3_exp(&e_2, &z, k);
            bls12_4split_G3_exp(&e_4, &z, k);
            if (Fp12_cmp(&e_plain, &e_2) != 0) bad++;
            if (Fp12_cmp(&e_plain, &e_4) != 0) bad++;

            Fp12_clear(&e_plain); Fp12_clear(&e_2); Fp12_clear(&e_4);
            EFp12_clear(&r_plain); EFp12_clear(&r_2); EFp12_clear(&r_4);
            mpz_clear(k);
        }
        ok(bad == 0, "split-scalar G1/G2 SCM and G3 exp agree with the plain versions");
    }

    Fp12_clear(&z); EFp12_clear(&P); EFp12_clear(&Q);
    gmp_randclear(st);

    printf("%s: %s (%d failure%s)\n", name, failures ? "FAILED" : "OK",
           failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
