/*
 * Example 1 -- the pairing itself.
 *
 * What a pairing gives you: a map e(P, Q) from a point of G1 and a point of G2
 * into G_T, which is bilinear. Scalars slide across it:
 *
 *     e([a]P, [b]Q) = e(P, Q)^(a*b)
 *
 * That one identity is what every pairing-based protocol is built on. This
 * example demonstrates it, shows the input validation, and times the pieces.
 *
 * Build:  cmake -B build && cmake --build build
 * Run:    ./build/examples/elips_example_pairing
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <time.h>

#include "elips/pairing.h"
#include "elips/random.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

int main(void)
{
    printf("== ELiPS pairing example, curve %s ==\n\n", ELIPS_CURVE_NAME);

    /* The generators. On BLS12-381 these are the ones the specification
     * publishes, so they match any other implementation. */
    ep_t  P;
    ep2_t Q;
    ep_generator(&P);
    ep2_generator(&Q);

    /* Two secret scalars from the operating system's CSPRNG. Every routine
     * here returns a status; a silently non-random scalar is a private key an
     * attacker can guess, so the check is not optional. */
    limb_t a[ELIPS_ORDER_LIMBS], b[ELIPS_ORDER_LIMBS];
    if (elips_random_scalar(a) != 0 || elips_random_scalar(b) != 0) {
        fprintf(stderr, "no system entropy available\n");
        return 1;
    }

    ep_t  aP;
    ep2_t bQ;
    ep_mul(&aP, &P, a, ELIPS_ORDER_BITS);
    ep2_mul(&bQ, &Q, b, ELIPS_ORDER_BITS);

    /* elips_pairing validates both inputs: it returns 0, and leaves the result
     * at one, if either point is the identity or outside its order-r subgroup.
     * Ignoring the return value gets you a useless answer, not a subtly wrong
     * one -- but check it anyway. */
    fp12_t base, lhs, rhs;
    if (!elips_pairing(base, &P, &Q))   { fprintf(stderr, "generators rejected\n"); return 1; }
    if (!elips_pairing(lhs, &aP, &bQ))  { fprintf(stderr, "multiples rejected\n");  return 1; }

    /* e(P,Q)^(a*b), computed as two exponentiations so the example does not
     * have to implement multiplication modulo r. */
    fp12_exp(rhs, base, a, ELIPS_ORDER_BITS);
    fp12_exp(rhs, rhs, b, ELIPS_ORDER_BITS);

    printf("bilinearity   e([a]P, [b]Q) == e(P, Q)^(ab)   %s\n",
           fp12_eq(lhs, rhs) ? "holds" : "FAILED");

    fp12_t one, chk;
    fp12_set_one(one);
    printf("non-degenerate   e(P, Q) != 1                 %s\n",
           !fp12_eq(base, one) ? "holds" : "FAILED");

    fp12_exp(chk, base, ELIPS_ORDER, ELIPS_ORDER_BITS);
    printf("in mu_r          e(P, Q)^r == 1               %s\n",
           fp12_eq(chk, one) ? "holds" : "FAILED");

    /* The validation is the interesting part for a caller: bad input is
     * refused rather than processed. */
    ep_t O;
    ep_set_infinity(&O);
    printf("\nvalidation:\n");
    printf("  the identity in G1 is           %s\n",
           elips_pairing(base, &O, &Q) ? "ACCEPTED (bug)" : "refused");

    /* Timing. Nothing here is a benchmark of the machine; it is a sense of
     * where the cost sits. */
    const int n = 50;
    double t0 = now();
    for (int i = 0; i < n; i++) elips_pairing(base, &P, &Q);
    double whole = (now() - t0) / n * 1e6;

    fp_t px, py;
    fp2_t qx, qy;
    ep_to_affine(px, py, &P);
    ep2_to_affine(qx, qy, &Q);
    fp12_t f;
    t0 = now();
    for (int i = 0; i < n; i++) pairing_miller(f, qx, qy, px, py);
    double miller = (now() - t0) / n * 1e6;

    t0 = now();
    for (int i = 0; i < n; i++) pairing_final_exp_fast(base, f);
    double fexp = (now() - t0) / n * 1e6;

    printf("\ntiming, microseconds, %d iterations:\n", n);
    printf("  Miller loop                  %8.1f\n", miller);
    printf("  final exponentiation (fast)  %8.1f\n", fexp);
    printf("  elips_pairing (both, plus subgroup checks on the inputs) %8.1f\n", whole);
    printf("\n  The gap is the two subgroup checks, each a full scalar\n"
           "  multiplication by r. They are what makes an attacker-supplied\n"
           "  point safe to pair.\n");

#ifdef ELIPS_FAMILY_BLS12
    printf("\nNote: on BLS12 the value returned is e^3, not e. That is a\n"
           "property of the standard final-exponentiation chain, which RELIC\n"
           "and the original library share. Since gcd(3, r) = 1 it is still\n"
           "bilinear and non-degenerate, so protocols that compare pairings\n"
           "are unaffected. Use pairing_final_exp_plain if you need e exactly.\n");
#endif
    return 0;
}
