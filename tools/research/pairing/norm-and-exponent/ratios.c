/* Primitive cost ratios on this machine: decides whether Toom-style
 * interpolation (many Fp adds, few Fp muls) can beat the current tower. */
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "elips/fp.h"
#include "elips/fpx.h"

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#define REPS 2000000

static uint64_t rng_state = 0x9e3779b97f4a7c15ull;
static uint64_t nextr(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}
static void rand_fp(fp_t a) {
    limb_t plain[FP_LIMBS];
    for (int i = 0; i < FP_LIMBS; i++) plain[i] = (limb_t)nextr();
    fp_from_limbs(a, plain);
}
static void rand_fp2(fp2_t a) { rand_fp(a[0]); rand_fp(a[1]); }

int main(void) {
    fp_t a, b, c;
    fp2_t a2, b2, c2;
    rand_fp(a); rand_fp(b); rand_fp(c);
    rand_fp2(a2); rand_fp2(b2); rand_fp2(c2);

    double t0, t1;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp_add(c, a, b); fp_add(a, c, b); }
    t1 = now();
    double add_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp_mul(c, a, b); fp_mul(a, c, b); }
    t1 = now();
    double mul_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp_sqr(c, a); fp_sqr(a, c); }
    t1 = now();
    double sqr_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp2_add(c2, a2, b2); fp2_add(a2, c2, b2); }
    t1 = now();
    double a2_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp2_mul(c2, a2, b2); fp2_mul(a2, c2, b2); }
    t1 = now();
    double m2_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    t0 = now();
    for (long i = 0; i < REPS; i++) { fp2_sqr(c2, a2); fp2_sqr(a2, c2); }
    t1 = now();
    double s2_ns = (t1 - t0) / (2.0 * REPS) * 1e9;

    printf("fp_add : %7.2f ns\n", add_ns);
    printf("fp_mul : %7.2f ns\n", mul_ns);
    printf("fp_sqr : %7.2f ns\n", sqr_ns);
    printf("fp2_add: %7.2f ns\n", a2_ns);
    printf("fp2_mul: %7.2f ns\n", m2_ns);
    printf("fp2_sqr: %7.2f ns\n", s2_ns);
    printf("M/A ratio: %.1f   M2 measured: %.2f Mp   S2 measured: %.2f Mp\n",
           mul_ns / add_ns, m2_ns / mul_ns, s2_ns / mul_ns);
    return 0;
}
