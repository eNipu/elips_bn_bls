/*
 * Repeated timing with a reported spread, so a regression can be told apart
 * from a noisy machine.
 *
 * The measurement this replaces was one sample per number, taken by hand. That
 * is not enough to make a claim with. During this work the same routine on the
 * same machine read 1.53x and then 1.35x across two runs, minutes apart, which
 * is larger than several of the speedups the repository claims. A single sample
 * cannot distinguish those two situations; that is the whole reason this file
 * exists.
 *
 * Measurements are INTERLEAVED, and that is the whole design.
 *
 * The obvious structure -- all reps of A, then all reps of B -- was tried
 * first and is wrong. Its reported spread understates the real error badly,
 * because reps taken back to back share a CPU frequency and a cache state
 * while separate runs do not. Measured on the machine this was written on,
 * eight of twelve operations disagreed between two consecutive runs by MORE
 * than the spread each run reported for itself, by up to 27.8%. A spread that
 * does not bound the run-to-run difference is worse than no spread, because it
 * invites confidence.
 *
 * So one rep now times every operation in turn before the next rep starts. A
 * frequency change or a noisy neighbour then lands on all of them together,
 * which leaves the RATIOS between operations trustworthy even when the
 * absolute numbers wander. Ratios are what the questions here actually need:
 * is GLV faster than the plain ladder, is this commit slower than the last.
 *
 * What it reports, per operation:
 *
 *   median   the headline number, over interleaved reps. Robust to a single
 *            descheduled sample in a way the mean is not.
 *   min      the cleanest sample seen. Closest to the cost of the code rather
 *            than of the machine around it.
 *   spread   (max - min) / median, as a percent. Read it as a health check on
 *            the machine, not as a confidence interval: above roughly 10% the
 *            absolute numbers are not worth quoting, though the ratios may
 *            still be.
 *
 * Not a correctness test. It calls routines the test suite already covers, and
 * only sinks results so the optimiser cannot delete the work.
 */
/* clock_gettime and CLOCK_MONOTONIC are POSIX, not C11, and the project builds
 * with -std=c11 which hides them behind this. Same reason test/dudect_test.c
 * does it. */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "elips/ec.h"
#include "elips/fpx.h"
#include "elips/pairing.h"
#include "elips/random.h"
#include "elips/hash_to_curve.h"

#ifndef BENCH_REPS
#define BENCH_REPS 9        /* odd, so the median is an observed sample */
#endif
#ifndef BENCH_INNER
#define BENCH_INNER 60      /* per timed run, to swamp clock granularity */
#endif
#define BENCH_MAX_REPS 999

static double now_us(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1e6 + (double)t.tv_nsec * 1e-3;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* --- what is being measured -------------------------------------------- */

static ep_t   g1, sink_ep;
static ep2_t  g2, sink_ep2;
static fp12_t sink_fp12, gt_base;
static limb_t scalars[BENCH_MAX_REPS][ELIPS_ORDER_LIMBS];
static int    scalar_idx;

static const limb_t *next_scalar(void)
{
    const limb_t *k = scalars[scalar_idx];
    scalar_idx = (scalar_idx + 1) % BENCH_MAX_REPS;
    return k;
}

static void b_ep_mul(void)      { ep_mul(&sink_ep, &g1, next_scalar(), ELIPS_ORDER_BITS); }
static void b_ep_mul_glv(void)  { ep_mul_glv(&sink_ep, &g1, next_scalar(), ELIPS_ORDER_BITS); }
static void b_ep2_mul(void)     { ep2_mul(&sink_ep2, &g2, next_scalar(), ELIPS_ORDER_BITS); }
static void b_ep2_mul_glv(void) { ep2_mul_glv(&sink_ep2, &g2, next_scalar(), ELIPS_ORDER_BITS); }
static void b_gt_exp(void)      { fp12_exp(sink_fp12, gt_base, next_scalar(), ELIPS_ORDER_BITS); }
static void b_gt_exp_ct(void)   { fp12_exp_gt(sink_fp12, gt_base, next_scalar(), ELIPS_ORDER_BITS); }
static void b_pairing(void)     { elips_pairing(sink_fp12, &g1, &g2); }

/* Miller and the final exponentiation separately, because they are optimised
 * independently and a change to one is invisible in the total. */
static fp2_t qx, qy;
static fp_t  px, py;
static void b_miller(void)      { pairing_miller(sink_fp12, qx, qy, px, py); }
static void b_final_exp(void)   { pairing_final_exp_fast(sink_fp12, gt_base); }
static void b_ep_in_sub(void)   { if (ep_in_subgroup(&g1)) sink_ep.z[0] ^= 1; }
static void b_ep2_in_sub(void)  { if (ep2_in_subgroup(&g2)) sink_ep2.z[0][0] ^= 1; }

static uint8_t h2c_msg[32];
static void b_hash_g1(void) { elips_hash_to_g1(&sink_ep,  h2c_msg, sizeof h2c_msg, (const uint8_t *)"BENCH", 5); }
static void b_hash_g2(void) { elips_hash_to_g2(&sink_ep2, h2c_msg, sizeof h2c_msg, (const uint8_t *)"BENCH", 5); }

typedef struct { const char *name; void (*fn)(void); int inner; } bench_t;

static const bench_t BENCHES[] = {
    { "ep_mul",          b_ep_mul,      BENCH_INNER },
    { "ep_mul_glv",      b_ep_mul_glv,  BENCH_INNER },
    { "ep2_mul",         b_ep2_mul,     BENCH_INNER },
    { "ep2_mul_glv",     b_ep2_mul_glv, BENCH_INNER },
    { "gt_exp",          b_gt_exp,      BENCH_INNER },
    { "gt_exp_ct",       b_gt_exp_ct,   BENCH_INNER },
    { "miller",          b_miller,      BENCH_INNER },
    { "final_exp",       b_final_exp,   BENCH_INNER },
    { "pairing",         b_pairing,     BENCH_INNER },
    { "ep_in_subgroup",  b_ep_in_sub,   BENCH_INNER },
    { "ep2_in_subgroup", b_ep2_in_sub,  BENCH_INNER },
    { "hash_to_g1",      b_hash_g1,     BENCH_INNER },
    { "hash_to_g2",      b_hash_g2,     BENCH_INNER },
};
#define N_BENCH ((int)(sizeof BENCHES / sizeof BENCHES[0]))

int main(int argc, char **argv)
{
    int reps = BENCH_REPS, as_json = 0;
    const char *filter = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--json"))      as_json = 1;
        else if (!strcmp(argv[i], "--reps") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (argv[i][0] != '-')          filter = argv[i];
        else {
            fprintf(stderr, "usage: %s [--json] [--reps N] [name]\n", argv[0]);
            return 2;
        }
    }
    if (reps < 3) reps = 3;
    if (reps > BENCH_MAX_REPS) reps = BENCH_MAX_REPS;

    ep_generator(&g1);
    ep2_generator(&g2);
    for (int i = 0; i < BENCH_MAX_REPS; i++) elips_random_scalar(scalars[i]);
    memset(h2c_msg, 0xa5, sizeof h2c_msg);
    /* elips_pairing returns 0 on FAILURE, not on success. */
    if (elips_pairing(gt_base, &g1, &g2) == 0) {
        fprintf(stderr, "bench: pairing failed during setup\n");
        return 1;
    }
    if (ep2_to_affine(qx, qy, &g2) == 0 || ep_to_affine(px, py, &g1) == 0) {
        fprintf(stderr, "bench: a generator is infinity\n");
        return 1;
    }

    if (as_json)
        printf("{\n  \"curve\": \"%s\",\n  \"reps\": %d,\n  \"inner\": %d,\n  \"units\": \"microseconds\",\n  \"results\": {\n",
               ELIPS_CURVE_NAME, reps, BENCH_INNER);
    else
        printf("%s  %d reps of %d calls\n\n%-18s %10s %10s %8s\n",
               ELIPS_CURVE_NAME, reps, BENCH_INNER, "operation", "median", "min", "spread");

    int selected[N_BENCH], n_sel = 0;
    for (int b = 0; b < N_BENCH; b++)
        if (!filter || !strcmp(filter, BENCHES[b].name)) selected[n_sel++] = b;
    if (!n_sel) { fprintf(stderr, "bench: no benchmark matched\n"); return 2; }

    double *samples = malloc((size_t)n_sel * (size_t)reps * sizeof *samples);
    if (!samples) return 1;
#define SAMPLE(sel, rep) samples[(size_t)(sel) * (size_t)reps + (size_t)(rep)]

    /* One untimed pass each: first-touch page faults and cold caches belong to
     * the machine, not to the routine. */
    for (int j = 0; j < n_sel; j++) {
        const bench_t *bn = &BENCHES[selected[j]];
        for (int i = 0; i < bn->inner; i++) bn->fn();
    }

    /* Interleaved: every operation is timed once per rep, in the same order,
     * so machine drift lands on all of them alike. See the note at the top. */
    for (int r = 0; r < reps; r++) {
        for (int j = 0; j < n_sel; j++) {
            const bench_t *bn = &BENCHES[selected[j]];
            double t0 = now_us();
            for (int i = 0; i < bn->inner; i++) bn->fn();
            SAMPLE(j, r) = (now_us() - t0) / bn->inner;
        }
    }

    double *one = malloc((size_t)reps * sizeof *one);
    if (!one) { free(samples); return 1; }
    int printed = 0;

    for (int j = 0; j < n_sel; j++) {
        for (int r = 0; r < reps; r++) one[r] = SAMPLE(j, r);
        qsort(one, (size_t)reps, sizeof *one, cmp_double);

        double med = one[reps / 2], lo = one[0], hi = one[reps - 1];
        double spread = med > 0 ? (hi - lo) / med * 100.0 : 0.0;
        const char *nm = BENCHES[selected[j]].name;

        if (as_json)
            printf("%s    \"%s\": { \"median\": %.3f, \"min\": %.3f, \"max\": %.3f, \"spread_pct\": %.2f }",
                   printed ? ",\n" : "", nm, med, lo, hi, spread);
        else
            printf("%-18s %10.2f %10.2f %7.1f%%\n", nm, med, lo, spread);
        printed++;
    }
    free(one);
#undef SAMPLE

    if (as_json) printf("\n  }\n}\n");
    free(samples);
    return 0;
}
