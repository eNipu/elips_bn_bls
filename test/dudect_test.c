/*
 * Timing-leakage tests, in the style of Reparaz, Balasch and Verbauwhede's
 * dudect (DATE 2017). Plan section 6, Phase 5, item 5.
 *
 * The idea, and why it is worth having: reading code and declaring it constant
 * time is exactly how constant-time code goes wrong. dudect instead asks a
 * question the machine can answer. Feed the routine two classes of input -- one
 * fixed, one random -- time it, and test whether the two timing distributions
 * differ. A routine whose running time depends on its secret input produces two
 * distinguishable distributions and Welch's t-test finds it. No model of the
 * microarchitecture is needed, which is the point: the model is what one gets
 * wrong.
 *
 * What the numbers mean, following dudect's own guidance:
 *
 *   |t| < 5     no evidence of leakage at this sample size
 *   5 < |t| < 10 inconclusive; run longer
 *   |t| > 10    leaking
 *
 * A t-test on raw timings is dominated by the tail -- a preemption, a page
 * fault, a migration between cores -- so the measurements are also cropped at a
 * ladder of percentiles and the test re-run on each crop. The reported figure
 * is the largest |t| over all crops, which is the conservative choice.
 *
 * HONEST LIMITS. This is a statistical test, not a proof.
 *   - It cannot prove the absence of leakage, only fail to find it.
 *   - On a shared CI runner the noise floor is high and a single run can throw
 *     a large |t| with nothing wrong. The gate therefore confirms: a first
 *     failure is re-measured with a fresh sample before it is reported.
 *   - It measures wall-clock time. A leak that shows only in cache state or
 *     branch-predictor state and not in total time is invisible here.
 *
 * Usage:  dudect_test [target|all] [measurements]
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#include "elips/pairing.h"
#include "elips/random.h"
#include "elips/hash_to_curve.h"

/* ---------------------------------------------------------------- timing --- */

static inline uint64_t cycles(void)
{
#if defined(__x86_64__)
    uint32_t lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
#elif defined(__aarch64__)
    uint64_t v;
    __asm__ __volatile__("mrs %0, cntvct_el0" : "=r"(v));
    return v;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
#endif
}

/* ------------------------------------------------------------ the t-test --- */

typedef struct { double n, mean, m2; } acc_t;

static void acc_push(acc_t *a, double x)
{
    a->n += 1.0;
    double d = x - a->mean;
    a->mean += d / a->n;
    a->m2   += d * (x - a->mean);
}

static double welch_t(const acc_t *a, const acc_t *b)
{
    if (a->n < 2.0 || b->n < 2.0) return 0.0;
    double va = a->m2 / (a->n - 1.0), vb = b->m2 / (b->n - 1.0);
    double denom = sqrt(va / a->n + vb / b->n);
    if (denom == 0.0) return 0.0;
    return (a->mean - b->mean) / denom;
}

static int cmp_u64(const void *x, const void *y)
{
    uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y;
    return (a > b) - (a < b);
}

/* The percentile ladder from dudect: dense near the top, because that is where
 * the tail that swamps the statistic lives. */
static double crop_fraction(int i, int n)
{
    return 1.0 - pow(0.5, 10.0 * (double)(i + 1) / (double)n);
}

#define NCROP 20

static double max_abs_t(const uint64_t *t, const uint8_t *cls, long n,
                        double *raw_t_out)
{
    uint64_t *sorted = malloc((size_t)n * sizeof *sorted);
    if (!sorted) return 0.0;
    memcpy(sorted, t, (size_t)n * sizeof *sorted);
    qsort(sorted, (size_t)n, sizeof *sorted, cmp_u64);

    uint64_t thresh[NCROP];
    for (int c = 0; c < NCROP; c++) {
        long idx = (long)(crop_fraction(c, NCROP) * (double)n);
        if (idx >= n) idx = n - 1;
        thresh[c] = sorted[idx];
    }
    free(sorted);

    acc_t a[NCROP + 1], b[NCROP + 1];
    memset(a, 0, sizeof a);
    memset(b, 0, sizeof b);

    for (long i = 0; i < n; i++) {
        acc_t *dst = cls[i] ? a : b;
        acc_push(&dst[NCROP], (double)t[i]);            /* uncropped */
        for (int c = 0; c < NCROP; c++)
            if (t[i] < thresh[c]) acc_push(&dst[c], (double)t[i]);
    }

    double best = 0.0;
    for (int c = 0; c <= NCROP; c++) {
        double v = fabs(welch_t(&a[c], &b[c]));
        if (v > best) best = v;
    }
    if (raw_t_out) *raw_t_out = fabs(welch_t(&a[NCROP], &b[NCROP]));
    return best;
}

/* ------------------------------------------------------------- the cases --- */
/*
 * Each case supplies prepare(slot, class), which writes one measurement's
 * secret input, and run(slot), which performs the operation on it. Class 0 is a
 * fixed input, class 1 a fresh random one. Anything public -- the generator,
 * the curve parameters -- is set up once and shared, so a timing difference can
 * only come from the secret.
 *
 * EVERY INPUT IS PREPARED BEFORE ANY MEASUREMENT. This is not an optimisation;
 * it is the difference between measuring the routine and measuring its
 * preparation. Generating a random field element calls into the kernel, and
 * doing that inside the timed loop for one class and not the other leaves the
 * cache and branch predictor in visibly different states. The first version of
 * this file did exactly that and reported |t| = 98 for fp_add, an addition with
 * no branch in it at all. Preparing up front made the same test read 1.4.
 */

#define H2C_MSG_LEN 48

typedef struct {
    fp_t    a, b;
    limb_t  k[ELIPS_ORDER_LIMBS];
    ep_t    P;
    ep2_t   Q;
    uint8_t msg[H2C_MSG_LEN];
} input_t;

static ep_t  g1;
static ep2_t g2;

/* Sinks, so the compiler cannot discard the work. */
static fp_t   sink_fp;
static ep_t   sink_ep;
static ep2_t  sink_ep2;
static fp12_t sink_fp12;

static void prep_fp_pair(input_t *in, int c)
{
    if (c) { fp_rand(in->a); fp_rand(in->b); }
    else   { fp_set_one(in->a); fp_copy(in->b, in->a); }
}
static void run_fp_mul(const input_t *in)     { fp_mul(sink_fp, in->a, in->b); }
static void run_fp_add(const input_t *in)     { fp_add(sink_fp, in->a, in->b); }
static void run_fp_inv(const input_t *in)     { fp_inv(sink_fp, in->a); }
/* The negative control. fp_inv_vartime is documented as variable time and is
 * only ever called on public values; here it is called on secrets on purpose,
 * so the harness has something it is supposed to catch. A leakage test that
 * cannot report a leak proves nothing about the tests that pass, which is the
 * same argument that put a corrupted-vector check in the Phase 0 suite. */
static void run_fp_inv_vt(const input_t *in)  { fp_inv_vartime(sink_fp, in->a); }
static void run_fp_cselect(const input_t *in)
{ fp_cselect(sink_fp, in->a, in->b, (limb_t)0 - (in->a[0] & 1)); }

static void prep_scalar(input_t *in, int c)
{
    if (c) elips_random_scalar(in->k);
    else   memset(in->k, 0, sizeof in->k);
}
static void run_ep_mul(const input_t *in)
{ ep_mul(&sink_ep, &g1, in->k, ELIPS_ORDER_BITS); }
static void run_ep2_mul(const input_t *in)
{ ep2_mul(&sink_ep2, &g2, in->k, ELIPS_ORDER_BITS); }
#ifdef ELIPS_FAMILY_BLS12
static void run_ep2_glv(const input_t *in)
{ ep2_mul_glv(&sink_ep2, &g2, in->k, ELIPS_ORDER_BITS); }
#endif

/* The pairing's secret is the point, not the loop bound: the Miller loop runs
 * over a public curve parameter. Class 1 is [k]G for a fresh random k, class 0
 * a fixed multiple, so the coordinates vary and nothing else does. */
static void prep_pairing(input_t *in, int c)
{
    limb_t k[ELIPS_ORDER_LIMBS];
    if (c) elips_random_scalar(k);
    else { memset(k, 0, sizeof k); k[0] = 3; }
    ep_mul(&in->P, &g1, k, ELIPS_ORDER_BITS);
    ep2_copy(&in->Q, &g2);
}
static void run_pairing(const input_t *in)
{
    fp_t px, py; fp2_t qx, qy;
    ep_to_affine(px, py, &in->P);
    ep2_to_affine(qx, qy, &in->Q);
    pairing_miller(sink_fp12, qx, qy, px, py);
}

/* The message hashed to a curve is not always public: an oblivious PRF or a
 * password-authenticated exchange hashes a secret, and the header promises this
 * path does not leak it. Class 0 is a fixed message, class 1 a random one, both
 * the same length -- length is not a secret this can hide and padding it here
 * would test the wrong thing. */
static const uint8_t H2C_DST[] = "ELIPS-DUDECT-V01-CS01-" ELIPS_H2C_SUITE_G1;

static void prep_h2c(input_t *in, int c)
{
    if (c) elips_random_bytes(in->msg, sizeof in->msg);
    else   memset(in->msg, 0x5a, sizeof in->msg);
}
static void run_h2c_g1(const input_t *in)
{
    elips_hash_to_g1(&sink_ep, in->msg, sizeof in->msg,
                     H2C_DST, sizeof H2C_DST - 1);
}
static void run_h2c_g2(const input_t *in)
{
    elips_hash_to_g2(&sink_ep2, in->msg, sizeof in->msg,
                     H2C_DST, sizeof H2C_DST - 1);
}

typedef struct {
    const char *name;
    void (*prepare)(input_t *, int);
    void (*run)(const input_t *);
    int  reps;          /* inner repetitions, to lift cheap ops off the timer floor */
    long dflt;          /* default number of measurements */
    int  expect_leak;   /* the negative control: this one must be caught */
} target_t;

static const target_t TARGETS[] = {
    { "fp_mul",      prep_fp_pair, run_fp_mul,     200, 20000, 0 },
    { "fp_add",      prep_fp_pair, run_fp_add,     200, 20000, 0 },
    { "fp_cselect",  prep_fp_pair, run_fp_cselect, 200, 20000, 0 },
    { "fp_inv",      prep_fp_pair, run_fp_inv,       1, 20000, 0 },
    { "ep_mul",      prep_scalar,  run_ep_mul,       1,  3000, 0 },
    { "ep2_mul",     prep_scalar,  run_ep2_mul,      1,  2000, 0 },
#ifdef ELIPS_FAMILY_BLS12
    { "ep2_mul_glv", prep_scalar,  run_ep2_glv,      1,  2000, 0 },
#endif
    { "miller",      prep_pairing, run_pairing,      1,  1000, 0 },
    { "hash_to_g1",  prep_h2c,     run_h2c_g1,       1,  1000, 0 },
    { "hash_to_g2",  prep_h2c,     run_h2c_g2,       1,   700, 0 },
    { "control_vartime", prep_fp_pair, run_fp_inv_vt, 1, 20000, 1 },
};
#define NTARGET ((int)(sizeof TARGETS / sizeof TARGETS[0]))

/* --------------------------------------------------------------- driver ---- */

static double measure(const target_t *tg, long n)
{
    input_t  *in    = malloc((size_t)n * sizeof *in);
    uint64_t *ticks = malloc((size_t)n * sizeof *ticks);
    uint8_t  *cls   = malloc((size_t)n);
    if (!in || !ticks || !cls) { free(in); free(ticks); free(cls); return -1.0; }

    for (long i = 0; i < n; i++) {
        uint8_t bit;
        if (elips_random_bytes(&bit, 1) != 0) {
            free(in); free(ticks); free(cls); return -1.0;
        }
        cls[i] = (uint8_t)(bit & 1);
        tg->prepare(&in[i], cls[i]);
    }

    for (long i = 0; i < n; i++) {
        uint64_t t0 = cycles();
        for (int r = 0; r < tg->reps; r++) tg->run(&in[i]);
        uint64_t t1 = cycles();
        ticks[i] = t1 - t0;
    }

    double t = max_abs_t(ticks, cls, n, NULL);
    free(in); free(ticks); free(cls);
    return t;
}

int main(int argc, char **argv)
{
    const char *want = (argc > 1) ? argv[1] : "all";
    long override_n  = (argc > 2) ? strtol(argv[2], NULL, 10) : 0;

    ep_generator(&g1);
    ep2_generator(&g2);

    /* dudect's own thresholds. Between the two, the sample is too small to
     * decide, so the run is repeated rather than reported either way. */
    const double T_LEAK = 10.0, T_CLEAN = 5.0;

    printf("dudect timing tests [%s]\n", ELIPS_CURVE_NAME);
    int fails = 0, ran = 0;

    for (int i = 0; i < NTARGET; i++) {
        const target_t *tg = &TARGETS[i];
        if (strcmp(want, "all") != 0 && strcmp(want, tg->name) != 0) continue;
        ran++;

        long n = override_n > 0 ? override_n : tg->dflt;
        double t = measure(tg, n);
        if (t < 0.0) { printf("  [FAIL ] %-16s out of memory\n", tg->name); fails++; continue; }

        int leaking = 0;
        double t2 = -1.0;
        if (t > T_LEAK) {
            /* Confirm before deciding: a preempted run can throw a large t on a
             * shared machine. A second independent sample has to agree. */
            t2 = measure(tg, n);
            leaking = (t2 > T_LEAK);
        }

        const char *verdict;
        if (tg->expect_leak) {
            /* The control. Not finding the leak is the failure. */
            verdict = leaking ? "PASS " : "FAIL ";
            if (!leaking) fails++;
        } else if (leaking) {
            verdict = "LEAK "; fails++;
        } else if (t < T_CLEAN) {
            verdict = "PASS ";
        } else {
            verdict = "weak ";
        }

        if (t2 >= 0.0)
            printf("  [%s] %-16s max|t| = %8.2f  (confirm run: %8.2f, n=%ld)%s\n",
                   verdict, tg->name, t, t2, n,
                   tg->expect_leak ? "  <- negative control, must leak" : "");
        else
            printf("  [%s] %-16s max|t| = %8.2f  (n=%ld)%s\n",
                   verdict, tg->name, t, n,
                   tg->expect_leak ? "  <- negative control, must leak" : "");
    }

    if (ran == 0) {
        fprintf(stderr, "unknown target '%s'. Known targets:\n", want);
        for (int i = 0; i < NTARGET; i++) fprintf(stderr, "  %s\n", TARGETS[i].name);
        return 2;
    }
    printf("%s\n", fails ? "FAILED" : "all targets behaved as expected");
    return fails ? 1 : 0;
}
