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
 * A t-test on raw timings is dominated by the tail -- a preemption, a page
 * fault, a migration between cores -- so the measurements are also cropped at a
 * ladder of percentiles and the test re-run on each crop. The reported figure
 * is the largest |t| over all crops, which is the conservative choice.
 *
 * What the numbers mean here. dudect's own ladder is |t| < 5 clean, 5 to 10
 * inconclusive, above 10 leaking -- but that is for ONE t-test, and the figure
 * above is the largest of 21 correlated ones. Reaching 10 to 16 by chance is
 * ordinary for a maximum over 21 tests on a shared runner, so above 10 is
 * treated as a suspicion, re-measured with four times the data, and judged
 * against T_CONFIRM. See the note on T_CONFIRM in main() for how that number
 * was measured rather than chosen.
 *
 * TWO CONTROLS, at opposite ends. control_vartime calls a documented
 * variable-time inversion on secrets and reads in the hundreds: it proves the
 * harness can find an obvious leak. sensitivity plants ONE extra field
 * multiplication on one bit of the secret and reads 30 to 60: it proves the
 * harness can find a small one, and it is what keeps T_CONFIRM honest. Both
 * must report a leak or the run fails.
 *
 * HONEST LIMITS. This is a statistical test, not a proof.
 *   - It cannot prove the absence of leakage, only fail to find it.
 *   - The floor is set by the runner, not by this file. A leak much smaller
 *     than the one the sensitivity control plants would sit inside the noise
 *     band on shared CI hardware and would not be reported.
 *   - It measures wall-clock time. A leak that shows only in cache state or
 *     branch-predictor state and not in total time is invisible here.
 *   - The fixed class must be a representative secret. It was once k = 0,
 *     which held the accumulator at infinity and compared all-zero operands
 *     against random ones -- a hardware question, not a control-flow one.
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

/* Internal, not public API: the reduction that replaced GMP's mpn_sec_div_r.
 * It is reachable because src/ is on the include path for this target. */
#include "arith/wide.h"

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
    limb_t  wide[2 * FP_LIMBS + 2];
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

/* elips_mod_wide, at the width hash_to_field calls it with.
 *
 * This is the routine that replaced GMP's mpn_sec_div_r, and it is the one
 * place in the library where a secret is reduced modulo a public modulus:
 * hash_to_curve runs it on a value derived from the message, which is usually
 * the thing the caller is hiding. If its running time depended on that value,
 * every hash-to-curve in the library would leak.
 *
 * The fixed class is a pinned pseudorandom value, not zero and not one. A
 * degenerate fixed class measures the hardware multiplier rather than the
 * control flow; PROGRESS.md records what that cost to learn. */
#define MODWIDE_NN ((int)((ELIPS_H2C_L + 7) / 8))

static void prep_wide(input_t *in, int c)
{
    if (c) {
        elips_random_bytes(in->wide, (size_t)MODWIDE_NN * sizeof(limb_t));
    } else {
        uint64_t z = 0;
        for (int i = 0; i < MODWIDE_NN; i++) {
            uint64_t x = (z += 0x9e3779b97f4a7c15ULL);
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            in->wide[i] = (limb_t)(x ^ (x >> 31));
        }
    }
}

static limb_t sink_wide[2 * FP_LIMBS + 2];

static void run_mod_wide(const input_t *in)
{
    memcpy(sink_wide, in->wide, (size_t)MODWIDE_NN * sizeof(limb_t));
    elips_mod_wide(sink_wide, MODWIDE_NN, FP_MODULUS, (int)FP_LIMBS);
}

/* The fixed class must be a REPRESENTATIVE secret, not a degenerate one.
 *
 * This used to be k = 0, and that is wrong in a way that took a macOS Release
 * failure to expose. With k = 0 every digit of the GLV decomposition is zero,
 * the ladder selects the identity at every step, and the accumulator stays at
 * infinity from the first doubling to the last. So every field multiplication
 * in class 0 runs on all-zero operands while class 1 runs on random ones.
 *
 * The instruction counts are still identical -- every loop bound in
 * ep2_mul_glv comes from a public bit length, and the 16-entry table scan
 * masks every entry -- so what the t-test was picking up is operand-value
 * dependence in the hardware multiplier and memory system, not a branch in the
 * library. It reported |t| = 57 on Apple silicon at -O2 and stayed under 2.3 on
 * x86-64, which is the signature of a microarchitectural effect rather than a
 * control-flow one.
 *
 * That comparison is not the threat model. The question is whether timing
 * distinguishes one realistic scalar from another, so both classes have to be
 * realistic. A fixed pseudorandom value pinned here gives a reproducible
 * class 0 with an ordinary bit pattern, and a genuine data-dependent branch
 * would still separate it from class 1.
 */
static void fixed_scalar(limb_t *k)
{
    /* splitmix64 from a pinned seed: same value every run, every platform. */
    uint64_t s = 0;
    for (int i = 0; i < ELIPS_ORDER_LIMBS; i++) {
        uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        k[i] = (limb_t)(z ^ (z >> 31));
    }
    /* Clear bit ORDER_BITS-1 and everything above it. The result is less than
     * 2^(ORDER_BITS-1), which is at most r, so the scalar is always in range
     * without needing a reduction. */
    const int top = ELIPS_ORDER_BITS - 1;
    for (int i = 0; i < ELIPS_ORDER_LIMBS; i++) {
        const int lo = i * 64;
        if (lo >= top)            k[i] = 0;
        else if (lo + 64 > top)   k[i] &= ((limb_t)1 << (top - lo)) - 1;
    }
}

/* The SENSITIVITY control, and the calibration for T_CONFIRM.
 *
 * fp_inv_vartime is the other control, and it leaks enormously -- it reads in
 * the hundreds. Passing it only proves the harness can find a leak that nobody
 * would miss. This one is the opposite end: a single extra field
 * multiplication, on one bit of the secret. That is about as small as a
 * data-dependent branch in this library could plausibly be.
 *
 * It reads 30 to 60 and confirms at 39 to 51, against a measured noise ceiling
 * of 15.76 on the CI runner. That gap is what T_CONFIRM = 25 sits in, so this
 * target is not decoration: if the threshold is ever raised past what a
 * one-multiply leak produces, this test fails and says so. */
static void run_sensitivity_probe(const input_t *in)
{
    int extra = (int)(in->k[0] & 1u);
    fp_t acc; fp_copy(acc, in->a);
    for (int i = 0; i < extra; i++) fp_mul(acc, acc, acc);
    fp_copy(sink_fp, acc);
}

static void prep_scalar_and_fp(input_t *in, int c)
{
    if (c) { elips_random_scalar(in->k); fp_rand(in->a); }
    else   { fixed_scalar(in->k); fp_set_one(in->a); }
}

static void prep_scalar(input_t *in, int c)
{
    if (c) elips_random_scalar(in->k);
    else   fixed_scalar(in->k);
}
static void run_ep_mul(const input_t *in)
{ ep_mul(&sink_ep, &g1, in->k, ELIPS_ORDER_BITS); }
static void run_ep2_mul(const input_t *in)
{ ep2_mul(&sink_ep2, &g2, in->k, ELIPS_ORDER_BITS); }
static void run_ep2_glv(const input_t *in)
{ ep2_mul_glv(&sink_ep2, &g2, in->k, ELIPS_ORDER_BITS); }
static void run_ep_glv(const input_t *in)
{ ep_mul_glv(&sink_ep, &g1, in->k, ELIPS_ORDER_BITS); }

/* fp12_exp_gt is the constant-time G_T exponentiation. Its sibling fp12_exp
 * branches on the exponent and is documented as public-exponent only, so it is
 * deliberately NOT a target here: it would leak by design and the harness
 * would be right to say so. */
static fp12_t gt_base, sink_fp12;
static void run_gt_exp(const input_t *in)
{ fp12_exp_gt(sink_fp12, gt_base, in->k, ELIPS_ORDER_BITS); }

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
    { "mod_wide",    prep_wide,    run_mod_wide,    20, 20000, 0 },
    { "ep_mul",      prep_scalar,  run_ep_mul,       1,  3000, 0 },
    { "ep2_mul",     prep_scalar,  run_ep2_mul,      1,  2000, 0 },
    { "ep2_mul_glv", prep_scalar,  run_ep2_glv,      1,  2000, 0 },
    { "ep_mul_glv",  prep_scalar,  run_ep_glv,       1,  3000, 0 },
    { "fp12_exp_gt", prep_scalar,  run_gt_exp,       1,  2000, 0 },
    { "miller",      prep_pairing, run_pairing,      1,  1000, 0 },
    { "hash_to_g1",  prep_h2c,     run_h2c_g1,       1,  1000, 0 },
    { "hash_to_g2",  prep_h2c,     run_h2c_g2,       1,   700, 0 },
    /* n is large because the probed operation is tiny: one field multiply
     * against a scalar multiplication's ~900 us. At n=2000 the reading
     * swung 19 to 79 on a loaded machine and the control itself became
     * flaky, which is the failure it exists to prevent. More samples
     * steady the reading without changing what it detects. */
    { "sensitivity",     prep_scalar_and_fp, run_sensitivity_probe, 1, 40000, 1 },
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

    /* The verdict is taken from a second, four-times-larger sample, against a
     * threshold set from measurement rather than from dudect's nominal 10.
     *
     * Why 10 is the wrong number HERE. dudect's ladder is for a single t-test.
     * max_abs_t reports the largest |t| over 21 correlated crops, and the
     * maximum of 21 tests clears 10 under the null far more often than one
     * test does. On the shared macOS runner this test read 15.76 with a
     * same-size confirm of 10.18 and failed CI, then passed on the very same
     * commit in the next run. Nothing was wrong with the code.
     *
     * Why 25. Both numbers below are measured, not guessed. A deliberately
     * planted leak of ONE field multiplication, on one bit of the secret,
     * inside a ~900 us scalar multiplication reads 40 to 72 -- six out of six
     * runs. The observed noise ceiling is the 15.76 above. 25 sits in that gap,
     * so the gate still catches leaks well under one multiply while the noise
     * band cannot reach it.
     *
     * An escalating-sample rule was tried first and rejected: a real leak's t
     * grows as sqrt(n) only while the effect is small, and these leaks are
     * already saturated at n=2000, so requiring growth reported two genuine
     * leaks out of six as clean. A gate that misses real leaks is worse than
     * one that occasionally cries wolf, so magnitude decides and the larger
     * sample is used only to make that magnitude a steadier reading. */
    const long   CONFIRM_SCALE = 4;
    const double T_CONFIRM = 25.0;

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
            /* Anything over T_LEAK is only a suspicion. Re-measure with four
             * times the data and judge that reading, which is both steadier
             * than the first and taken against a threshold calibrated to this
             * statistic. Repeating at the same n was what used to happen, and
             * it is much weaker than it sounds: a second sample taken a second
             * later on a shared runner sees the same noise burst as the first.
             * See the note on T_CONFIRM above. */
            t2 = measure(tg, n * CONFIRM_SCALE);
            leaking = (t2 > T_CONFIRM);
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
            printf("  [%s] %-16s max|t| = %8.2f (n=%ld)  ->  %8.2f (n=%ld)%s\n",
                   verdict, tg->name, t, n, t2, n * CONFIRM_SCALE,
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
