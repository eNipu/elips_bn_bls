#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "kernels.h"

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); exit(1); \
} } while (0)

typedef struct { fp2_t B, C; } normalized_line;
typedef struct { normalized_line l[ELIPS_MILLER_LINES]; } normalized_table;
typedef void (*square_fn)(fp12_t, const fp12_t);
enum { MODES = 4, REPS = 31 };
static square_fn methods[] = {fp12_sqr, gaussian_sqr, cubic_sqr, fft4_sqr};
static const char *names[] = {
    "repository-square", "gaussian-scaled4", "cubic-scaled2", "fft4-scaled4"
};
static uint64_t rng_state = UINT64_C(0x12dbba82178a449);
static volatile limb_t sink;

static uint64_t random_word(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}
static void random_fp(fp_t r)
{
    limb_t a[FP_LIMBS];
    int less;
    do {
        for (int i = 0; i < FP_LIMBS; i++) a[i] = random_word();
        a[FP_LIMBS-1] &= UINT64_MAX >> (64*FP_LIMBS-FP_BITS);
        less = 0;
        for (int i = FP_LIMBS-1; i >= 0; i--)
            if (a[i] != FP_MODULUS[i]) { less = a[i] < FP_MODULUS[i]; break; }
    } while (!less);
    fp_from_limbs(r, a);
}
static void random12(fp12_t r)
{
    for (int i = 0; i < 2; i++) for (int j = 0; j < 3; j++)
        for (int k = 0; k < 2; k++) random_fp(r[i][j][k]);
}
static void normalize(normalized_table *out, const ep2_prec_t *pc)
{
    fp2_t prefix[ELIPS_MILLER_LINES], inv, ai;
    CHECK(!fp2_is_zero(pc->l[0].a));
    fp2_copy(prefix[0], pc->l[0].a);
    for (int i = 1; i < ELIPS_MILLER_LINES; i++) {
        CHECK(!fp2_is_zero(pc->l[i].a));
        fp2_mul(prefix[i], prefix[i-1], pc->l[i].a);
    }
    fp2_inv(inv, prefix[ELIPS_MILLER_LINES-1]);
    for (int i = ELIPS_MILLER_LINES-1; i >= 0; i--) {
        if (i) fp2_mul(ai, inv, prefix[i-1]); else fp2_copy(ai, inv);
        fp2_mul(out->l[i].B, pc->l[i].c, ai);
        fp2_mul(out->l[i].C, pc->l[i].b, ai);
        if (i) fp2_mul(inv, inv, pc->l[i].a);
    }
}
static void apply(fp12_t f, const normalized_line *line,
                  const fp_t iy, const fp_t xy)
{
    fp2_t B, C;
    fp2_mul_fp(B, line->B, iy); fp2_mul_fp(C, line->C, xy);
    normalized8(f, B, C);
}
static void miller(fp12_t f, const normalized_table *nt,
                   const fp_t *px, const fp_t *py, int n, int mode)
{
    fp_t iy[8], xy[8];
    CHECK(n >= 0 && n <= 8);
    if (!n) { fp12_set_one(f); return; }
    for (int j = 0; j < n; j++) {
        fp_inv(iy[j], py[j]); fp_mul(xy[j], px[j], iy[j]);
    }
    int pos = 0;
    fp12_set_one(f);
    for (int i = ELIPS_LOOP_TOP-1; i >= 0; i--) {
        methods[mode](f, f);
        for (int j = 0; j < n; j++) {
            apply(f, &nt[j].l[pos], iy[j], xy[j]);
            if (ELIPS_LOOP[i]) apply(f, &nt[j].l[pos+1], iy[j], xy[j]);
        }
        pos += 1 + (ELIPS_LOOP[i] != 0);
    }
    CHECK(pos == ELIPS_MILLER_LINES);
}
static double clock_us(void)
{
    struct timespec ts;
    CHECK(clock_gettime(CLOCK_MONOTONIC, &ts) == 0);
    return ts.tv_sec*1e6 + ts.tv_nsec/1e3;
}
static int compare_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static void print_times(double times[MODES][REPS], const char *kind, int n)
{
    for (int rep = 0; rep < REPS; rep++)
        for (int mode = 0; mode < MODES; mode++)
            printf("SAMPLE,%s,%d,%d,%s,%.6f\n", kind, n, rep,
                   names[mode], times[mode][rep]);
    for (int mode = 0; mode < MODES; mode++) {
        qsort(times[mode], REPS, sizeof(double), compare_double);
        printf("BENCH %s n=%d %s median=%.3f us min=%.3f max=%.3f\n",
               kind, n, names[mode], times[mode][REPS/2],
               times[mode][0], times[mode][REPS-1]);
    }
}
static void bench_squares(void)
{
    double times[MODES][REPS];
    fp12_t seed, f;
    random12(seed);
    for (int rep = -2; rep < REPS; rep++) for (int k = 0; k < MODES; k++) {
        int mode = (k + rep + MODES) % MODES;
        fp12_copy(f, seed);
        double begin = clock_us();
        for (int i = 0; i < 4000; i++) methods[mode](f, f);
        double elapsed = (clock_us() - begin) / 4000;
        sink ^= f[0][0][0][0];
        if (rep >= 0) times[mode][rep] = elapsed;
    }
    print_times(times, "square", 0);
}
static void bench_miller(const normalized_table *nt, const fp_t *px,
                         const fp_t *py, int n, int final)
{
    double times[MODES][REPS];
    fp12_t f, out;
    for (int rep = -2; rep < REPS; rep++) for (int k = 0; k < MODES; k++) {
        int mode = (k + rep + MODES) % MODES;
        double begin = clock_us();
        for (int i = 0; i < 30; i++) {
            miller(f, nt, px, py, n, mode);
            if (final) { pairing_final_exp_fast(out, f); sink ^= out[0][0][0][0]; }
            else sink ^= f[0][0][0][0];
        }
        double elapsed = (clock_us() - begin) / 30;
        if (rep >= 0) times[mode][rep] = elapsed;
    }
    print_times(times, final ? "miller+final" : "miller", n);
}
int main(int argc, char **argv)
{
    (void)argv;
    fp12_t f, expected, out, inplace;
    for (int i = 0; i < 4000; i++) {
        random12(f);
        if (i % 101 == 0) fp12_set_zero(f);
        if (i % 103 == 0) fp12_set_one(f);
        if (i % 107 == 0) {
            fp2_copy(f[0][2], f[0][0]); fp2_neg(f[1][0], f[0][0]);
            fp2_set_zero(f[1][2]);
        }
        fp12_mul(expected, f, f);
        for (int mode = 0; mode < MODES; mode++) {
            fp12_t scaled;
            fp12_copy(scaled, expected);
            if (mode) fp12_add(scaled, scaled, scaled);
            if (mode == 1 || mode == 3) fp12_add(scaled, scaled, scaled);
            methods[mode](out, f); CHECK(fp12_eq(out, scaled));
            fp12_copy(inplace, f); methods[mode](inplace, inplace);
            CHECK(fp12_eq(inplace, scaled));
        }
    }
    puts("PASS 4000 full-width/edge square checks, all methods, both alias modes");
    ep_t G, P;
    ep2_t Q, R;
    ep_generator(&G); ep2_generator(&Q);
    ep2_prec_t pc[8];
    normalized_table nt[8];
    fp_t px[8], py[8];
    fp2_t qx[8], qy[8];
    for (int trial = 0; trial < 3; trial++) {
        for (int j = 0; j < 8; j++) {
            /* Full-width public test scalars, reduced by group multiplication. */
            limb_t kp[4], kq[4];
            for (int i = 0; i < 4; i++) { kp[i] = random_word(); kq[i] = random_word(); }
            ep_mul(&P, &G, kp, 256); ep2_mul(&R, &Q, kq, 256);
            CHECK(ep_to_affine(px[j], py[j], &P));
            CHECK(ep2_to_affine(qx[j], qy[j], &R));
            CHECK(ep2_precompute(&pc[j], &R)); normalize(&nt[j], &pc[j]);
        }
        for (int n = 0; n <= 8; n++) {
            pairing_miller_multi(f, qx, qy, px, py, (size_t)n);
            pairing_final_exp_fast(expected, f);
            for (int mode = 0; mode < MODES; mode++) {
                miller(f, nt, px, py, n, mode);
                pairing_final_exp_fast(out, f); CHECK(fp12_eq(out, expected));
            }
        }
    }
    puts("PASS 27 multi-pairing comparisons n=0..8, all methods");
    pairing_miller_multi(f, qx, qy, px, py, 2);
    pairing_final_exp_plain(expected, f);
    for (int mode = 0; mode < MODES; mode++) {
        miller(f, nt, px, py, 2, mode);
        pairing_final_exp_plain(out, f); CHECK(fp12_eq(out, expected));
    }
    puts("PASS exact final-exponent comparison, all methods");
    if (argc > 1) {
        bench_squares();
        const int counts[] = {1, 2, 8};
        for (int i = 0; i < 3; i++)
            for (int final = 0; final < 2; final++)
                bench_miller(nt, px, py, counts[i], final);
    }
    return 0;
}
