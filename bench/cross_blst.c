/* P0: like-for-like ELiPS vs blst.
 *
 * The question this exists to answer: the recorded 2.6x pairing ratio compared
 * elips_pairing, which validates both input points, against a blst number that
 * may not have. Validation is 21% of the ELiPS pairing, so the ratio may be
 * measuring an API difference rather than an arithmetic one.
 *
 * Interleaved reps and a median, copied from bench/bench.c for the reason
 * stated there: back-to-back reps of one case share a frequency and a cache
 * state, so a non-interleaved spread understates the real error badly.
 *
 * Every case goes through the same function-pointer table and the same timing
 * loop. That is deliberate: an earlier measurement in this project put five
 * cases behind one switch and penalised the cheapest of them.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "elips/ec.h"
#include "elips/fpx.h"
#include "elips/pairing.h"
#include "blst.h"

#define REPS 21
#define INNER 30

static double now_us(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1e6 + (double)t.tv_nsec * 1e-3;
}
static int cmpd(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* --- ELiPS state --- */
static ep_t   e_g1;
static ep2_t  e_g2;
static fp2_t  e_qx, e_qy;
static fp_t   e_px, e_py;
static fp12_t e_sink, e_gt;
static volatile unsigned e_acc;

static void e_miller(void)    { pairing_miller(e_sink, e_qx, e_qy, e_px, e_py); }
static void e_fexp(void)      { pairing_final_exp_fast(e_sink, e_gt); }
static void e_arith(void)     { fp12_t f; pairing_miller(f, e_qx, e_qy, e_px, e_py);
                                pairing_final_exp_fast(e_sink, f); }
static void e_full(void)      { e_acc += (unsigned)elips_pairing(e_sink, &e_g1, &e_g2); }
static void e_in_g1(void)     { e_acc += (unsigned)ep_in_subgroup(&e_g1); }
static void e_in_g2(void)     { e_acc += (unsigned)ep2_in_subgroup(&e_g2); }

/* --- blst state --- */
static blst_p1_affine b_p;
static blst_p2_affine b_q;
static blst_fp12 b_sink, b_gt;
static volatile unsigned b_acc;

static void b_miller(void)    { blst_miller_loop(&b_sink, &b_q, &b_p); }
static void b_fexp(void)      { blst_final_exp(&b_sink, &b_gt); }
static void b_arith(void)     { blst_fp12 f; blst_miller_loop(&f, &b_q, &b_p);
                                blst_final_exp(&b_sink, &f); }
static void b_in_g1(void)     { b_acc += (unsigned)blst_p1_affine_in_g1(&b_p); }
static void b_in_g2(void)     { b_acc += (unsigned)blst_p2_affine_in_g2(&b_q); }
static void b_full(void)      { blst_fp12 f;
                                b_acc += (unsigned)blst_p1_affine_in_g1(&b_p);
                                b_acc += (unsigned)blst_p2_affine_in_g2(&b_q);
                                blst_miller_loop(&f, &b_q, &b_p);
                                blst_final_exp(&b_sink, &f); }

typedef struct { const char *name; void (*fn)(void); } case_t;
static const case_t CASES[] = {
    { "elips.miller",          e_miller },
    { "elips.final_exp",       e_fexp   },
    { "elips.miller+fexp",     e_arith  },
    { "elips.pairing_full",    e_full   },
    { "elips.in_g1",           e_in_g1  },
    { "elips.in_g2",           e_in_g2  },
    { "blst.miller",           b_miller },
    { "blst.final_exp",        b_fexp   },
    { "blst.miller+fexp",      b_arith  },
    { "blst.pairing_full",     b_full   },
    { "blst.in_g1",            b_in_g1  },
    { "blst.in_g2",            b_in_g2  },
};
#define NC ((int)(sizeof CASES / sizeof CASES[0]))

int main(int argc, char **argv)
{
    int reps = REPS, as_json = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--json")) as_json = 1;
        else if (!strcmp(argv[i], "--reps") && i + 1 < argc) reps = atoi(argv[++i]);
    }

    ep_generator(&e_g1);
    ep2_generator(&e_g2);
    if (!ep_to_affine(e_px, e_py, &e_g1) || !ep2_to_affine(e_qx, e_qy, &e_g2)) {
        fprintf(stderr, "elips: affine conversion failed\n"); return 1;
    }
    if (elips_pairing(e_gt, &e_g1, &e_g2) == 0) {
        fprintf(stderr, "elips: pairing failed in setup\n"); return 1;
    }
    b_p = *blst_p1_affine_generator();
    b_q = *blst_p2_affine_generator();
    blst_miller_loop(&b_gt, &b_q, &b_p);
    blst_final_exp(&b_gt, &b_gt);

    static double s[NC][512];
    for (int r = 0; r < reps; r++)
        for (int c = 0; c < NC; c++) {
            double t0 = now_us();
            for (int k = 0; k < INNER; k++) CASES[c].fn();
            s[c][r] = (now_us() - t0) / INNER;
        }

    if (as_json) printf("{\n");
    for (int c = 0; c < NC; c++) {
        qsort(s[c], (size_t)reps, sizeof(double), cmpd);
        double med = s[c][reps / 2], lo = s[c][0], hi = s[c][reps - 1];
        if (as_json)
            printf("  \"%s\": {\"median\": %.3f, \"min\": %.3f, \"spread_pct\": %.1f}%s\n",
                   CASES[c].name, med, lo, 100.0 * (hi - lo) / med, c + 1 < NC ? "," : "");
        else
            printf("%-22s median %9.3f us   min %9.3f   spread %5.1f%%\n",
                   CASES[c].name, med, lo, 100.0 * (hi - lo) / med);
    }
    if (as_json) printf("}\n");
    return (int)(e_acc & b_acc & 0u);
}
