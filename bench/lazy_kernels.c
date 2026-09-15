#define _POSIX_C_SOURCE 200809L
/* MEASUREMENT ARTIFACT, not shipped code. See PAIRING-PERFORMANCE.md, "Lazy
 * reduction, measured and rejected".
 *
 * Kept in its own translation unit on purpose. With the probe and the kernels
 * in one file gcc hoists the pure calls out of the timing loop and reports a
 * Montgomery multiply at a few cycles; that happened in this repository during
 * the IFMA work and nearly shipped as a conclusion.
 *
 * Lazy Fp6 multiply prototype. Not integrated: this exists to answer one
 * question, whether amortising reductions 3:1 at Fp6 beats the 13% that
 * splitting the fused CIOS multiply costs.
 *
 * Every wide value is kept in [0, p*R) by conditionally adding or subtracting
 * p*R after each wide add or sub, which is what blst's add_mod_384x384 and
 * sub_mod_384x384 do. That makes the bias bookkeeping trivial and keeps
 * redcw's precondition true by construction: working modulo p*R is exact,
 * because REDC(w + p*R) = w*R^-1 + p == w*R^-1 (mod p). */
#include <string.h>
#include "elips/fp.h"
#include "elips/fpx.h"

void elips_fp_mulw_6_x86_64(limb_t w[12], const limb_t a[6], const limb_t b[6]);
void elips_fp_redcw_6_x86_64(limb_t r[6], const limb_t w[12], const limb_t p[6], limb_t n0);

typedef limb_t fpw_t[12];
static limb_t PR[12];                       /* p * 2^384 */

void fpw_init(void)
{
    memset(PR, 0, sizeof PR);
    for (int i = 0; i < 6; i++) PR[i+6] = FP_MODULUS[i];
}

/* p*R has six ZERO low limbs, so reducing mod p*R can never borrow out of the
 * low half: only the top six limbs ever need the conditional add or subtract.
 * That turns a 768-bit modular add into a 768-bit add plus a 384-bit cmov,
 * which is what makes the whole scheme affordable. */
static inline void fpw_addm(fpw_t r, const fpw_t a, const fpw_t b)
{
    limb_t c = 0, br = 0, hi[6], cand[6];
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)a[i] + b[i] + c;
        r[i] = (limb_t)s; c = (limb_t)(s >> 64);
    }
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)a[i+6] + b[i+6] + c;
        hi[i] = (limb_t)s; c = (limb_t)(s >> 64);
    }
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)hi[i] - FP_MODULUS[i] - br;
        cand[i] = (limb_t)s; br = (limb_t)((s >> 64) & 1);
    }
    limb_t keep = (limb_t)0 - (limb_t)(br & (c ^ 1));   /* below p*R: keep hi */
    for (int i = 0; i < 6; i++) r[i+6] = (hi[i] & keep) | (cand[i] & ~keep);
}

static inline void fpw_subm(fpw_t r, const fpw_t a, const fpw_t b)
{
    limb_t br = 0, c = 0, hi[6], cand[6];
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)a[i] - b[i] - br;
        r[i] = (limb_t)s; br = (limb_t)((s >> 64) & 1);
    }
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)a[i+6] - b[i+6] - br;
        hi[i] = (limb_t)s; br = (limb_t)((s >> 64) & 1);
    }
    for (int i = 0; i < 6; i++) {
        unsigned __int128 s = (unsigned __int128)hi[i] + FP_MODULUS[i] + c;
        cand[i] = (limb_t)s; c = (limb_t)(s >> 64);
    }
    limb_t take = (limb_t)0 - br;                       /* went negative */
    for (int i = 0; i < 6; i++) r[i+6] = (cand[i] & take) | (hi[i] & ~take);
}

/* wide Fp2: two 768-bit coefficients */
typedef fpw_t fp2w_t[2];

static void fp2_mulw(fp2w_t c, const fp2_t a, const fp2_t b)
{
    fpw_t t0, t1, s; fp_t as, bs;
    elips_fp_mulw_6_x86_64(t0, a[0], b[0]);
    elips_fp_mulw_6_x86_64(t1, a[1], b[1]);
    fp_add(as, a[0], a[1]); fp_add(bs, b[0], b[1]);
    elips_fp_mulw_6_x86_64(s, as, bs);
    fpw_subm(c[0], t0, t1);
    fpw_subm(s, s, t0); fpw_subm(c[1], s, t1);
}

static inline void fp2w_add(fp2w_t r, const fp2w_t a, const fp2w_t b)
{ fpw_addm(r[0], a[0], b[0]); fpw_addm(r[1], a[1], b[1]); }
static inline void fp2w_sub(fp2w_t r, const fp2w_t a, const fp2w_t b)
{ fpw_subm(r[0], a[0], b[0]); fpw_subm(r[1], a[1], b[1]); }
/* xi = 1 + u */
static inline void fp2w_mul_xi(fp2w_t r, const fp2w_t a)
{ fpw_t t; fpw_subm(t, a[0], a[1]); fpw_addm(r[1], a[0], a[1]); memcpy(r[0], t, sizeof t); }

static inline void fp2w_redc(fp2_t r, const fp2w_t a)
{ elips_fp_redcw_6_x86_64(r[0], a[0], FP_MODULUS, FP_MONT_N0);
  elips_fp_redcw_6_x86_64(r[1], a[1], FP_MODULUS, FP_MONT_N0); }

void fp6_mul_lazy(fp6_t r, const fp6_t a, const fp6_t b)
{
    fp2w_t v0, v1, v2, w01, w02, w12, t, u;
    fp2_t s1, s2;

    fp2_mulw(v0, a[0], b[0]);
    fp2_mulw(v1, a[1], b[1]);
    fp2_mulw(v2, a[2], b[2]);
    fp2_add(s1, a[1], a[2]); fp2_add(s2, b[1], b[2]); fp2_mulw(w12, s1, s2);
    fp2_add(s1, a[0], a[1]); fp2_add(s2, b[0], b[1]); fp2_mulw(w01, s1, s2);
    fp2_add(s1, a[0], a[2]); fp2_add(s2, b[0], b[2]); fp2_mulw(w02, s1, s2);

    fp2w_sub(t, w12, v1); fp2w_sub(t, t, v2);
    fp2w_mul_xi(u, t); fp2w_add(t, u, v0);
    fp2w_redc(r[0], t);

    fp2w_sub(t, w01, v0); fp2w_sub(t, t, v1);
    fp2w_mul_xi(u, v2); fp2w_add(t, t, u);
    fp2w_redc(r[1], t);

    fp2w_sub(t, w02, v0); fp2w_sub(t, t, v2);
    fp2w_add(t, t, v1);
    fp2w_redc(r[2], t);
}

/* ---- component timing, so the cost is attributed rather than inferred ---- */
#include <time.h>
static double ns_(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec*1e9+(double)t.tv_nsec;}
void fp6_components(const fp6_t a, const fp6_t b, long n,
                    double *addm, double *mulw3, double *fp2mulw, double *redc2)
{
    fpw_t x, y, z; fp2w_t c; fp2_t r; volatile limb_t s=0; double t0;
    elips_fp_mulw_6_x86_64(x, a[0][0], b[0][0]);
    elips_fp_mulw_6_x86_64(y, a[1][0], b[1][0]);
    t0=ns_(); for(long i=0;i<n;i++){ fpw_addm(z,x,y); s^=z[0]; } *addm=(ns_()-t0)/n;
    t0=ns_(); for(long i=0;i<n;i++){ elips_fp_mulw_6_x86_64(x,a[0][0],b[0][0]);
                                     elips_fp_mulw_6_x86_64(y,a[0][1],b[0][1]);
                                     elips_fp_mulw_6_x86_64(z,a[1][0],b[1][0]); s^=x[0]^y[0]^z[0]; }
    *mulw3=(ns_()-t0)/n;
    t0=ns_(); for(long i=0;i<n;i++){ fp2_mulw(c,a[0],b[0]); s^=c[0][0]; } *fp2mulw=(ns_()-t0)/n;
    fp2_mulw(c,a[0],b[0]);
    t0=ns_(); for(long i=0;i<n;i++){ fp2w_redc(r,c); s^=r[0][0]; } *redc2=(ns_()-t0)/n;
    (void)s;
}

void fp6_compare(const fp6_t a, const fp6_t b, long n, double *eager, double *lazy)
{
    fp6_t x, y; volatile limb_t s = 0; double t0;
    fp6_mul(x, a, b); fp6_mul_lazy(y, a, b);
    t0 = ns_(); for (long i = 0; i < n; i++) { fp6_mul(x, a, b); s ^= x[0][0][0]; }
    *eager = (ns_() - t0) / (double)n;
    t0 = ns_(); for (long i = 0; i < n; i++) { fp6_mul_lazy(y, a, b); s ^= y[0][0][0]; }
    *lazy = (ns_() - t0) / (double)n;
    (void)s;
}

double bench_fp_mul(const fp_t a, const fp_t b, long n)
{
    fp_t r; volatile limb_t s = 0;
    fp_mul(r, a, b);
    double t0 = ns_();
    for (long i = 0; i < n; i++) { fp_mul(r, a, b); s ^= r[0]; }
    (void)s;
    return (ns_() - t0) / (double)n;
}
