/* MEASUREMENT ARTIFACT, not shipped. See PAIRING-PERFORMANCE.md, "The assembly
 * path, as far as it got".
 *
 * A complete lazy Fp6 multiply for the real library: six wide Fp2 products from
 * bench/fp2_mulx2_x86_64.S, combined at 768 bits, reduced once per output
 * coefficient instead of once per multiplication. Drop this over fp6_mul in
 * src/arith/fpx.c, add both .S files to CMakeLists.txt, build with
 * -DELIPS_LAZY_FP6, and the full test suite passes.
 *
 * It measures at PARITY with the eager tower, not better. What is missing is
 * the combination below in assembly; it is about 2,211 instructions per call
 * against 4,873 for the six multiplies it wraps. See the document.
 */
#if defined(ELIPS_LAZY_FP6) && (FP_LIMBS == 6) && defined(__x86_64__) && !defined(__ILP32__)
/* EXPERIMENT: Fp6 multiply with the reductions amortised 3:1.
 * Six wide Fp2 products, combined at 768 bits, reduced once per output
 * coefficient. Algebra identical to bench/lazy_kernels.c, which is
 * differential-checked against the eager form. */
void elips_fp_mulw_6_x86_64(limb_t w[12], const limb_t a[6], const limb_t b[6]);
void elips_fp_redcw_6_x86_64(limb_t r[6], const limb_t w[12], const limb_t p[6], limb_t n0);

typedef limb_t fpw_t[12];

/* p*R has six zero low limbs, so the conditional only touches the high half. */
/* adc/sbb directly. The __int128 borrow-extraction these replaced compiled to
 * roughly 5,000 instructions per fp6_mul for the combination alone, more than
 * the six wide multiplies it wraps. */
#include <x86intrin.h>

static inline void fpw_addm(fpw_t r, const fpw_t a, const fpw_t b)
{
    unsigned char c = 0, br = 0;
    limb_t hi[6], cand[6];
    for (int i = 0; i < 6; i++) c = _addcarry_u64(c, a[i], b[i], (unsigned long long *)&r[i]);
    for (int i = 0; i < 6; i++) c = _addcarry_u64(c, a[i+6], b[i+6], (unsigned long long *)&hi[i]);
    for (int i = 0; i < 6; i++) br = _subborrow_u64(br, hi[i], FP_MODULUS[i], (unsigned long long *)&cand[i]);
    limb_t keep = (limb_t)0 - (limb_t)(br & (c ^ 1));
    for (int i = 0; i < 6; i++) r[i+6] = (hi[i] & keep) | (cand[i] & ~keep);
}

static inline void fpw_subm(fpw_t r, const fpw_t a, const fpw_t b)
{
    unsigned char br = 0, c = 0;
    limb_t hi[6], cand[6];
    for (int i = 0; i < 6; i++) br = _subborrow_u64(br, a[i], b[i], (unsigned long long *)&r[i]);
    for (int i = 0; i < 6; i++) br = _subborrow_u64(br, a[i+6], b[i+6], (unsigned long long *)&hi[i]);
    for (int i = 0; i < 6; i++) c = _addcarry_u64(c, hi[i], FP_MODULUS[i], (unsigned long long *)&cand[i]);
    limb_t take = (limb_t)0 - (limb_t)br;
    for (int i = 0; i < 6; i++) r[i+6] = (cand[i] & take) | (hi[i] & ~take);
}

void elips_fp2_mulx2_6_x86_64(limb_t c0[12], limb_t c1[12],
                              const limb_t a[2][6], const limb_t b[2][6],
                              const limb_t p[6]);

/* One assembly routine, not three calls plus a C combination: 83 ns against
 * 128, and against 137 for the eager fp2_mul it replaces. */
static inline void fp2_mulw(fpw_t c0, fpw_t c1, const fp2_t a, const fp2_t b)
{
    elips_fp2_mulx2_6_x86_64(c0, c1, (const limb_t (*)[6])a,
                             (const limb_t (*)[6])b, FP_MODULUS);
}

void fp6_mul(fp6_t r, const fp6_t a, const fp6_t b)
{
    fpw_t v0[2], v1[2], v2[2], w01[2], w02[2], w12[2], t[2], u[2];
    fp2_t s1, s2;

    fp2_mulw(v0[0], v0[1], a[0], b[0]);
    fp2_mulw(v1[0], v1[1], a[1], b[1]);
    fp2_mulw(v2[0], v2[1], a[2], b[2]);
    fp2_add(s1, a[1], a[2]); fp2_add(s2, b[1], b[2]); fp2_mulw(w12[0], w12[1], s1, s2);
    fp2_add(s1, a[0], a[1]); fp2_add(s2, b[0], b[1]); fp2_mulw(w01[0], w01[1], s1, s2);
    fp2_add(s1, a[0], a[2]); fp2_add(s2, b[0], b[2]); fp2_mulw(w02[0], w02[1], s1, s2);

    /* xi = 1 + u : (x0, x1) -> (x0 - x1, x0 + x1) */
    for (int k = 0; k < 2; k++) { fpw_subm(t[k], w12[k], v1[k]); fpw_subm(t[k], t[k], v2[k]); }
    fpw_subm(u[0], t[0], t[1]); fpw_addm(u[1], t[0], t[1]);
    for (int k = 0; k < 2; k++) fpw_addm(t[k], u[k], v0[k]);
    elips_fp_redcw_6_x86_64(r[0][0], t[0], FP_MODULUS, FP_MONT_N0);
    elips_fp_redcw_6_x86_64(r[0][1], t[1], FP_MODULUS, FP_MONT_N0);

    for (int k = 0; k < 2; k++) { fpw_subm(t[k], w01[k], v0[k]); fpw_subm(t[k], t[k], v1[k]); }
    fpw_subm(u[0], v2[0], v2[1]); fpw_addm(u[1], v2[0], v2[1]);
    for (int k = 0; k < 2; k++) fpw_addm(t[k], t[k], u[k]);
    elips_fp_redcw_6_x86_64(r[1][0], t[0], FP_MODULUS, FP_MONT_N0);
    elips_fp_redcw_6_x86_64(r[1][1], t[1], FP_MODULUS, FP_MONT_N0);

    for (int k = 0; k < 2; k++) { fpw_subm(t[k], w02[k], v0[k]); fpw_subm(t[k], t[k], v2[k]);
                                  fpw_addm(t[k], t[k], v1[k]); }
    elips_fp_redcw_6_x86_64(r[2][0], t[0], FP_MODULUS, FP_MONT_N0);
    elips_fp_redcw_6_x86_64(r[2][1], t[1], FP_MODULUS, FP_MONT_N0);
}
