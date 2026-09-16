/* Does lazy reduction pay in this library? Measured, twice over.
 *
 * Run it:
 *   gcc -O3 -std=c11 -Iinclude -DELIPS_CURVE_BLS12_381 \
 *       bench/lazy_probe.c bench/lazy_kernels.c src/arith/fp_lazy_x86_64.S \
 *       -o lazy_probe build/libelips_arith_BLS12_381.a -lm && ./lazy_probe
 *
 * The answer is in PAIRING-PERFORMANCE.md. In short: the algorithm is worth
 * up to 50% at Fp6, and none of it is reachable from C, because the 768-bit
 * intermediates spill and the spill costs more than twice the reductions
 * saved. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elips/fp.h"
#include "elips/fpx.h"
#include "elips/random.h"

void fpw_init(void);
void fp6_mul_lazy(fp6_t r, const fp6_t a, const fp6_t b);
void fp6_components(const fp6_t a, const fp6_t b, long n,
                    double *addm, double *mulw3, double *fp2mulw, double *redc2);
void fp6_compare(const fp6_t a, const fp6_t b, long n, double *eager, double *lazy);
double bench_fp_mul(const fp_t a, const fp_t b, long n);

static void rf(fp_t x){ limb_t k[ELIPS_ORDER_LIMBS]; elips_random_scalar(k); fp_from_limbs(x,k); }

int main(int argc, char **argv)
{
    long n = (argc > 1) ? atol(argv[1]) : 20000;
    fp6_t a, b, want, got;
    long bad = 0;

    fpw_init();

    /* Correctness first. A faster wrong kernel is not a result. */
    for (long i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) { rf(a[j][0]); rf(a[j][1]); rf(b[j][0]); rf(b[j][1]); }
        fp6_mul(want, a, b);
        fp6_mul_lazy(got, a, b);
        if (memcmp(want, got, sizeof(fp6_t))) { bad++; if (bad < 3) printf("FAIL at %ld\n", i); }
    }
    printf("%s  %ld differential rounds, %ld mismatches\n\n",
           bad ? "FAILED" : "PASS", n, bad);
    if (bad) return 1;

    double addm, mulw3, fp2mulw, redc2, eager, lazy;
    fp6_components(a, b, 3000000, &addm, &mulw3, &fp2mulw, &redc2);
    double fused = bench_fp_mul(a[0][0], b[0][0], 3000000);
    fp6_compare(a, b, 2000000, &eager, &lazy);

    printf("fused fp_mul (asm, mul+reduce) %7.2f ns   100%%\n", fused);
    printf("mulw  raw product              %7.2f ns   %3.0f%%\n", mulw3/3, 100*(mulw3/3)/fused);
    printf("redcw 768 -> 384 reduction     %7.2f ns   %3.0f%%\n", redc2/2, 100*(redc2/2)/fused);
    printf("  the two halves sum to %.0f%% of the fused multiply they replace\n\n",
           100*((mulw3/3)+(redc2/2))/fused);
    printf("fpw_addm 768-bit modular add   %7.2f ns\n", addm);
    printf("fp2_mulw 3 mulw + 3 fpw_subm   %7.2f ns  (parts sum to %.2f)\n\n",
           fp2mulw, mulw3 + 3*addm);

    double ideal = 18*(mulw3/3) + 6*(redc2/2) + 40*addm;
    printf("fp6_mul eager                  %7.1f ns\n", eager);
    printf("fp6_mul lazy, ideal            %7.1f ns   %+.1f%%\n", ideal, 100*(ideal-eager)/eager);
    printf("fp6_mul lazy, measured in C    %7.1f ns   %+.1f%%\n", lazy, 100*(lazy-eager)/eager);
    printf("\noverhead above ideal           %7.1f ns, more than the whole eager multiply\n",
           lazy - ideal);
    return 0;
}
