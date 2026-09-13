/*
 * Does an AVX-512 IFMA Montgomery multiply beat the scalar one? As a drop-in,
 * no, by a factor of four. Batched four at a time with the radix conversion
 * removed, it edges ahead. This is the measurement that says so, kept so step
 * 3 of issue #36 starts from numbers rather than an argument.
 *
 * NOT A BACKEND. Nothing in the library calls it. The drop-in failed the §6
 * gate outright, and under that gate a routine that does not win is deleted
 * rather than carried. What is carried is the probe.
 *
 * WHAT IT FOUND, this machine, -O3, separate processes (issue #40 has the
 * write-up):
 *
 *   fp_mul, mulx/adcx/adox, what ships         97.6 cycles
 *   IFMA drop-in: one multiply + conversion   389.7          4.0x slower
 *
 *   reduction loop, 1 multiply                303.4
 *   reduction loop, 2 in flight               152.5   per multiply
 *   reduction loop, 4 in flight                76.2   per multiply
 *   COMPLETE multiply, 4 in flight, no conv    82.3   per multiply   1.18x
 *
 * THE LOOP IS LATENCY-BOUND, which is the finding that matters. Every CIOS
 * step depends on the previous accumulator, so one multiplication cannot fill
 * the pipeline: two in flight halves the per-multiply cost and four quarters
 * it. Batching is therefore not an optimisation on top of an IFMA backend, it
 * is the thing that makes one possible at all.
 *
 * And even then it is close. 82.3 against 97.6 is 1.18x, which does not clear
 * the §6 bar of 1.2x, for a rewrite of the entire extension-field arithmetic
 * into a redundant representation.
 *
 * READ THE ABSOLUTE VECTOR NUMBERS WITH SUSPICION. The host this ran on does
 * NOT advertise AVX512_IFMA -- CPUID leaf 7.0 EBX bit 21 is clear and
 * /proc/cpuinfo does not list it -- yet vpmadd52luq executes and returns
 * correct results, so a hypervisor is masking the bit. Whether this silicon
 * runs IFMA at full speed is unknown. Hardware that admits to the feature may
 * do much better, and step 3 should be measured there before anything is
 * decided. That is now the first task, not a footnote.
 *
 * ONE NUMBER THAT DID NOT SETTLE. The conversion row reads 77 to 79 cycles
 * here and 37 to 40 in a second harness running the same kernel at the same
 * flags, both stable across runs, and the difference was not explained. It is
 * reported rather than quietly dropped, and it is not load-bearing: the
 * drop-in already loses four to one on the reduction loop alone, whichever
 * conversion figure is right.
 *
 * ONE CORRECTNESS SUBTLETY, easy to get wrong while still returning plausible
 * numbers. The library's Montgomery domain is R = 2^384; eight 52-bit CIOS
 * steps divide by 2^416. The missing 2^32 is folded into the radix conversion
 * of one operand, a bit offset in a repacking loop rather than a multiply:
 *
 *     CIOS(a << 32, b) = a*2^32*b*2^-416 = a*b*2^-384
 *
 * a<<32 < 2^413 < 2^416 and b < p, so the result stays below 2p and one
 * conditional subtraction still suffices. Checked against GMP over 300,000
 * random and boundary inputs before any timing was taken.
 *
 * HOW THIS FILE LEARNED TO MEASURE, since three separate readings were fiction
 * before they were facts. gcc infers that these kernels are pure, and of an
 * unchanged input a pure function need only be called once: it hoisted the
 * whole call out of the timing loop and reported 6 cycles for a Montgomery
 * multiply. neither noinline nor a memory clobber stopped it; varying the
 * operand each iteration did. Separately, timing all five cases through one
 * switch inside one loop cost the cheapest case the most. And the first
 * decomposition was taken at -O2 while the project builds at -O3, which made
 * the conversion look three times slower and four-way batching look two times
 * slower than they are.
 *
 *   cmake --build build --target ifma_probe && ./build/ifma_probe
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <cpuid.h>
#include "ifma_kernels.h"

static inline uint64_t rdtscp_now(void)
{
    unsigned lo, hi;
    __asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi) :: "rcx");
    return ((uint64_t)hi << 32) | lo;
}

static int have_ifma(void)
{
    unsigned a, b, c, d;
    if (__get_cpuid_max(0, NULL) < 7) return 0;
    if (!__get_cpuid_count(7, 0, &a, &b, &c, &d)) return 0;
    return (b & (1u << 21)) != 0;           /* leaf 7.0 EBX bit 21 */
}

/* Not const: the operand is varied each iteration so the compiler cannot treat
 * the call as loop-invariant. It only touches the low 16 bits, so the value
 * stays far below p and the CIOS bound holds. */
static uint64_t OPA[6] = {
    0x3f7a1c9e05b2d481ull, 0x91c4e7a2b8d6035full, 0x2e85b0c7d1936a4full,
    0x7c3d92e8a45f01b6ull, 0x18b6f34c9d27e50aull, 0x0a4e17d3b8c95f26ull };

static volatile uint64_t sink;

#define REPS  21
#define INNER 300000

int main(int argc, char **argv)
{
    int force = (argc > 1 && strcmp(argv[1], "--force") == 0);

    /* CPUID can lie, and on a virtual machine it lies downward. A dispatching
     * backend must believe it, because it cannot risk SIGILL. A probe must be
     * able to disbelieve it, or it cannot measure on the machines that are
     * interesting. */
    if (!have_ifma() && !force) {
        printf("CPUID reports no AVX512_IFMA.\n"
               "On a VM the bit may be masked while the instruction works.\n"
               "Re-run with --force to try it; a SIGILL means it really is absent.\n");
        return 0;
    }
    if (force && !have_ifma())
        printf("(CPUID says no IFMA; forced, so the instruction is tried anyway.\n"
               " Treat the absolute vector timings as suspect -- see the header.)\n\n");

    uint64_t a52[8], b52[8], base = OPA[0], a52_0;
    part_to52(a52, OPA, 32);
    part_to52(b52, OPA, 0);
    a52_0 = a52[0];

    printf("AVX-512 IFMA probe, BLS12-381, best of %d x %d\n\n", REPS, INNER);

    /* One tight loop per case. Putting all five behind a switch inside a
     * single loop looks tidier and is wrong: the branch and the other four
     * kernels' code stay in the loop body, which costs the cheapest case the
     * most. Measured that way the conversion read 78.8 cycles; in its own loop
     * it reads 37.8. The expensive cases barely moved, which is exactly how
     * that error hides. */
#define TIME(label, per, expr)                                                \
    do {                                                                      \
        uint64_t best = ~(uint64_t)0;                                         \
        for (int t = 0; t < REPS; t++) {                                      \
            uint64_t s0 = rdtscp_now();                                       \
            for (int i = 0; i < INNER; i++) {                                 \
                /* The operand must change, or gcc hoists these pure calls    \
                 * clean out of the loop and the probe reports single-digit   \
                 * cycles for a Montgomery multiply. Only the low 16 bits     \
                 * move, so the value stays far below p. */                   \
                OPA[0] = base ^ (uint64_t)(i & 0xffff);                       \
                a52[0] = a52_0 ^ (uint64_t)(i & 0xffff);                      \
                sink += (expr);                                               \
            }                                                                 \
            uint64_t e = rdtscp_now() - s0;                                    \
            if (e < best) best = e;                                           \
        }                                                                     \
        printf("  %-32s %8.2f cycles%s\n", label,                             \
               (double)best / INNER / (per), (per) > 1.0 ? "  per multiply" : ""); \
    } while (0)

    TIME("radix conversion alone",         1.0, part_convert(OPA));
    TIME("reduction loop, 1 multiply",     1.0, part_loop(a52, b52));
    TIME("reduction loop, 2 in flight",    2.0, part_loop_x2(a52, b52));
    TIME("reduction loop, 4 in flight",    4.0, part_loop_x4(a52, b52));
    TIME("COMPLETE multiply, 4 in flight", 4.0, part_full_x4(a52, b52));

    printf("\nsink 0x%016llx (so none of the above is optimised away)\n",
           (unsigned long long)sink);
    printf("\nCompare against the shipping fp_mul with\n"
           "  ./build/bench_BLS12_381 --reps 21\n"
           "in a SEPARATE process: AVX-512 pulls the core clock down, and timing\n"
           "both in one process flatters the vector code by slowing its rival.\n");
    return 0;
}
