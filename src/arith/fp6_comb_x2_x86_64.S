/*
 * The Fp6 Karatsuba combination at 768 bits, in one routine.
 *
 *   void elips_fp6_comb_x2_6_x86_64(limb_t out[6][12],
 *                                   const limb_t prod[12][12],
 *                                   const limb_t p[6]);
 *
 * prod holds the six wide Fp2 products from elips_fp2_mulx2, component-major:
 *   0,1 = v0   2,3 = v1   4,5 = v2   6,7 = w01   8,9 = w02   10,11 = w12
 * out receives the three wide Fp2 coefficients, ready for redcw.
 *
 *   c0 = v0 + xi(w12 - v1 - v2)
 *   c1 = (w01 - v0 - v1) + xi(v2)
 *   c2 = (w02 - v0 - v2) + v1        with xi(x0,x1) = (x0 - x1, x0 + x1)
 *
 * WHY THIS EXISTS. With the wide multiply already in assembly, this
 * combination was still 2,211 instructions per fp6_mul against 4,873 for the
 * six multiplies it wraps -- and that is AFTER replacing __int128 borrow
 * extraction with adc/sbb intrinsics, which had it at over 5,000. Twenty-two
 * calls' worth of prologue and twelve-limb round trips, done once here.
 *
 * The wide add and subtract it is built from are in arith/wide_x86_64.h,
 * shared with the sparse combination.
 */
#if defined(__x86_64__) && !defined(__ILP32__)

#if defined(__APPLE__)
#  define SYM(x)         _##x
#  define FUNC_HIDDEN(x)
#  define FUNC_TYPE(x)
#  define FUNC_SIZE(x)
#else
#  define SYM(x)         x
#  define FUNC_HIDDEN(x) .hidden x
#  define FUNC_TYPE(x)   .type x,@function
#  define FUNC_SIZE(x)   .size x, .-x
#endif

/* scratch */
#define T0 0
#define T1 96
#define U0 192
#define U1 288
#define SC 384
#define PM 432
#define FRAME 480

.text

#include "arith/wide_x86_64.h"

/* prod offsets, component-major */
#define V0_0    0
#define V0_1   96
#define V1_0  192
#define V1_1  288
#define V2_0  384
#define V2_1  480
#define W01_0 576
#define W01_1 672
#define W02_0 768
#define W02_1 864
#define W12_0 960
#define W12_1 1056
/* out offsets */
#define O0_0    0
#define O0_1   96
#define O1_0  192
#define O1_1  288
#define O2_0  384
#define O2_1  480

    .globl SYM(elips_fp6_comb_x2_6_x86_64)
    FUNC_HIDDEN(elips_fp6_comb_x2_6_x86_64)
    FUNC_TYPE(elips_fp6_comb_x2_6_x86_64)
    .p2align 4
SYM(elips_fp6_comb_x2_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $FRAME, %rsp

    movq    %rdx, %rcx                  /* p */
    /* rdi = out, rsi = prod */

    /* ---- c0 = v0 + xi(w12 - v1 - v2) ---- */
    WSUB %rsp, T0, %rsi, W12_0, %rsi, V1_0
    WSUB %rsp, T0, %rsp, T0,    %rsi, V2_0
    WSUB %rsp, T1, %rsi, W12_1, %rsi, V1_1
    WSUB %rsp, T1, %rsp, T1,    %rsi, V2_1
    WSUB %rsp, U0, %rsp, T0,    %rsp, T1
    WADD %rsp, U1, %rsp, T0,    %rsp, T1
    WADD %rdi, O0_0, %rsp, U0,  %rsi, V0_0
    WADD %rdi, O0_1, %rsp, U1,  %rsi, V0_1

    /* ---- c1 = (w01 - v0 - v1) + xi(v2) ---- */
    WSUB %rsp, T0, %rsi, W01_0, %rsi, V0_0
    WSUB %rsp, T0, %rsp, T0,    %rsi, V1_0
    WSUB %rsp, T1, %rsi, W01_1, %rsi, V0_1
    WSUB %rsp, T1, %rsp, T1,    %rsi, V1_1
    WSUB %rsp, U0, %rsi, V2_0,  %rsi, V2_1
    WADD %rsp, U1, %rsi, V2_0,  %rsi, V2_1
    WADD %rdi, O1_0, %rsp, T0,  %rsp, U0
    WADD %rdi, O1_1, %rsp, T1,  %rsp, U1

    /* ---- c2 = (w02 - v0 - v2) + v1 ---- */
    WSUB %rsp, T0, %rsi, W02_0, %rsi, V0_0
    WSUB %rsp, T0, %rsp, T0,    %rsi, V2_0
    WSUB %rsp, T1, %rsi, W02_1, %rsi, V0_1
    WSUB %rsp, T1, %rsp, T1,    %rsi, V2_1
    WADD %rdi, O2_0, %rsp, T0,  %rsi, V1_0
    WADD %rdi, O2_1, %rsp, T1,  %rsi, V1_1

    addq    $FRAME, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp6_comb_x2_6_x86_64)

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif
#endif
