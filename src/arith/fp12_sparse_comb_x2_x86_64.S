/*
 * The sparse line multiply's combination at 768 bits, in one routine.
 *
 *   void elips_fp12_sparse_comb_x2_6_x86_64(limb_t out[12][12],
 *                                           const limb_t t0[6][12],
 *                                           const limb_t m[12][12],
 *                                           const limb_t s[6][12],
 *                                           const limb_t p[6]);
 *
 * Everything is component-major: a wide Fp2 value is two slots, real then
 * imaginary, and a wide Fp6 value is three such pairs.
 *
 *   t0 = f0 * (c0,0,0), three wide Fp2 products
 *   m  = the six products of f1 * (0,c3,c5):
 *        0,1 = a1*c5   2,3 = a2*c3   4,5 = a0*c3
 *        6,7 = a2*c5   8,9 = a0*c5  10,11 = a1*c3
 *   s  = (f0+f1) * (c0,c3,c5), already through elips_fp6_comb_x2
 *
 * This computes, with v^3 = xi and xi(x0,x1) = (x0 - x1, x0 + x1),
 *
 *   t1 = ( xi(a1c5 + a2c3), a0c3 + xi(a2c5), a0c5 + a1c3 )
 *   out[0..5]  = f0 = t0 + v*t1,  and v*t1 = (xi(t1_2), t1_0, t1_1)
 *   out[6..11] = f1 = s - t0 - t1
 *
 * ready for twelve redcw. Thirty wide operations, one prologue: the same
 * reason elips_fp6_comb_x2 exists, and see PAIRING-PERFORMANCE.md for what the
 * same thirty cost when written in C with adc/sbb intrinsics.
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

/* scratch: t1's six wide coefficients, then two working slots */
#define TT0     0
#define TT1    96
#define TT2   192
#define TT3   288
#define TT4   384
#define TT5   480
#define WK0   576
#define WK1   672
#define SC    768
#define PM    816
#define FRAME 864

.text

#include "arith/wide_x86_64.h"

/* m offsets */
#define M1_0     0
#define M1_1    96
#define M2_0   192
#define M2_1   288
#define M3_0   384
#define M3_1   480
#define M4_0   576
#define M4_1   672
#define M5_0   768
#define M5_1   864
#define M6_0   960
#define M6_1  1056
/* t0 and s offsets */
#define C0_0     0
#define C0_1    96
#define C1_0   192
#define C1_1   288
#define C2_0   384
#define C2_1   480
/* out offsets: f0 then f1 */
#define F0_0     0
#define F0_1    96
#define F0_2   192
#define F0_3   288
#define F0_4   384
#define F0_5   480
#define F1_0   576
#define F1_1   672
#define F1_2   768
#define F1_3   864
#define F1_4   960
#define F1_5  1056

    .globl SYM(elips_fp12_sparse_comb_x2_6_x86_64)
    FUNC_HIDDEN(elips_fp12_sparse_comb_x2_6_x86_64)
    FUNC_TYPE(elips_fp12_sparse_comb_x2_6_x86_64)
    .p2align 4
SYM(elips_fp12_sparse_comb_x2_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $FRAME, %rsp

    movq    %rcx, %rbx                  /* s */
    movq    %r8,  %rcx                  /* p */
    /* rdi = out, rsi = t0, rdx = m */

    /* ---- t1_0 = xi(a1c5 + a2c3) ---- */
    WADD %rsp, WK0, %rdx, M1_0, %rdx, M2_0
    WADD %rsp, WK1, %rdx, M1_1, %rdx, M2_1
    WSUB %rsp, TT0, %rsp, WK0,  %rsp, WK1
    WADD %rsp, TT1, %rsp, WK0,  %rsp, WK1

    /* ---- t1_1 = a0c3 + xi(a2c5) ---- */
    WSUB %rsp, WK0, %rdx, M4_0, %rdx, M4_1
    WADD %rsp, WK1, %rdx, M4_0, %rdx, M4_1
    WADD %rsp, TT2, %rdx, M3_0, %rsp, WK0
    WADD %rsp, TT3, %rdx, M3_1, %rsp, WK1

    /* ---- t1_2 = a0c5 + a1c3 ---- */
    WADD %rsp, TT4, %rdx, M5_0, %rdx, M6_0
    WADD %rsp, TT5, %rdx, M5_1, %rdx, M6_1

    /* ---- f1 = s - t0 - t1 ---- */
    WSUB %rsp, WK0,   %rbx, C0_0, %rsi, C0_0
    WSUB %rdi, F1_0,  %rsp, WK0,  %rsp, TT0
    WSUB %rsp, WK0,   %rbx, C0_1, %rsi, C0_1
    WSUB %rdi, F1_1,  %rsp, WK0,  %rsp, TT1
    WSUB %rsp, WK0,   %rbx, C1_0, %rsi, C1_0
    WSUB %rdi, F1_2,  %rsp, WK0,  %rsp, TT2
    WSUB %rsp, WK0,   %rbx, C1_1, %rsi, C1_1
    WSUB %rdi, F1_3,  %rsp, WK0,  %rsp, TT3
    WSUB %rsp, WK0,   %rbx, C2_0, %rsi, C2_0
    WSUB %rdi, F1_4,  %rsp, WK0,  %rsp, TT4
    WSUB %rsp, WK0,   %rbx, C2_1, %rsi, C2_1
    WSUB %rdi, F1_5,  %rsp, WK0,  %rsp, TT5

    /* ---- f0 = t0 + v*t1,  v*t1 = (xi(t1_2), t1_0, t1_1) ---- */
    WSUB %rsp, WK0,  %rsp, TT4,  %rsp, TT5
    WADD %rsp, WK1,  %rsp, TT4,  %rsp, TT5
    WADD %rdi, F0_0, %rsi, C0_0, %rsp, WK0
    WADD %rdi, F0_1, %rsi, C0_1, %rsp, WK1
    WADD %rdi, F0_2, %rsi, C1_0, %rsp, TT0
    WADD %rdi, F0_3, %rsi, C1_1, %rsp, TT1
    WADD %rdi, F0_4, %rsi, C2_0, %rsp, TT2
    WADD %rdi, F0_5, %rsi, C2_1, %rsp, TT3

    addq    $FRAME, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp12_sparse_comb_x2_6_x86_64)

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif
#endif
