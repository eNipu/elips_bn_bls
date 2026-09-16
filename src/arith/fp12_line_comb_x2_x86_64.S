/*
 * The rotated line multiply's combination at 768 bits, in one routine.
 *
 *   void elips_fp12_line_comb_x2_6_x86_64(limb_t out[12][12],
 *                                         const limb_t prod[26][12],
 *                                         const limb_t p[6]);
 *
 * Everything is component-major: a wide Fp2 value is two slots, real then
 * imaginary. prod holds the thirteen wide Fp2 products of
 * fp12_mul_sparse035w, in three blocks.
 *
 *   t0 = f0 * (A, 0, B)          0,1 = P1   2,3 = P2   4,5 = P3
 *                                6,7 = P4   8,9 = P5
 *   s  = (f0+f1) * (E, 0, B)   10,11 = R1 12,13 = R2 14,15 = R3
 *                              16,17 = R4 18,19 = R5
 *   t1 = f1 * (C, 0, 0)        20,21 = Q0 22,23 = Q1 24,25 = Q2
 *
 * With v^3 = xi and xi(x0,x1) = (x0 - x1, x0 + x1), the two five-product
 * blocks combine as
 *
 *   c0 = P1 + xi(P2),  c1 = P3 + xi(P4),  c2 = P5 - P1 - P4
 *
 * and then
 *
 *   out[0..5]  = f0 = t0 + v*t1,  v*t1 = (xi(t1_2), t1_0, t1_1)
 *   out[6..11] = f1 = s - t0 - t1
 *
 * ready for twelve redcw. Forty-four wide operations under one prologue,
 * against the fifty-two the unrotated form needed across two routines.
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

/* scratch: t0's six wide coefficients, s's six, then two working slots */
#define TA0      0
#define TA1     96
#define TA2    192
#define TA3    288
#define TA4    384
#define TA5    480
#define SB0    576
#define SB1    672
#define SB2    768
#define SB3    864
#define SB4    960
#define SB5   1056
#define WK0   1152
#define WK1   1248
#define SC    1344
#define PM    1392
#define FRAME 1440

.text

#include "arith/wide_x86_64.h"

/* prod offsets */
#define P1_0     0
#define P1_1    96
#define P2_0   192
#define P2_1   288
#define P3_0   384
#define P3_1   480
#define P4_0   576
#define P4_1   672
#define P5_0   768
#define P5_1   864
#define R1_0   960
#define R1_1  1056
#define R2_0  1152
#define R2_1  1248
#define R3_0  1344
#define R3_1  1440
#define R4_0  1536
#define R4_1  1632
#define R5_0  1728
#define R5_1  1824
#define Q0_0  1920
#define Q0_1  2016
#define Q1_0  2112
#define Q1_1  2208
#define Q2_0  2304
#define Q2_1  2400
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

    .globl SYM(elips_fp12_line_comb_x2_6_x86_64)
    FUNC_HIDDEN(elips_fp12_line_comb_x2_6_x86_64)
    FUNC_TYPE(elips_fp12_line_comb_x2_6_x86_64)
    .p2align 4
SYM(elips_fp12_line_comb_x2_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $FRAME, %rsp

    movq    %rdx, %rcx                  /* p */
    /* rdi = out, rsi = prod */

    /* ---- t0_0 = P1 + xi(P2) ---- */
    WSUB %rsp, WK0, %rsi, P2_0, %rsi, P2_1
    WADD %rsp, WK1, %rsi, P2_0, %rsi, P2_1
    WADD %rsp, TA0, %rsi, P1_0, %rsp, WK0
    WADD %rsp, TA1, %rsi, P1_1, %rsp, WK1
    /* ---- t0_1 = P3 + xi(P4) ---- */
    WSUB %rsp, WK0, %rsi, P4_0, %rsi, P4_1
    WADD %rsp, WK1, %rsi, P4_0, %rsi, P4_1
    WADD %rsp, TA2, %rsi, P3_0, %rsp, WK0
    WADD %rsp, TA3, %rsi, P3_1, %rsp, WK1
    /* ---- t0_2 = P5 - P1 - P4 ---- */
    WSUB %rsp, WK0, %rsi, P5_0, %rsi, P1_0
    WSUB %rsp, TA4, %rsp, WK0,  %rsi, P4_0
    WSUB %rsp, WK0, %rsi, P5_1, %rsi, P1_1
    WSUB %rsp, TA5, %rsp, WK0,  %rsi, P4_1

    /* ---- s_0 = R1 + xi(R2) ---- */
    WSUB %rsp, WK0, %rsi, R2_0, %rsi, R2_1
    WADD %rsp, WK1, %rsi, R2_0, %rsi, R2_1
    WADD %rsp, SB0, %rsi, R1_0, %rsp, WK0
    WADD %rsp, SB1, %rsi, R1_1, %rsp, WK1
    /* ---- s_1 = R3 + xi(R4) ---- */
    WSUB %rsp, WK0, %rsi, R4_0, %rsi, R4_1
    WADD %rsp, WK1, %rsi, R4_0, %rsi, R4_1
    WADD %rsp, SB2, %rsi, R3_0, %rsp, WK0
    WADD %rsp, SB3, %rsi, R3_1, %rsp, WK1
    /* ---- s_2 = R5 - R1 - R4 ---- */
    WSUB %rsp, WK0, %rsi, R5_0, %rsi, R1_0
    WSUB %rsp, SB4, %rsp, WK0,  %rsi, R4_0
    WSUB %rsp, WK0, %rsi, R5_1, %rsi, R1_1
    WSUB %rsp, SB5, %rsp, WK0,  %rsi, R4_1

    /* ---- f1 = s - t0 - t1,  t1 = (Q0, Q1, Q2) ---- */
    WSUB %rsp, WK0,  %rsp, SB0, %rsp, TA0
    WSUB %rdi, F1_0, %rsp, WK0, %rsi, Q0_0
    WSUB %rsp, WK0,  %rsp, SB1, %rsp, TA1
    WSUB %rdi, F1_1, %rsp, WK0, %rsi, Q0_1
    WSUB %rsp, WK0,  %rsp, SB2, %rsp, TA2
    WSUB %rdi, F1_2, %rsp, WK0, %rsi, Q1_0
    WSUB %rsp, WK0,  %rsp, SB3, %rsp, TA3
    WSUB %rdi, F1_3, %rsp, WK0, %rsi, Q1_1
    WSUB %rsp, WK0,  %rsp, SB4, %rsp, TA4
    WSUB %rdi, F1_4, %rsp, WK0, %rsi, Q2_0
    WSUB %rsp, WK0,  %rsp, SB5, %rsp, TA5
    WSUB %rdi, F1_5, %rsp, WK0, %rsi, Q2_1

    /* ---- f0 = t0 + (xi(Q2), Q0, Q1) ---- */
    WSUB %rsp, WK0,  %rsi, Q2_0, %rsi, Q2_1
    WADD %rsp, WK1,  %rsi, Q2_0, %rsi, Q2_1
    WADD %rdi, F0_0, %rsp, TA0,  %rsp, WK0
    WADD %rdi, F0_1, %rsp, TA1,  %rsp, WK1
    WADD %rdi, F0_2, %rsp, TA2,  %rsi, Q0_0
    WADD %rdi, F0_3, %rsp, TA3,  %rsi, Q0_1
    WADD %rdi, F0_4, %rsp, TA4,  %rsi, Q1_0
    WADD %rdi, F0_5, %rsp, TA5,  %rsi, Q1_1

    addq    $FRAME, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp12_line_comb_x2_6_x86_64)

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif
#endif
