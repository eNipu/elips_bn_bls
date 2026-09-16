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
 * Everything stays mod p*R, and p*R has six zero low limbs, so a borrow can
 * never leave the low half: the low six limbs are a plain subtract and only the
 * top six take the conditional add of p. Same fact the C used.
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

/* dst = a - b  mod p*R.  rcx = p.  Clobbers rax, r8-r13. */
.macro WSUB dbase, doff, abase, aoff, bbase, boff
    movq    \aoff+0(\abase),  %r8
    movq    \aoff+8(\abase),  %r9
    movq    \aoff+16(\abase), %r10
    movq    \aoff+24(\abase), %r11
    movq    \aoff+32(\abase), %r12
    movq    \aoff+40(\abase), %r13
    subq    \boff+0(\bbase),  %r8
    sbbq    \boff+8(\bbase),  %r9
    sbbq    \boff+16(\bbase), %r10
    sbbq    \boff+24(\bbase), %r11
    sbbq    \boff+32(\bbase), %r12
    sbbq    \boff+40(\bbase), %r13
    movq    %r8,  \doff+0(\dbase)
    movq    %r9,  \doff+8(\dbase)
    movq    %r10, \doff+16(\dbase)
    movq    %r11, \doff+24(\dbase)
    movq    %r12, \doff+32(\dbase)
    movq    %r13, \doff+40(\dbase)
    movq    \aoff+48(\abase), %r8
    movq    \aoff+56(\abase), %r9
    movq    \aoff+64(\abase), %r10
    movq    \aoff+72(\abase), %r11
    movq    \aoff+80(\abase), %r12
    movq    \aoff+88(\abase), %r13
    sbbq    \boff+48(\bbase), %r8
    sbbq    \boff+56(\bbase), %r9
    sbbq    \boff+64(\bbase), %r10
    sbbq    \boff+72(\bbase), %r11
    sbbq    \boff+80(\bbase), %r12
    sbbq    \boff+88(\bbase), %r13
    sbbq    %rax, %rax                  /* rax = -(borrow) */
    /* mask p into scratch first: andq sets flags and would break the chain */
    MASKP
    addq    PM+0(%rsp),  %r8
    adcq    PM+8(%rsp),  %r9
    adcq    PM+16(%rsp), %r10
    adcq    PM+24(%rsp), %r11
    adcq    PM+32(%rsp), %r12
    adcq    PM+40(%rsp), %r13
    movq    %r8,  \doff+48(\dbase)
    movq    %r9,  \doff+56(\dbase)
    movq    %r10, \doff+64(\dbase)
    movq    %r11, \doff+72(\dbase)
    movq    %r12, \doff+80(\dbase)
    movq    %r13, \doff+88(\dbase)
.endm

/* dst = a + b  mod p*R. */
.macro WADD dbase, doff, abase, aoff, bbase, boff
    movq    \aoff+0(\abase),  %r8
    movq    \aoff+8(\abase),  %r9
    movq    \aoff+16(\abase), %r10
    movq    \aoff+24(\abase), %r11
    movq    \aoff+32(\abase), %r12
    movq    \aoff+40(\abase), %r13
    addq    \boff+0(\bbase),  %r8
    adcq    \boff+8(\bbase),  %r9
    adcq    \boff+16(\bbase), %r10
    adcq    \boff+24(\bbase), %r11
    adcq    \boff+32(\bbase), %r12
    adcq    \boff+40(\bbase), %r13
    movq    %r8,  \doff+0(\dbase)
    movq    %r9,  \doff+8(\dbase)
    movq    %r10, \doff+16(\dbase)
    movq    %r11, \doff+24(\dbase)
    movq    %r12, \doff+32(\dbase)
    movq    %r13, \doff+40(\dbase)
    movq    \aoff+48(\abase), %r8
    movq    \aoff+56(\abase), %r9
    movq    \aoff+64(\abase), %r10
    movq    \aoff+72(\abase), %r11
    movq    \aoff+80(\abase), %r12
    movq    \aoff+88(\abase), %r13
    adcq    \boff+48(\bbase), %r8
    adcq    \boff+56(\bbase), %r9
    adcq    \boff+64(\bbase), %r10
    adcq    \boff+72(\bbase), %r11
    adcq    \boff+80(\bbase), %r12
    adcq    \boff+88(\bbase), %r13
    sbbq    %r14, %r14                  /* r14 = -(carry out of the top) */
    /* subtract p when the sum reached p*R: either it carried out, or it did
     * not borrow when p is taken off the high half. */
    movq    %r8,  SC+0(%rsp)
    movq    %r9,  SC+8(%rsp)
    movq    %r10, SC+16(%rsp)
    movq    %r11, SC+24(%rsp)
    movq    %r12, SC+32(%rsp)
    movq    %r13, SC+40(%rsp)
    subq    0(%rcx),  %r8
    sbbq    8(%rcx),  %r9
    sbbq    16(%rcx), %r10
    sbbq    24(%rcx), %r11
    sbbq    32(%rcx), %r12
    sbbq    40(%rcx), %r13
    sbbq    $0, %r14                    /* borrow against the carry */
    cmovcq  SC+0(%rsp),  %r8
    cmovcq  SC+8(%rsp),  %r9
    cmovcq  SC+16(%rsp), %r10
    cmovcq  SC+24(%rsp), %r11
    cmovcq  SC+32(%rsp), %r12
    cmovcq  SC+40(%rsp), %r13
    movq    %r8,  \doff+48(\dbase)
    movq    %r9,  \doff+56(\dbase)
    movq    %r10, \doff+64(\dbase)
    movq    %r11, \doff+72(\dbase)
    movq    %r12, \doff+80(\dbase)
    movq    %r13, \doff+88(\dbase)
.endm

.macro MASKP
    movq    0(%rcx),  %r14
    andq    %rax, %r14
    movq    %r14, PM+0(%rsp)
    movq    8(%rcx),  %r14
    andq    %rax, %r14
    movq    %r14, PM+8(%rsp)
    movq    16(%rcx), %r14
    andq    %rax, %r14
    movq    %r14, PM+16(%rsp)
    movq    24(%rcx), %r14
    andq    %rax, %r14
    movq    %r14, PM+24(%rsp)
    movq    32(%rcx), %r14
    andq    %rax, %r14
    movq    %r14, PM+32(%rsp)
    movq    40(%rcx), %r14
    andq    %rax, %r14
    movq    %r14, PM+40(%rsp)
.endm

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
