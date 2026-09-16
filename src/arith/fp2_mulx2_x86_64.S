/*
 * Fp2 multiply into double width: c0, c1 each 768 bits, unreduced.
 * blst calls its equivalent mulx_382x.
 *
 *   void elips_fp2_mulx2_6_x86_64(limb_t c0[12], limb_t c1[12],
 *                                 const limb_t a[2][6], const limb_t b[2][6],
 *                                 const limb_t p[6]);
 *
 * WHY. Three separate calls to the raw-product routine measure 137 ns for the
 * Karatsuba, the same as the eager fp2_mul they are meant to beat, although the
 * three products alone are 52 ns. The rest is C: each call saves and restores
 * five callee-saved registers, and each 96-byte product is stored and then
 * immediately reloaded by the combination. Inlining all three into one routine
 * removes the prologues and lets the stores and the loads be scheduled
 * together.
 *
 * Karatsuba, with the sums reduced so both operands stay below p:
 *     t0 = a0*b0,  t1 = a1*b1,  s = (a0+a1)(b0+b1)
 *     c0 = t0 - t1,  c1 = s - t0 - t1,   all mod p*R
 * p*R has six zero low limbs, so a subtraction mod p*R can only borrow out of
 * the high half: the low six limbs are a plain subtract and only the top six
 * need the conditional add of p. Same fact the C version uses.
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

.text

/* stack layout: t0 at 0, t1 at 96, s at 192, as at 288, bs at 336 */
#define T0   0
#define T1   96
#define S    192
#define AS   288
#define BS   336
#define STASH 384
#define PAD   392
#define FRAME 448

/* One round of the raw product: T += a*b[i], store the settled low word.
 * rsi = a, rbx = b, rdi = out. Identical to MULROUND6 in fp_lazy_x86_64.S. */
.macro PROD6 t0,t1,t2,t3,t4,t5,t6, off
    movq    \off(%rbx), %rdx
    xorl    %eax, %eax
    mulxq   0(%rsi), %rax, %r15
    adcxq   %rax, \t0
    adoxq   %r15, \t1
    mulxq   8(%rsi), %rax, %r15
    adcxq   %rax, \t1
    adoxq   %r15, \t2
    mulxq   16(%rsi), %rax, %r15
    adcxq   %rax, \t2
    adoxq   %r15, \t3
    mulxq   24(%rsi), %rax, %r15
    adcxq   %rax, \t3
    adoxq   %r15, \t4
    mulxq   32(%rsi), %rax, %r15
    adcxq   %rax, \t4
    adoxq   %r15, \t5
    mulxq   40(%rsi), %rax, %r15
    adcxq   %rax, \t5
    adoxq   %r15, \t6
    movl    $0, %eax
    adcxq   %rax, \t6
    adoxq   %rax, \t6
    movq    \t0, \off(%rdi)
    xorl    %eax, %eax
    movq    %rax, \t0
.endm

/* a full 384x384 -> 768 product, rsi=a rbx=b rdi=out */
.macro PRODUCT
    xorl    %r8d,  %r8d
    xorl    %r9d,  %r9d
    xorl    %r10d, %r10d
    xorl    %r11d, %r11d
    xorl    %r12d, %r12d
    xorl    %r13d, %r13d
    xorl    %r14d, %r14d
    PROD6  %r8,  %r9,  %r10, %r11, %r12, %r13, %r14,  0
    PROD6  %r9,  %r10, %r11, %r12, %r13, %r14, %r8,   8
    PROD6  %r10, %r11, %r12, %r13, %r14, %r8,  %r9,   16
    PROD6  %r11, %r12, %r13, %r14, %r8,  %r9,  %r10,  24
    PROD6  %r12, %r13, %r14, %r8,  %r9,  %r10, %r11,  32
    PROD6  %r13, %r14, %r8,  %r9,  %r10, %r11, %r12,  40
    movq    %r14, 48(%rdi)
    movq    %r8,  56(%rdi)
    movq    %r9,  64(%rdi)
    movq    %r10, 72(%rdi)
    movq    %r11, 80(%rdi)
    movq    %r12, 88(%rdi)
.endm

/* r = a - b mod p*R, 12 limbs. rcx = p. Clobbers rax, r8-r14. */
.macro SUBPR dst, srcA, srcB
    movq    \srcA+0(%rsp),  %r8
    movq    \srcA+8(%rsp),  %r9
    movq    \srcA+16(%rsp), %r10
    movq    \srcA+24(%rsp), %r11
    movq    \srcA+32(%rsp), %r12
    movq    \srcA+40(%rsp), %r13
    subq    \srcB+0(%rsp),  %r8
    sbbq    \srcB+8(%rsp),  %r9
    sbbq    \srcB+16(%rsp), %r10
    sbbq    \srcB+24(%rsp), %r11
    sbbq    \srcB+32(%rsp), %r12
    sbbq    \srcB+40(%rsp), %r13
    movq    %r8,  0(\dst)
    movq    %r9,  8(\dst)
    movq    %r10, 16(\dst)
    movq    %r11, 24(\dst)
    movq    %r12, 32(\dst)
    movq    %r13, 40(\dst)

    movq    \srcA+48(%rsp), %r8
    movq    \srcA+56(%rsp), %r9
    movq    \srcA+64(%rsp), %r10
    movq    \srcA+72(%rsp), %r11
    movq    \srcA+80(%rsp), %r12
    movq    \srcA+88(%rsp), %r13
    sbbq    \srcB+48(%rsp), %r8
    sbbq    \srcB+56(%rsp), %r9
    sbbq    \srcB+64(%rsp), %r10
    sbbq    \srcB+72(%rsp), %r11
    sbbq    \srcB+80(%rsp), %r12
    sbbq    \srcB+88(%rsp), %r13
    sbbq    %r14, %r14                  /* r14 = -(borrow) */

    /* Mask p into scratch FIRST. andq sets flags, so interleaving it with the
     * add chain below destroys the carry between adcq's -- which is exactly
     * the bug the first version of this routine had. */
    movq    0(%rcx),  %rax
    andq    %r14, %rax
    movq    %rax, PAD+0(%rsp)
    movq    8(%rcx),  %rax
    andq    %r14, %rax
    movq    %rax, PAD+8(%rsp)
    movq    16(%rcx), %rax
    andq    %r14, %rax
    movq    %rax, PAD+16(%rsp)
    movq    24(%rcx), %rax
    andq    %r14, %rax
    movq    %rax, PAD+24(%rsp)
    movq    32(%rcx), %rax
    andq    %r14, %rax
    movq    %rax, PAD+32(%rsp)
    movq    40(%rcx), %rax
    andq    %r14, %rax
    movq    %rax, PAD+40(%rsp)

    addq    PAD+0(%rsp),  %r8
    adcq    PAD+8(%rsp),  %r9
    adcq    PAD+16(%rsp), %r10
    adcq    PAD+24(%rsp), %r11
    adcq    PAD+32(%rsp), %r12
    adcq    PAD+40(%rsp), %r13

    movq    %r8,  48(\dst)
    movq    %r9,  56(\dst)
    movq    %r10, 64(\dst)
    movq    %r11, 72(\dst)
    movq    %r12, 80(\dst)
    movq    %r13, 88(\dst)
.endm


    .globl SYM(elips_fp2_mulx2_6_x86_64)
    FUNC_HIDDEN(elips_fp2_mulx2_6_x86_64)
    FUNC_TYPE(elips_fp2_mulx2_6_x86_64)
    .p2align 4
SYM(elips_fp2_mulx2_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $FRAME, %rsp

    movq    %rdi, %rbp                  /* keep c0 */
    movq    %rsi, %rax
    movq    %rax, STASH(%rsp)           /* stash c1 clear of BS */
    movq    %rdx, %rsi                  /* a */
    movq    %rcx, %rbx                  /* b */
    movq    %r8,  %rcx                  /* p */

    /* as = a0 + a1 mod p, bs = b0 + b1 mod p */
    movq    0(%rsi),  %r8
    movq    8(%rsi),  %r9
    movq    16(%rsi), %r10
    movq    24(%rsi), %r11
    movq    32(%rsi), %r12
    movq    40(%rsi), %r13
    addq    48(%rsi), %r8
    adcq    56(%rsi), %r9
    adcq    64(%rsi), %r10
    adcq    72(%rsi), %r11
    adcq    80(%rsi), %r12
    adcq    88(%rsi), %r13
    movq    %r8,  AS+0(%rsp)
    movq    %r9,  AS+8(%rsp)
    movq    %r10, AS+16(%rsp)
    movq    %r11, AS+24(%rsp)
    movq    %r12, AS+32(%rsp)
    movq    %r13, AS+40(%rsp)
    subq    0(%rcx),  %r8
    sbbq    8(%rcx),  %r9
    sbbq    16(%rcx), %r10
    sbbq    24(%rcx), %r11
    sbbq    32(%rcx), %r12
    sbbq    40(%rcx), %r13
    movq    AS+0(%rsp), %rax
    cmovncq %r8, %rax
    movq    %rax, AS+0(%rsp)
    movq    AS+8(%rsp), %rax
    cmovncq %r9, %rax
    movq    %rax, AS+8(%rsp)
    movq    AS+16(%rsp), %rax
    cmovncq %r10, %rax
    movq    %rax, AS+16(%rsp)
    movq    AS+24(%rsp), %rax
    cmovncq %r11, %rax
    movq    %rax, AS+24(%rsp)
    movq    AS+32(%rsp), %rax
    cmovncq %r12, %rax
    movq    %rax, AS+32(%rsp)
    movq    AS+40(%rsp), %rax
    cmovncq %r13, %rax
    movq    %rax, AS+40(%rsp)

    movq    0(%rbx),  %r8
    movq    8(%rbx),  %r9
    movq    16(%rbx), %r10
    movq    24(%rbx), %r11
    movq    32(%rbx), %r12
    movq    40(%rbx), %r13
    addq    48(%rbx), %r8
    adcq    56(%rbx), %r9
    adcq    64(%rbx), %r10
    adcq    72(%rbx), %r11
    adcq    80(%rbx), %r12
    adcq    88(%rbx), %r13
    movq    %r8,  BS+0(%rsp)
    movq    %r9,  BS+8(%rsp)
    movq    %r10, BS+16(%rsp)
    movq    %r11, BS+24(%rsp)
    movq    %r12, BS+32(%rsp)
    movq    %r13, BS+40(%rsp)
    subq    0(%rcx),  %r8
    sbbq    8(%rcx),  %r9
    sbbq    16(%rcx), %r10
    sbbq    24(%rcx), %r11
    sbbq    32(%rcx), %r12
    sbbq    40(%rcx), %r13
    movq    BS+0(%rsp), %rax
    cmovncq %r8, %rax
    movq    %rax, BS+0(%rsp)
    movq    BS+8(%rsp), %rax
    cmovncq %r9, %rax
    movq    %rax, BS+8(%rsp)
    movq    BS+16(%rsp), %rax
    cmovncq %r10, %rax
    movq    %rax, BS+16(%rsp)
    movq    BS+24(%rsp), %rax
    cmovncq %r11, %rax
    movq    %rax, BS+24(%rsp)
    movq    BS+32(%rsp), %rax
    cmovncq %r12, %rax
    movq    %rax, BS+32(%rsp)
    movq    BS+40(%rsp), %rax
    cmovncq %r13, %rax
    movq    %rax, BS+40(%rsp)

    /* t0 = a0*b0 */
    leaq    T0(%rsp), %rdi
    PRODUCT
    /* t1 = a1*b1 */
    addq    $48, %rsi
    addq    $48, %rbx
    leaq    T1(%rsp), %rdi
    PRODUCT
    /* s = as*bs */
    leaq    AS(%rsp), %rsi
    leaq    BS(%rsp), %rbx
    leaq    S(%rsp),  %rdi
    PRODUCT

    /* c0 = t0 - t1 ; then s = s - t0 ; c1 = s - t1 */
    SUBPR   %rbp, T0, T1                /* c0 = t0 - t1 */
    leaq    S(%rsp), %r15
    SUBPR   %r15, S, T0                 /* s -= t0, in place */
    movq    STASH(%rsp), %rdi
    SUBPR   %rdi, S, T1                 /* c1 = s - t1 */

    addq    $FRAME, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp2_mulx2_6_x86_64)

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif
#endif
