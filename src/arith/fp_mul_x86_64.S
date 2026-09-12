/*
 * Montgomery multiplication for x86-64 with BMI2 and ADX.
 *
 * Same CIOS algorithm as the portable C in src/arith/fp.c, so the two can be
 * compared word for word. What assembly adds is the thing C cannot express:
 * TWO independent carry chains. adcx uses the carry flag and adox uses the
 * overflow flag, so the low halves and the high halves of the partial products
 * accumulate at the same time instead of queueing behind one flag.
 *
 *   void elips_fp_mul_6_x86_64(limb_t r[6], const limb_t a[6],
 *                              const limb_t b[6], const limb_t p[6],
 *                              limb_t n0);
 *
 * System V AMD64: rdi=r, rsi=a, rdx=b, rcx=p, r8=n0.
 *
 * Two facts the code relies on, both proved rather than assumed:
 *
 *   1. a < p and T < 2p on entry to every round, so T stays below 2p and ONE
 *      conditional subtraction is enough at the end. This is the standard CIOS
 *      invariant and needs a < p, which every caller satisfies because
 *      Montgomery representatives are reduced.
 *
 *   2. The n+1 word accumulator never overflows, so the adox chain's carry out
 *      of the top word is always zero and is not added anywhere. This needs a
 *      SPARSE modulus, p < 2^(64n-1): 381 < 383 for BLS12-381, and 461 and 462
 *      are both < 511 for the 8 limb curves. src/arith/fp.c static-asserts it.
 *
 * Rounds rotate register names instead of moving data. t0 is exactly zero
 * after its round's reduction step, which is what the multiplier m is chosen
 * to arrange, so the word that leaves the bottom is already the zero the top
 * needs. The shift by one limb costs nothing.
 *
 * Constant time by construction: no branches at all, and the only conditional
 * is cmov, whose operands live at fixed addresses.
 */

#if defined(__x86_64__) && !defined(__ILP32__)

/* Symbol decoration, one directive per SOURCE line.
 *
 * It has to be one per line. A cpp macro cannot contain a newline, so the
 * obvious thing is to separate the directives with semicolons -- and that
 * silently produces an object file with no symbols in it on arm64 Darwin,
 * where ";" begins a comment. The .globl survives, everything after it on the
 * line does not, and the link fails with the function "referenced from" its
 * only caller and never defined. It assembles without a word of complaint.
 * So: each macro expands to at most one directive, and the call sites stack
 * them. */
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

/* ----------------------------------------------------------- 6 limbs ---- */

/* rsi=a  rbx=b  rcx=p  rbp=n0  rdi=r  t0..t6 = the accumulator
 * rax and r15 are the mulx destinations. */
.macro ROUND6 t0,t1,t2,t3,t4,t5,t6, off

    movq    \off(%rbx), %rdx            /* b[i] */
    xorl    %eax, %eax                  /* clears CF and OF together */

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

    movl    $0, %eax                    /* mov leaves CF and OF alone */
    adcxq   %rax, \t6                   /* carry out of t5, weight 2^384 */
                                        /* the adox carry is zero: see (2) */

    movq    \t0, %rdx
    imulq   %rbp, %rdx                  /* m = t0 * n0 mod 2^64 */
    xorl    %eax, %eax

    mulxq   0(%rcx), %rax, %r15
    adcxq   %rax, \t0                   /* t0 becomes zero, by choice of m */
    adoxq   %r15, \t1
    mulxq   8(%rcx), %rax, %r15
    adcxq   %rax, \t1
    adoxq   %r15, \t2
    mulxq   16(%rcx), %rax, %r15
    adcxq   %rax, \t2
    adoxq   %r15, \t3
    mulxq   24(%rcx), %rax, %r15
    adcxq   %rax, \t3
    adoxq   %r15, \t4
    mulxq   32(%rcx), %rax, %r15
    adcxq   %rax, \t4
    adoxq   %r15, \t5
    mulxq   40(%rcx), %rax, %r15
    adcxq   %rax, \t5
    adoxq   %r15, \t6

    movl    $0, %eax
    adcxq   %rax, \t6
.endm

    .globl SYM(elips_fp_mul_6_x86_64)
    FUNC_HIDDEN(elips_fp_mul_6_x86_64)
    FUNC_TYPE(elips_fp_mul_6_x86_64)
    .p2align 4
SYM(elips_fp_mul_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $56, %rsp                   /* 48 for the copy of T, 8 to pad */

    movq    %rdx, %rbx                  /* b, out of rdx before mulx needs it */
    movq    %r8, %rbp                   /* n0 */

    xorl    %r8d, %r8d                  /* T = 0 */
    xorl    %r9d, %r9d
    xorl    %r10d, %r10d
    xorl    %r11d, %r11d
    xorl    %r12d, %r12d
    xorl    %r13d, %r13d
    xorl    %r14d, %r14d

    ROUND6  %r8,  %r9,  %r10, %r11, %r12, %r13, %r14,  0
    ROUND6  %r9,  %r10, %r11, %r12, %r13, %r14, %r8,   8
    ROUND6  %r10, %r11, %r12, %r13, %r14, %r8,  %r9,   16
    ROUND6  %r11, %r12, %r13, %r14, %r8,  %r9,  %r10,  24
    ROUND6  %r12, %r13, %r14, %r8,  %r9,  %r10, %r11,  32
    ROUND6  %r13, %r14, %r8,  %r9,  %r10, %r11, %r12,  40

    /* T is now (r14, r8, r9, r10, r11, r12) and r13 is zero. */
    movq    %r14, 0(%rsp)
    movq    %r8,  8(%rsp)
    movq    %r9,  16(%rsp)
    movq    %r10, 24(%rsp)
    movq    %r11, 32(%rsp)
    movq    %r12, 40(%rsp)

    subq    0(%rcx),  %r14
    sbbq    8(%rcx),  %r8
    sbbq    16(%rcx), %r9
    sbbq    24(%rcx), %r10
    sbbq    32(%rcx), %r11
    sbbq    40(%rcx), %r12

    /* CF set means T < p, so put the saved copy back. cmov leaves flags. */
    cmovcq  0(%rsp),  %r14
    cmovcq  8(%rsp),  %r8
    cmovcq  16(%rsp), %r9
    cmovcq  24(%rsp), %r10
    cmovcq  32(%rsp), %r11
    cmovcq  40(%rsp), %r12

    movq    %r14, 0(%rdi)
    movq    %r8,  8(%rdi)
    movq    %r9,  16(%rdi)
    movq    %r10, 24(%rdi)
    movq    %r11, 32(%rdi)
    movq    %r12, 40(%rdi)

    addq    $56, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_mul_6_x86_64)

/* ----------------------------------------------------------- 8 limbs ---- */

/* Nine accumulator words plus two mulx destinations plus three pointers is
 * exactly the fifteen general registers x86-64 has, so the result pointer and
 * n0 spend the round loop on the stack. n0 is only ever an imul operand and
 * imul reads memory, so that costs nothing.
 *
 * rsi=a  rdi=b  rcx=p  0(%rsp)=r  8(%rsp)=n0  rax,rbp = mulx destinations */
.macro ROUND8 t0,t1,t2,t3,t4,t5,t6,t7,t8, off

    movq    \off(%rdi), %rdx            /* b[i] */
    xorl    %eax, %eax                  /* clears CF and OF together */

    mulxq   0(%rsi), %rax, %rbp
    adcxq   %rax, \t0
    adoxq   %rbp, \t1
    mulxq   8(%rsi), %rax, %rbp
    adcxq   %rax, \t1
    adoxq   %rbp, \t2
    mulxq   16(%rsi), %rax, %rbp
    adcxq   %rax, \t2
    adoxq   %rbp, \t3
    mulxq   24(%rsi), %rax, %rbp
    adcxq   %rax, \t3
    adoxq   %rbp, \t4
    mulxq   32(%rsi), %rax, %rbp
    adcxq   %rax, \t4
    adoxq   %rbp, \t5
    mulxq   40(%rsi), %rax, %rbp
    adcxq   %rax, \t5
    adoxq   %rbp, \t6
    mulxq   48(%rsi), %rax, %rbp
    adcxq   %rax, \t6
    adoxq   %rbp, \t7
    mulxq   56(%rsi), %rax, %rbp
    adcxq   %rax, \t7
    adoxq   %rbp, \t8

    movl    $0, %eax                    /* mov leaves CF and OF alone */
    adcxq   %rax, \t8                   /* the adox carry is zero: see (2) */

    movq    \t0, %rdx
    imulq   8(%rsp), %rdx               /* m = t0 * n0 mod 2^64 */
    xorl    %eax, %eax

    mulxq   0(%rcx), %rax, %rbp
    adcxq   %rax, \t0                   /* t0 becomes zero, by choice of m */
    adoxq   %rbp, \t1
    mulxq   8(%rcx), %rax, %rbp
    adcxq   %rax, \t1
    adoxq   %rbp, \t2
    mulxq   16(%rcx), %rax, %rbp
    adcxq   %rax, \t2
    adoxq   %rbp, \t3
    mulxq   24(%rcx), %rax, %rbp
    adcxq   %rax, \t3
    adoxq   %rbp, \t4
    mulxq   32(%rcx), %rax, %rbp
    adcxq   %rax, \t4
    adoxq   %rbp, \t5
    mulxq   40(%rcx), %rax, %rbp
    adcxq   %rax, \t5
    adoxq   %rbp, \t6
    mulxq   48(%rcx), %rax, %rbp
    adcxq   %rax, \t6
    adoxq   %rbp, \t7
    mulxq   56(%rcx), %rax, %rbp
    adcxq   %rax, \t7
    adoxq   %rbp, \t8

    movl    $0, %eax
    adcxq   %rax, \t8
.endm

    .globl SYM(elips_fp_mul_8_x86_64)
    FUNC_HIDDEN(elips_fp_mul_8_x86_64)
    FUNC_TYPE(elips_fp_mul_8_x86_64)
    .p2align 4
SYM(elips_fp_mul_8_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $88, %rsp                   /* 0: r, 8: n0, 16..79: copy of T */

    movq    %rdi, 0(%rsp)               /* r, not needed until the end */
    movq    %r8,  8(%rsp)               /* n0 */
    movq    %rdx, %rdi                  /* b, out of rdx before mulx wants it */

    xorl    %r8d,  %r8d                 /* T = 0 */
    xorl    %r9d,  %r9d
    xorl    %r10d, %r10d
    xorl    %r11d, %r11d
    xorl    %r12d, %r12d
    xorl    %r13d, %r13d
    xorl    %r14d, %r14d
    xorl    %r15d, %r15d
    xorl    %ebx,  %ebx

    ROUND8  %r8,  %r9,  %r10, %r11, %r12, %r13, %r14, %r15, %rbx,  0
    ROUND8  %r9,  %r10, %r11, %r12, %r13, %r14, %r15, %rbx, %r8,   8
    ROUND8  %r10, %r11, %r12, %r13, %r14, %r15, %rbx, %r8,  %r9,   16
    ROUND8  %r11, %r12, %r13, %r14, %r15, %rbx, %r8,  %r9,  %r10,  24
    ROUND8  %r12, %r13, %r14, %r15, %rbx, %r8,  %r9,  %r10, %r11,  32
    ROUND8  %r13, %r14, %r15, %rbx, %r8,  %r9,  %r10, %r11, %r12,  40
    ROUND8  %r14, %r15, %rbx, %r8,  %r9,  %r10, %r11, %r12, %r13,  48
    ROUND8  %r15, %rbx, %r8,  %r9,  %r10, %r11, %r12, %r13, %r14,  56

    /* T is now (rbx, r8..r14) and r15 is zero. */
    movq    %rbx, 16(%rsp)
    movq    %r8,  24(%rsp)
    movq    %r9,  32(%rsp)
    movq    %r10, 40(%rsp)
    movq    %r11, 48(%rsp)
    movq    %r12, 56(%rsp)
    movq    %r13, 64(%rsp)
    movq    %r14, 72(%rsp)

    subq    0(%rcx),  %rbx
    sbbq    8(%rcx),  %r8
    sbbq    16(%rcx), %r9
    sbbq    24(%rcx), %r10
    sbbq    32(%rcx), %r11
    sbbq    40(%rcx), %r12
    sbbq    48(%rcx), %r13
    sbbq    56(%rcx), %r14

    /* CF set means T < p, so put the saved copy back. cmov leaves flags. */
    cmovcq  16(%rsp), %rbx
    cmovcq  24(%rsp), %r8
    cmovcq  32(%rsp), %r9
    cmovcq  40(%rsp), %r10
    cmovcq  48(%rsp), %r11
    cmovcq  56(%rsp), %r12
    cmovcq  64(%rsp), %r13
    cmovcq  72(%rsp), %r14

    movq    0(%rsp), %rdi
    movq    %rbx, 0(%rdi)
    movq    %r8,  8(%rdi)
    movq    %r9,  16(%rdi)
    movq    %r10, 24(%rdi)
    movq    %r11, 32(%rdi)
    movq    %r12, 40(%rdi)
    movq    %r13, 48(%rdi)
    movq    %r14, 56(%rdi)

    addq    $88, %rsp
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_mul_8_x86_64)

#if defined(__linux__)
.section .note.GNU-stack,"",@progbits
#endif

#endif /* __x86_64__ */
