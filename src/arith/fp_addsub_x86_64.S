/*
 * Modular addition and subtraction for x86-64.
 *
 *   void elips_fp_add_6_x86_64(limb_t r[6], const limb_t a[6],
 *                              const limb_t b[6], const limb_t p[6]);
 *   void elips_fp_sub_6_x86_64(...);   and the 8 limb pair.
 *
 * System V AMD64: rdi=r, rsi=a, rdx=b, rcx=p.
 *
 * WHY THIS EXISTS. The portable C makes three passes over the limbs -- add,
 * subtract p, select -- and builds its carry chain out of comparisons, because
 * C has no adc. Measured on an Ice Lake Xeon that is 24.18 cycles for a
 * 384-bit modular add and 26.43 for the subtract, against 70.21 for a whole
 * Montgomery multiply. A pairing spends 50% of its time in these two and 39%
 * in the multiply, so they were the larger half and the unoptimised one. See
 * issue #36 for the profile and #37 for this change.
 *
 * NO CPUID CHECK, deliberately, and this differs from fp_mul. adc, sbb and
 * cmov are baseline x86-64; mulx/adcx/adox are not. Gating this behind the
 * existing BMI2+ADX test would switch it off on precisely the older machines
 * that have no fast multiply either.
 *
 * CONSTANT TIME. No branches. fp_add selects with cmov between two live
 * registers; fp_sub adds p masked by 0 or ~0. Neither reads a secret address.
 * tools/verify/ct_branch_scan.py covers both shapes it knows about, and the
 * dudect targets fp_add and fp_sub run on every build.
 *
 * ALIASING. Every input limb is loaded before any output limb is stored, so
 * r == a, r == b and r == a == b are all safe. test/fp_difftest.c checks each
 * of them; a value-based test would not.
 */

#if defined(__x86_64__) && !defined(__ILP32__)

/* Symbol decoration, one directive per SOURCE line. A cpp macro cannot contain
 * a newline, and joining directives with ";" silently produces an object with
 * no symbols on arm64 Darwin, where ";" begins a comment. Same rule as
 * fp_mul_x86_64.S, and the same bug if it is broken. */
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

/* ----------------------------------------------------------- 6 limbs ----
 *
 * At six limbs the sum and the reduced candidate both fit in registers at the
 * same time, once a and b have been consumed and rsi and rdx freed. That is
 * 6 + 6 + carry + r + p = 15, which is exactly what is available.
 *
 *   sum       r8  r9  r10 r11 r12 r13
 *   candidate rax rsi rdx r14 r15 rbx
 *   rbp       carry out of the sum
 */

    .globl SYM(elips_fp_add_6_x86_64)
    FUNC_HIDDEN(elips_fp_add_6_x86_64)
    FUNC_TYPE(elips_fp_add_6_x86_64)
    .align 16
SYM(elips_fp_add_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    xorl    %ebp, %ebp                  /* also clears CF for the chain below */

    movq    0(%rsi),  %r8
    addq    0(%rdx),  %r8
    movq    8(%rsi),  %r9               /* mov does not touch the flags, so */
    adcq    8(%rdx),  %r9               /* the loads sit inside the chain   */
    movq    16(%rsi), %r10
    adcq    16(%rdx), %r10
    movq    24(%rsi), %r11
    adcq    24(%rdx), %r11
    movq    32(%rsi), %r12
    adcq    32(%rdx), %r12
    movq    40(%rsi), %r13
    adcq    40(%rdx), %r13
    adcq    $0, %rbp                    /* rbp = the 385th bit, 0 or 1 */

    /* a and b are fully read now, so their pointers become scratch. */
    movq    %r8,  %rax
    subq    0(%rcx),  %rax
    movq    %r9,  %rsi
    sbbq    8(%rcx),  %rsi
    movq    %r10, %rdx
    sbbq    16(%rcx), %rdx
    movq    %r11, %r14
    sbbq    24(%rcx), %r14
    movq    %r12, %r15
    sbbq    32(%rcx), %r15
    movq    %r13, %rbx
    sbbq    40(%rcx), %rbx

    /* rbp - 0 - borrow. It goes negative, so CF is set, exactly when the sum
     * did not overflow AND subtracting p borrowed -- which is exactly the case
     * where the sum was already below p and must be kept. */
    sbbq    $0, %rbp
    cmovcq  %r8,  %rax
    cmovcq  %r9,  %rsi
    cmovcq  %r10, %rdx
    cmovcq  %r11, %r14
    cmovcq  %r12, %r15
    cmovcq  %r13, %rbx

    movq    %rax, 0(%rdi)
    movq    %rsi, 8(%rdi)
    movq    %rdx, 16(%rdi)
    movq    %r14, 24(%rdi)
    movq    %r15, 32(%rdi)
    movq    %rbx, 40(%rdi)

    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_add_6_x86_64)

    .globl SYM(elips_fp_sub_6_x86_64)
    FUNC_HIDDEN(elips_fp_sub_6_x86_64)
    FUNC_TYPE(elips_fp_sub_6_x86_64)
    .align 16
SYM(elips_fp_sub_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    movq    0(%rsi),  %r8
    subq    0(%rdx),  %r8
    movq    8(%rsi),  %r9
    sbbq    8(%rdx),  %r9
    movq    16(%rsi), %r10
    sbbq    16(%rdx), %r10
    movq    24(%rsi), %r11
    sbbq    24(%rdx), %r11
    movq    32(%rsi), %r12
    sbbq    32(%rdx), %r12
    movq    40(%rsi), %r13
    sbbq    40(%rdx), %r13
    sbbq    %rbp, %rbp                  /* rbp = 0 - borrow: 0 or all ones */

    /* Mask p FIRST and add afterwards. and clears CF, so a masked limb
     * computed inside the carry chain would break it. */
    movq    0(%rcx),  %rax
    andq    %rbp, %rax
    movq    8(%rcx),  %rsi
    andq    %rbp, %rsi
    movq    16(%rcx), %rdx
    andq    %rbp, %rdx
    movq    24(%rcx), %r14
    andq    %rbp, %r14
    movq    32(%rcx), %r15
    andq    %rbp, %r15
    movq    40(%rcx), %rbx
    andq    %rbp, %rbx

    addq    %rax, %r8
    adcq    %rsi, %r9
    adcq    %rdx, %r10
    adcq    %r14, %r11
    adcq    %r15, %r12
    adcq    %rbx, %r13

    movq    %r8,  0(%rdi)
    movq    %r9,  8(%rdi)
    movq    %r10, 16(%rdi)
    movq    %r11, 24(%rdi)
    movq    %r12, 32(%rdi)
    movq    %r13, 40(%rdi)

    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_sub_6_x86_64)


/* ----------------------------------------------------------- 8 limbs ----
 *
 * Eight limbs do not leave room for the sum and the candidate at once: that
 * would be 16 values plus r and p, and there are 15 registers. So the sum
 * stays in r8..r15 and the candidate is written to r[] and selected back out
 * of it with a memory-source cmov. r[] is safe to use as scratch because a and
 * b are both fully read by then, which is also what makes aliasing work.
 *
 * A few cycles slower than the six limb path, and it affects only BLS12-461
 * and BN-462. BLS12-381 is six.
 */

    .globl SYM(elips_fp_add_8_x86_64)
    FUNC_HIDDEN(elips_fp_add_8_x86_64)
    FUNC_TYPE(elips_fp_add_8_x86_64)
    .align 16
SYM(elips_fp_add_8_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    xorl    %ebp, %ebp

    movq    0(%rsi),  %r8
    addq    0(%rdx),  %r8
    movq    8(%rsi),  %r9
    adcq    8(%rdx),  %r9
    movq    16(%rsi), %r10
    adcq    16(%rdx), %r10
    movq    24(%rsi), %r11
    adcq    24(%rdx), %r11
    movq    32(%rsi), %r12
    adcq    32(%rdx), %r12
    movq    40(%rsi), %r13
    adcq    40(%rdx), %r13
    movq    48(%rsi), %r14
    adcq    48(%rdx), %r14
    movq    56(%rsi), %r15
    adcq    56(%rdx), %r15
    adcq    $0, %rbp

    /* candidate = sum - p, parked in r[]. Stores do not touch the flags, so
     * they sit inside the borrow chain. */
    movq    %r8,  %rax
    subq    0(%rcx),  %rax
    movq    %rax, 0(%rdi)
    movq    %r9,  %rax
    sbbq    8(%rcx),  %rax
    movq    %rax, 8(%rdi)
    movq    %r10, %rax
    sbbq    16(%rcx), %rax
    movq    %rax, 16(%rdi)
    movq    %r11, %rax
    sbbq    24(%rcx), %rax
    movq    %rax, 24(%rdi)
    movq    %r12, %rax
    sbbq    32(%rcx), %rax
    movq    %rax, 32(%rdi)
    movq    %r13, %rax
    sbbq    40(%rcx), %rax
    movq    %rax, 40(%rdi)
    movq    %r14, %rax
    sbbq    48(%rcx), %rax
    movq    %rax, 48(%rdi)
    movq    %r15, %rax
    sbbq    56(%rcx), %rax
    movq    %rax, 56(%rdi)

    /* CF here means the sum was already below p, so keep it; otherwise take
     * the candidate back out of r[]. cmov always reads its memory operand. */
    sbbq    $0, %rbp
    cmovncq 0(%rdi),  %r8
    cmovncq 8(%rdi),  %r9
    cmovncq 16(%rdi), %r10
    cmovncq 24(%rdi), %r11
    cmovncq 32(%rdi), %r12
    cmovncq 40(%rdi), %r13
    cmovncq 48(%rdi), %r14
    cmovncq 56(%rdi), %r15

    movq    %r8,  0(%rdi)
    movq    %r9,  8(%rdi)
    movq    %r10, 16(%rdi)
    movq    %r11, 24(%rdi)
    movq    %r12, 32(%rdi)
    movq    %r13, 40(%rdi)
    movq    %r14, 48(%rdi)
    movq    %r15, 56(%rdi)

    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_add_8_x86_64)

    .globl SYM(elips_fp_sub_8_x86_64)
    FUNC_HIDDEN(elips_fp_sub_8_x86_64)
    FUNC_TYPE(elips_fp_sub_8_x86_64)
    .align 16
SYM(elips_fp_sub_8_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    movq    0(%rsi),  %r8
    subq    0(%rdx),  %r8
    movq    8(%rsi),  %r9
    sbbq    8(%rdx),  %r9
    movq    16(%rsi), %r10
    sbbq    16(%rdx), %r10
    movq    24(%rsi), %r11
    sbbq    24(%rdx), %r11
    movq    32(%rsi), %r12
    sbbq    32(%rdx), %r12
    movq    40(%rsi), %r13
    sbbq    40(%rdx), %r13
    movq    48(%rsi), %r14
    sbbq    48(%rdx), %r14
    movq    56(%rsi), %r15
    sbbq    56(%rdx), %r15
    sbbq    %rbp, %rbp                  /* 0, or all ones if it borrowed */

    /* candidate = difference + p, parked in r[]. The borrow is already safe
     * in rbp, so starting a fresh carry chain here costs nothing. */
    movq    %r8,  %rax
    addq    0(%rcx),  %rax
    movq    %rax, 0(%rdi)
    movq    %r9,  %rax
    adcq    8(%rcx),  %rax
    movq    %rax, 8(%rdi)
    movq    %r10, %rax
    adcq    16(%rcx), %rax
    movq    %rax, 16(%rdi)
    movq    %r11, %rax
    adcq    24(%rcx), %rax
    movq    %rax, 24(%rdi)
    movq    %r12, %rax
    adcq    32(%rcx), %rax
    movq    %rax, 32(%rdi)
    movq    %r13, %rax
    adcq    40(%rcx), %rax
    movq    %rax, 40(%rdi)
    movq    %r14, %rax
    adcq    48(%rcx), %rax
    movq    %rax, 48(%rdi)
    movq    %r15, %rax
    adcq    56(%rcx), %rax
    movq    %rax, 56(%rdi)

    /* rbp is zero exactly when nothing borrowed, and then the difference
     * stands. Otherwise take the one with p added back. */
    testq   %rbp, %rbp
    cmovnzq 0(%rdi),  %r8
    cmovnzq 8(%rdi),  %r9
    cmovnzq 16(%rdi), %r10
    cmovnzq 24(%rdi), %r11
    cmovnzq 32(%rdi), %r12
    cmovnzq 40(%rdi), %r13
    cmovnzq 48(%rdi), %r14
    cmovnzq 56(%rdi), %r15

    movq    %r8,  0(%rdi)
    movq    %r9,  8(%rdi)
    movq    %r10, 16(%rdi)
    movq    %r11, 24(%rdi)
    movq    %r12, 32(%rdi)
    movq    %r13, 40(%rdi)
    movq    %r14, 48(%rdi)
    movq    %r15, 56(%rdi)

    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbp
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_sub_8_x86_64)

#endif /* __x86_64__ */

/* OUTSIDE the architecture guard, deliberately. On another architecture this
 * file is an empty object, and an ELF object with no .note.GNU-stack tells GNU
 * ld the result needs an EXECUTABLE STACK -- which is how the Python extension
 * once ended up with RWE permissions. Same reason as fp_mul_x86_64.S. */
#if defined(__linux__)
.section .note.GNU-stack,"",@progbits
#endif
