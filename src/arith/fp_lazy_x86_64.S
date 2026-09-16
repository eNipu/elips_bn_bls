/*
 * Raw product and standalone Montgomery reduction: the two primitives the
 * lazy tower needs and the fused CIOS multiply cannot provide.
 *
 *   void elips_fp_mulw_6_x86_64 (limb_t w[12], const limb_t a[6],
 *                                const limb_t b[6]);
 *   void elips_fp_redcw_6_x86_64(limb_t r[6], const limb_t w[12],
 *                                const limb_t p[6], limb_t n0);
 *
 * WHY THESE EXIST. Lazy reduction replaces three fused multiplies with three
 * raw products and two reductions. Written in portable C that loses outright:
 * a raw 384x384 product in C measures 44.6 ns against 34.5 ns for the whole
 * fused assembly multiply-and-reduce, so three of them cost more than the
 * eager path spends in total. The algorithm only pays if both halves are as
 * well optimised as the fused routine they replace.
 *
 * Both routines are the halves of ROUND6 in fp_mul_x86_64.S, separated. Same
 * two carry chains, same register rotation, same sparse-modulus argument. Read
 * that file first; everything it says about adcx/adox and about why the top
 * word cannot overflow applies here unchanged.
 *
 * THE REDUCTION IS SPLIT, which is what makes it cheap. For w = lo + hi*R,
 *
 *     REDC(w) = w*R^-1 = lo*R^-1 + hi   (mod p)
 *
 * so reducing 768 bits is reducing the low 384 and adding the high 384. The
 * low half goes through six rounds of "m = t0*n0; t += m*p", which is ROUND6
 * without its multiply. blst splits it the same way, which is why its trace
 * shows __mulx_by_1_mont_384 and redcx_mont_384 in step.
 *
 * BOUNDS, and they are the whole correctness argument.
 *
 *   mulw: before round i the accumulator holds six significant words, so
 *   T < 2^384. a*b[i] <= (2^384-1)(2^64-1), and the sum is at most
 *   2^448 - 2^64, which fits the seven words exactly. Nothing leaves the top.
 *
 *   redcw: reducing the low half gives a value < p + 1, and the caller must
 *   supply w < p*R so that hi < p. The sum is then < 2p + 1 and ONE
 *   conditional subtraction lands in [0, p). Callers that bias a negative
 *   Karatsuba combination by a multiple of p^2 must keep the biased value
 *   under p*R for this to hold; p*2^384 is exactly p*R and is NOT a legal
 *   bias.
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

/* ------------------------------------------------ raw 384x384 -> 768 ---- */

/* T += a * b[i], then t0 is final: store it and let it rotate to the top as
 * the zero the next round needs. rsi=a rbx=b rdi=w */
.macro MULROUND6 t0,t1,t2,t3,t4,t5,t6, off

    movq    \off(%rbx), %rdx
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
    adcxq   %rax, \t6
    adoxq   %rax, \t6

    movq    \t0, \off(%rdi)             /* w[i] is final */
    xorl    %eax, %eax
    movq    %rax, \t0                   /* becomes the top on rotation */
.endm

    .globl SYM(elips_fp_mulw_6_x86_64)
    FUNC_HIDDEN(elips_fp_mulw_6_x86_64)
    FUNC_TYPE(elips_fp_mulw_6_x86_64)
    .p2align 4
SYM(elips_fp_mulw_6_x86_64):
    pushq   %rbx
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    movq    %rdx, %rbx                  /* b, out of rdx before mulx needs it */

    xorl    %r8d,  %r8d
    xorl    %r9d,  %r9d
    xorl    %r10d, %r10d
    xorl    %r11d, %r11d
    xorl    %r12d, %r12d
    xorl    %r13d, %r13d
    xorl    %r14d, %r14d

    MULROUND6  %r8,  %r9,  %r10, %r11, %r12, %r13, %r14,  0
    MULROUND6  %r9,  %r10, %r11, %r12, %r13, %r14, %r8,   8
    MULROUND6  %r10, %r11, %r12, %r13, %r14, %r8,  %r9,   16
    MULROUND6  %r11, %r12, %r13, %r14, %r8,  %r9,  %r10,  24
    MULROUND6  %r12, %r13, %r14, %r8,  %r9,  %r10, %r11,  32
    MULROUND6  %r13, %r14, %r8,  %r9,  %r10, %r11, %r12,  40

    /* the high half is (r14, r8, r9, r10, r11, r12) */
    movq    %r14, 48(%rdi)
    movq    %r8,  56(%rdi)
    movq    %r9,  64(%rdi)
    movq    %r10, 72(%rdi)
    movq    %r11, 80(%rdi)
    movq    %r12, 88(%rdi)

    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    popq    %rbx
    ret
    FUNC_SIZE(elips_fp_mulw_6_x86_64)

/* --------------------------------------------- 768 -> 384 reduction ---- */

/* One "m = t0*n0; T += m*p" round. This is ROUND6's second half verbatim.
 * rcx=p rbp=n0 */
.macro REDCROUND6 t0,t1,t2,t3,t4,t5,t6

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

    .globl SYM(elips_fp_redcw_6_x86_64)
    FUNC_HIDDEN(elips_fp_redcw_6_x86_64)
    FUNC_TYPE(elips_fp_redcw_6_x86_64)
    .p2align 4
SYM(elips_fp_redcw_6_x86_64):
    pushq   %rbx
    pushq   %rbp
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    subq    $56, %rsp

    movq    %rsi, %rbx                  /* w */
    movq    %rcx, %rbp                  /* n0 (4th arg) */
    movq    %rdx, %rcx                  /* p  (3rd arg) */

    movq    0(%rbx),  %r8               /* T = low half */
    movq    8(%rbx),  %r9
    movq    16(%rbx), %r10
    movq    24(%rbx), %r11
    movq    32(%rbx), %r12
    movq    40(%rbx), %r13
    xorl    %r14d, %r14d

    REDCROUND6  %r8,  %r9,  %r10, %r11, %r12, %r13, %r14
    REDCROUND6  %r9,  %r10, %r11, %r12, %r13, %r14, %r8
    REDCROUND6  %r10, %r11, %r12, %r13, %r14, %r8,  %r9
    REDCROUND6  %r11, %r12, %r13, %r14, %r8,  %r9,  %r10
    REDCROUND6  %r12, %r13, %r14, %r8,  %r9,  %r10, %r11
    REDCROUND6  %r13, %r14, %r8,  %r9,  %r10, %r11, %r12

    /* lo*R^-1 is (r14, r8, r9, r10, r11, r12); add the high half */
    addq    48(%rbx), %r14
    adcq    56(%rbx), %r8
    adcq    64(%rbx), %r9
    adcq    72(%rbx), %r10
    adcq    80(%rbx), %r11
    adcq    88(%rbx), %r12
    sbbq    %r13, %r13                  /* r13 = -(carry out) */

    movq    %r14, 0(%rsp)
    movq    %r8,  8(%rsp)
    movq    %r9,  16(%rsp)
    movq    %r10, 24(%rsp)
    movq    %r11, 32(%rsp)
    movq    %r12, 40(%rsp)

    subq    0(%rcx),  %r14              /* candidate = T - p */
    sbbq    8(%rcx),  %r8
    sbbq    16(%rcx), %r9
    sbbq    24(%rcx), %r10
    sbbq    32(%rcx), %r11
    sbbq    40(%rcx), %r12
    sbbq    $0, %r13                    /* borrow against the carry above */

    /* r13 == 0 means the subtraction was valid: keep the candidate. */
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
    FUNC_SIZE(elips_fp_redcw_6_x86_64)

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif

#endif /* __x86_64__ */
