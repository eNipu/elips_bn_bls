/*
 * Montgomery multiplication for AArch64.
 *
 * Same CIOS algorithm as the portable C in src/arith/fp.c and as the x86-64
 * assembly next to this file, so all three can be compared word for word.
 *
 *   void elips_fp_mul_6_aarch64(limb_t r[6], const limb_t a[6],
 *                               const limb_t b[6], const limb_t p[6],
 *                               limb_t n0);
 *
 * AAPCS64: x0=r, x1=a, x2=b, x3=p, x4=n0.
 *
 * AArch64 has no second carry chain to pair with the first, the way x86-64's
 * adcx and adox do, so the trick here is a different one. mul and umulh do not
 * touch the flags. That means the multiplies for the NEXT limbs can be issued
 * in between the adds and adcs that consume the previous ones, hiding the
 * multiply latency inside a carry chain that never breaks. The low halves go
 * in one pass and the high halves, shifted up one word, in a second.
 *
 * Two facts the code relies on, both proved rather than assumed:
 *
 *   1. a < p and T < 2p on entry to every round, so T stays below 2p and ONE
 *      conditional subtraction is enough at the end. This is the standard CIOS
 *      invariant and needs a < p, which every caller satisfies because
 *      Montgomery representatives are reduced.
 *
 *   2. The n+1 word accumulator never overflows, so the high-half pass never
 *      carries out of the top word and that carry is not added anywhere. This
 *      needs a SPARSE modulus, p < 2^(64n-1): 381 < 383 for BLS12-381, and 461
 *      and 462 are both < 511 for the 8 limb curves. src/arith/fp.c
 *      static-asserts it.
 *
 * Rounds rotate register names instead of moving data. t0 is exactly zero
 * after its round's reduction step, which is what the multiplier m is chosen
 * to arrange, so the word that leaves the bottom is already the zero the top
 * needs. The shift by one limb costs nothing.
 *
 * Constant time by construction: no branches at all, and the only conditional
 * is csel, which reads registers.
 */

#if defined(__aarch64__) || defined(__arm64__)

/* Apple's assembler wants a leading underscore and has no .type/.size. */
#if defined(__APPLE__)
#  define FUNC_BEGIN(x) .globl _##x; .p2align 4; _##x:
#  define FUNC_END(x)
#else
#  define FUNC_BEGIN(x) .globl x; .hidden x; .type x,@function; .p2align 4; x:
#  define FUNC_END(x) .size x, .-x
#endif

/* x18 is the platform register on Apple and Windows and is never touched. */

.text

/* ----------------------------------------------------------- 6 limbs ----
 *
 * a and p both fit in registers at this width, so the round body touches
 * memory only for b[i]. x12 holds b[i] and then m; x13 to x15 are the
 * multiply destinations, rotated so each is consumed before it is rewritten.
 *
 * a0..a5 = x16 x17 x19 x20 x21 x22     p0..p5 = x23 x24 x25 x26 x27 x28
 */
.macro ROUND6 t0,t1,t2,t3,t4,t5,t6, off

    ldr     x12, [x2, #\off]            /* b[i] */

    mul     x13, x16, x12
    mul     x14, x17, x12
    mul     x15, x19, x12
    adds    \t0, \t0, x13
    mul     x13, x20, x12
    adcs    \t1, \t1, x14
    mul     x14, x21, x12
    adcs    \t2, \t2, x15
    mul     x15, x22, x12
    adcs    \t3, \t3, x13
    adcs    \t4, \t4, x14
    adcs    \t5, \t5, x15
    adc     \t6, \t6, xzr

    umulh   x13, x16, x12
    umulh   x14, x17, x12
    umulh   x15, x19, x12
    adds    \t1, \t1, x13
    umulh   x13, x20, x12
    adcs    \t2, \t2, x14
    umulh   x14, x21, x12
    adcs    \t3, \t3, x15
    umulh   x15, x22, x12
    adcs    \t4, \t4, x13
    adcs    \t5, \t5, x14
    adcs    \t6, \t6, x15
                                        /* no carry out of t6: see (2) */

    mul     x12, \t0, x4                /* m = t0 * n0 mod 2^64 */

    mul     x13, x23, x12
    mul     x14, x24, x12
    mul     x15, x25, x12
    adds    \t0, \t0, x13               /* t0 becomes zero, by choice of m */
    mul     x13, x26, x12
    adcs    \t1, \t1, x14
    mul     x14, x27, x12
    adcs    \t2, \t2, x15
    mul     x15, x28, x12
    adcs    \t3, \t3, x13
    adcs    \t4, \t4, x14
    adcs    \t5, \t5, x15
    adc     \t6, \t6, xzr

    umulh   x13, x23, x12
    umulh   x14, x24, x12
    umulh   x15, x25, x12
    adds    \t1, \t1, x13
    umulh   x13, x26, x12
    adcs    \t2, \t2, x14
    umulh   x14, x27, x12
    adcs    \t3, \t3, x15
    umulh   x15, x28, x12
    adcs    \t4, \t4, x13
    adcs    \t5, \t5, x14
    adcs    \t6, \t6, x15
.endm

FUNC_BEGIN(elips_fp_mul_6_aarch64)
    stp     x19, x20, [sp, #-80]!
    stp     x21, x22, [sp, #16]
    stp     x23, x24, [sp, #32]
    stp     x25, x26, [sp, #48]
    stp     x27, x28, [sp, #64]

    ldp     x16, x17, [x1, #0]          /* a */
    ldp     x19, x20, [x1, #16]
    ldp     x21, x22, [x1, #32]
    ldp     x23, x24, [x3, #0]          /* p */
    ldp     x25, x26, [x3, #16]
    ldp     x27, x28, [x3, #32]

    mov     x5,  xzr                    /* T = 0 */
    mov     x6,  xzr
    mov     x7,  xzr
    mov     x8,  xzr
    mov     x9,  xzr
    mov     x10, xzr
    mov     x11, xzr

    ROUND6  x5,  x6,  x7,  x8,  x9,  x10, x11,  0
    ROUND6  x6,  x7,  x8,  x9,  x10, x11, x5,   8
    ROUND6  x7,  x8,  x9,  x10, x11, x5,  x6,   16
    ROUND6  x8,  x9,  x10, x11, x5,  x6,  x7,   24
    ROUND6  x9,  x10, x11, x5,  x6,  x7,  x8,   32
    ROUND6  x10, x11, x5,  x6,  x7,  x8,  x9,   40

    /* T is now (x11, x5, x6, x7, x8, x9) and x10 is zero. The a registers are
     * finished with, so they hold the difference. */
    subs    x12, x11, x23
    sbcs    x13, x5,  x24
    sbcs    x14, x6,  x25
    sbcs    x15, x7,  x26
    sbcs    x16, x8,  x27
    sbcs    x17, x9,  x28

    /* On AArch64 subs CLEARS the carry on borrow, so cc means T < p and the
     * original is the answer. csel leaves the flags alone. */
    csel    x11, x11, x12, cc
    csel    x5,  x5,  x13, cc
    csel    x6,  x6,  x14, cc
    csel    x7,  x7,  x15, cc
    csel    x8,  x8,  x16, cc
    csel    x9,  x9,  x17, cc

    stp     x11, x5, [x0, #0]
    stp     x6,  x7, [x0, #16]
    stp     x8,  x9, [x0, #32]

    ldp     x21, x22, [sp, #16]
    ldp     x23, x24, [sp, #32]
    ldp     x25, x26, [sp, #48]
    ldp     x27, x28, [sp, #64]
    ldp     x19, x20, [sp], #80
    ret
FUNC_END(elips_fp_mul_6_aarch64)

/* ----------------------------------------------------------- 8 limbs ----
 *
 * Nine accumulator words leave room to keep a in registers but not p, so p is
 * read back with ldp inside the round. ldp does not touch the flags, so the
 * loads sit inside the carry chain for free, the same way the multiplies do.
 *
 * a0..a7 = x19..x26     x14 = b[i] then m     x15 x16 x17 x27 x28 = temporaries
 */
.macro ROUND8 t0,t1,t2,t3,t4,t5,t6,t7,t8, off

    ldr     x14, [x2, #\off]            /* b[i] */

    mul     x15, x19, x14
    mul     x16, x20, x14
    mul     x17, x21, x14
    adds    \t0, \t0, x15
    mul     x15, x22, x14
    adcs    \t1, \t1, x16
    mul     x16, x23, x14
    adcs    \t2, \t2, x17
    mul     x17, x24, x14
    adcs    \t3, \t3, x15
    mul     x15, x25, x14
    adcs    \t4, \t4, x16
    mul     x16, x26, x14
    adcs    \t5, \t5, x17
    adcs    \t6, \t6, x15
    adcs    \t7, \t7, x16
    adc     \t8, \t8, xzr

    umulh   x15, x19, x14
    umulh   x16, x20, x14
    umulh   x17, x21, x14
    adds    \t1, \t1, x15
    umulh   x15, x22, x14
    adcs    \t2, \t2, x16
    umulh   x16, x23, x14
    adcs    \t3, \t3, x17
    umulh   x17, x24, x14
    adcs    \t4, \t4, x15
    umulh   x15, x25, x14
    adcs    \t5, \t5, x16
    umulh   x16, x26, x14
    adcs    \t6, \t6, x17
    adcs    \t7, \t7, x15
    adcs    \t8, \t8, x16
                                        /* no carry out of t8: see (2) */

    mul     x14, \t0, x4                /* m = t0 * n0 mod 2^64 */

    ldp     x27, x28, [x3, #0]
    mul     x15, x27, x14
    mul     x16, x28, x14
    ldp     x27, x28, [x3, #16]
    mul     x17, x27, x14
    adds    \t0, \t0, x15               /* t0 becomes zero, by choice of m */
    mul     x15, x28, x14
    adcs    \t1, \t1, x16
    ldp     x27, x28, [x3, #32]
    mul     x16, x27, x14
    adcs    \t2, \t2, x17
    mul     x17, x28, x14
    adcs    \t3, \t3, x15
    ldp     x27, x28, [x3, #48]
    mul     x15, x27, x14
    adcs    \t4, \t4, x16
    mul     x16, x28, x14
    adcs    \t5, \t5, x17
    adcs    \t6, \t6, x15
    adcs    \t7, \t7, x16
    adc     \t8, \t8, xzr

    ldp     x27, x28, [x3, #0]
    umulh   x15, x27, x14
    umulh   x16, x28, x14
    ldp     x27, x28, [x3, #16]
    umulh   x17, x27, x14
    adds    \t1, \t1, x15
    umulh   x15, x28, x14
    adcs    \t2, \t2, x16
    ldp     x27, x28, [x3, #32]
    umulh   x16, x27, x14
    adcs    \t3, \t3, x17
    umulh   x17, x28, x14
    adcs    \t4, \t4, x15
    ldp     x27, x28, [x3, #48]
    umulh   x15, x27, x14
    adcs    \t5, \t5, x16
    umulh   x16, x28, x14
    adcs    \t6, \t6, x17
    adcs    \t7, \t7, x15
    adcs    \t8, \t8, x16
.endm

FUNC_BEGIN(elips_fp_mul_8_aarch64)
    stp     x19, x20, [sp, #-80]!
    stp     x21, x22, [sp, #16]
    stp     x23, x24, [sp, #32]
    stp     x25, x26, [sp, #48]
    stp     x27, x28, [sp, #64]

    ldp     x19, x20, [x1, #0]          /* a */
    ldp     x21, x22, [x1, #16]
    ldp     x23, x24, [x1, #32]
    ldp     x25, x26, [x1, #48]

    mov     x5,  xzr                    /* T = 0 */
    mov     x6,  xzr
    mov     x7,  xzr
    mov     x8,  xzr
    mov     x9,  xzr
    mov     x10, xzr
    mov     x11, xzr
    mov     x12, xzr
    mov     x13, xzr

    ROUND8  x5,  x6,  x7,  x8,  x9,  x10, x11, x12, x13,  0
    ROUND8  x6,  x7,  x8,  x9,  x10, x11, x12, x13, x5,   8
    ROUND8  x7,  x8,  x9,  x10, x11, x12, x13, x5,  x6,   16
    ROUND8  x8,  x9,  x10, x11, x12, x13, x5,  x6,  x7,   24
    ROUND8  x9,  x10, x11, x12, x13, x5,  x6,  x7,  x8,   32
    ROUND8  x10, x11, x12, x13, x5,  x6,  x7,  x8,  x9,   40
    ROUND8  x11, x12, x13, x5,  x6,  x7,  x8,  x9,  x10,  48
    ROUND8  x12, x13, x5,  x6,  x7,  x8,  x9,  x10, x11,  56

    /* T is now (x13, x5..x11) and x12 is zero. The a registers are finished
     * with, so they take p for the one conditional subtraction. */
    ldp     x19, x20, [x3, #0]
    ldp     x21, x22, [x3, #16]
    ldp     x23, x24, [x3, #32]
    ldp     x25, x26, [x3, #48]

    subs    x15, x13, x19
    sbcs    x16, x5,  x20
    sbcs    x17, x6,  x21
    sbcs    x27, x7,  x22
    sbcs    x28, x8,  x23
    sbcs    x1,  x9,  x24
    sbcs    x2,  x10, x25
    sbcs    x4,  x11, x26

    /* On AArch64 subs CLEARS the carry on borrow, so cc means T < p and the
     * original is the answer. csel leaves the flags alone. */
    csel    x13, x13, x15, cc
    csel    x5,  x5,  x16, cc
    csel    x6,  x6,  x17, cc
    csel    x7,  x7,  x27, cc
    csel    x8,  x8,  x28, cc
    csel    x9,  x9,  x1,  cc
    csel    x10, x10, x2,  cc
    csel    x11, x11, x4,  cc

    stp     x13, x5,  [x0, #0]
    stp     x6,  x7,  [x0, #16]
    stp     x8,  x9,  [x0, #32]
    stp     x10, x11, [x0, #48]

    ldp     x21, x22, [sp, #16]
    ldp     x23, x24, [sp, #32]
    ldp     x25, x26, [sp, #48]
    ldp     x27, x28, [sp, #64]
    ldp     x19, x20, [sp], #80
    ret
FUNC_END(elips_fp_mul_8_aarch64)

#if defined(__linux__)
.section .note.GNU-stack,"",%progbits
#endif

#endif /* __aarch64__ */
