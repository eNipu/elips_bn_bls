/*
 * Addition and subtraction on 768-bit values modulo p*R, shared by the wide
 * combination routines of the tower.
 *
 * These are assembler macros, not C: the file is included from .S sources,
 * which the compiler runs through cpp first. They are here rather than copied
 * into each caller because they are the carry-critical part of the lazy path
 * and a mistake in one copy would be invisible in the other.
 *
 * Contract for the includer:
 *
 *   - %rcx holds p, the six-limb modulus.
 *   - SC and PM name two 48-byte scratch slots on the stack, in the frame the
 *     includer allocated.
 *   - %rax and %r8-%r14 are clobbered.
 *
 * Everything stays mod p*R, and p*R has six zero low limbs, so a borrow can
 * never leave the low half: the low six limbs are a plain subtract and only the
 * top six take the conditional add of p.
 */

/* p masked by %rax, parked in scratch. andq sets flags, so it cannot appear
 * between the addq and adcq of a carry chain; doing it here keeps the chain
 * contiguous. */
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

/* dst = a - b  mod p*R */
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

/* dst = a + b  mod p*R */
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
