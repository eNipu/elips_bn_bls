/*
 * The system CSPRNG, with no dependency on a compiled-in curve.
 *
 * Split out from elips/random.h so both layers can use it: the fixed-width
 * Phase 3 arithmetic is built once per curve, while the legacy mpz layer picks
 * its curve at run time and cannot include anything that needs FP_LIMBS.
 *
 * The source is the operating system's, not a userspace generator: getrandom(2)
 * on Linux, arc4random_buf on the BSDs and macOS, and a /dev/urandom read as
 * the fallback. There is no seeding, no state to reseed after a fork, and no
 * failure mode where two processes agree.
 */
#ifndef ELIPS_SYSRAND_H
#define ELIPS_SYSRAND_H

#include <stddef.h>

/* Fill buf with n bytes from the system CSPRNG. Returns 0 on success and -1 on
 * failure, and a failure leaves the buffer zeroed rather than partly filled.
 * Callers must check: a silently non-random key is worse than an error. */
int elips_random_bytes(void *buf, size_t n);

#endif /* ELIPS_SYSRAND_H */
