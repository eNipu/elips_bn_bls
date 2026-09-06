/*
 * The system CSPRNG. See include/elips/sysrand.h.
 *
 * Three sources, picked at compile time, in the order the platform prefers:
 *
 *   getrandom(2)      Linux 3.17+, glibc 2.25+. Blocks only until the pool is
 *                     initialised, and never after, so it is safe at start-up.
 *   arc4random_buf    macOS, the BSDs. Cannot fail and needs no descriptor.
 *   /dev/urandom      Everything else. Opened per call rather than cached,
 *                     because a cached descriptor is a file-descriptor leak in
 *                     a library and a use-after-close hazard across fork.
 *
 * There is deliberately no userspace generator and no seeding. Every failure
 * mode of the old time(NULL) design -- guessable seed, two processes agreeing,
 * a stale state after fork -- comes from having state to get wrong.
 */
#if defined(__linux__)
#  define _GNU_SOURCE 1
#endif

#include "elips/sysrand.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

#if defined(__linux__)
#  include <sys/random.h>
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || \
      defined(__NetBSD__) || defined(__DragonFly__)
#  include <stdlib.h>
#  define ELIPS_HAVE_ARC4RANDOM 1
#endif

int elips_random_bytes(void *buf, size_t n)
{
    unsigned char *p = (unsigned char *)buf;
    if (n == 0) return 0;

#if defined(ELIPS_HAVE_ARC4RANDOM)
    arc4random_buf(p, n);
    return 0;
#else
#  if defined(__linux__)
    size_t got = 0;
    while (got < n) {
        ssize_t rc = getrandom(p + got, n - got, 0);
        if (rc < 0) {
            if (errno == EINTR) continue;
            break;                      /* fall through to /dev/urandom */
        }
        got += (size_t)rc;
    }
    if (got == n) return 0;
#  endif
    /* Fallback: /dev/urandom. Also the path taken on a kernel too old for
     * getrandom, or inside a sandbox that blocks the syscall. */
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) { memset(p, 0, n); return -1; }
    size_t rd = fread(p, 1, n, f);
    fclose(f);
    if (rd != n) { memset(p, 0, n); return -1; }
    return 0;
#endif
}

