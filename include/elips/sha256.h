/*
 * SHA-256.
 *
 * Written here rather than linked from OpenSSL because RFC 9380's
 * expand_message_xmd is the only thing in this library that needs a hash, and
 * a whole TLS stack is a large dependency to acquire for one compression
 * function. It is also the only way the build stays a single library plus GMP.
 *
 * Pinned by the NIST vectors in test/h2c_test.c, and transitively by every
 * hash-to-curve vector.
 */
#ifndef ELIPS_SHA256_H
#define ELIPS_SHA256_H

#include <stddef.h>
#include <stdint.h>

#define ELIPS_SHA256_DIGEST 32
#define ELIPS_SHA256_BLOCK  64

typedef struct {
    uint32_t h[8];
    uint64_t len;                          /* message bytes absorbed so far */
    uint8_t  buf[ELIPS_SHA256_BLOCK];
    size_t   buf_len;
} elips_sha256_t;

void elips_sha256_init(elips_sha256_t *s);
void elips_sha256_update(elips_sha256_t *s, const void *data, size_t len);
void elips_sha256_final(uint8_t out[ELIPS_SHA256_DIGEST], elips_sha256_t *s);
void elips_sha256(uint8_t out[ELIPS_SHA256_DIGEST], const void *data, size_t len);

#endif /* ELIPS_SHA256_H */
