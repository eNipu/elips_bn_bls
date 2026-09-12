/*
 * HMAC-SHA256 and HKDF-SHA256, RFC 2104 and RFC 5869.
 *
 * Internal to the BLS layer, which needs them for exactly one thing: the
 * KeyGen of draft-irtf-cfrg-bls-signature derives a secret key from input
 * keying material with HKDF. Nothing else in the library hashes with a key.
 *
 * Built on elips/sha256.h, which is public, so this file is part of the
 * elips_bls target and the core library never sees it.
 *
 * Pinned by the RFC 5869 test vectors in test/bls_test.c. A key derivation
 * checked only against itself would agree with whatever it happens to do.
 */
#ifndef ELIPS_BLS_HKDF_H
#define ELIPS_BLS_HKDF_H

#include <stddef.h>
#include <stdint.h>

#include "elips/sha256.h"

#define ELIPS_HMAC_BYTES ELIPS_SHA256_DIGEST

void elips_hmac_sha256(uint8_t out[ELIPS_HMAC_BYTES],
                       const uint8_t *key, size_t key_len,
                       const uint8_t *msg, size_t msg_len);

/* PRK = HMAC(salt, ikm). RFC 5869 section 2.2. */
void elips_hkdf_extract(uint8_t prk[ELIPS_HMAC_BYTES],
                        const uint8_t *salt, size_t salt_len,
                        const uint8_t *ikm, size_t ikm_len);

/* OKM = the first out_len bytes of T(1) || T(2) || ... RFC 5869 section 2.3.
 * Returns -1 if out_len exceeds the 255*HashLen the construction allows. */
int elips_hkdf_expand(uint8_t *out, size_t out_len,
                      const uint8_t prk[ELIPS_HMAC_BYTES],
                      const uint8_t *info, size_t info_len);

#endif /* ELIPS_BLS_HKDF_H */
