/* See src/bls/hkdf.h. */
#include "hkdf.h"

#include <string.h>

void elips_hmac_sha256(uint8_t out[ELIPS_HMAC_BYTES],
                       const uint8_t *key, size_t key_len,
                       const uint8_t *msg, size_t msg_len)
{
    uint8_t k[ELIPS_SHA256_BLOCK];
    uint8_t pad[ELIPS_SHA256_BLOCK];
    uint8_t inner[ELIPS_SHA256_DIGEST];
    elips_sha256_t s;

    /* A key longer than the block is replaced by its hash, shorter is zero
     * padded. RFC 2104 section 2. */
    memset(k, 0, sizeof k);
    if (key_len > ELIPS_SHA256_BLOCK) {
        elips_sha256(k, key, key_len);
    } else if (key_len > 0) {
        memcpy(k, key, key_len);
    }

    for (size_t i = 0; i < ELIPS_SHA256_BLOCK; i++) pad[i] = k[i] ^ 0x36;
    elips_sha256_init(&s);
    elips_sha256_update(&s, pad, sizeof pad);
    elips_sha256_update(&s, msg, msg_len);
    elips_sha256_final(inner, &s);

    for (size_t i = 0; i < ELIPS_SHA256_BLOCK; i++) pad[i] = k[i] ^ 0x5c;
    elips_sha256_init(&s);
    elips_sha256_update(&s, pad, sizeof pad);
    elips_sha256_update(&s, inner, sizeof inner);
    elips_sha256_final(out, &s);

    /* The key schedule is derived from a secret in the one caller that
     * matters, so do not leave it on the stack. */
    memset(k, 0, sizeof k);
    memset(pad, 0, sizeof pad);
    memset(inner, 0, sizeof inner);
}

void elips_hkdf_extract(uint8_t prk[ELIPS_HMAC_BYTES],
                        const uint8_t *salt, size_t salt_len,
                        const uint8_t *ikm, size_t ikm_len)
{
    elips_hmac_sha256(prk, salt, salt_len, ikm, ikm_len);
}

int elips_hkdf_expand(uint8_t *out, size_t out_len,
                      const uint8_t prk[ELIPS_HMAC_BYTES],
                      const uint8_t *info, size_t info_len)
{
    if (out_len > 255u * ELIPS_HMAC_BYTES) return -1;

    uint8_t t[ELIPS_HMAC_BYTES];
    uint8_t block[ELIPS_HMAC_BYTES + 512 + 1];
    size_t done = 0;
    uint8_t counter = 1;

    if (info_len > 512) return -1;          /* no caller here comes close */

    while (done < out_len) {
        size_t n = 0;
        /* T(1) has no previous block; T(i) is prefixed by T(i-1). */
        if (counter > 1) { memcpy(block, t, sizeof t); n = sizeof t; }
        if (info_len) { memcpy(block + n, info, info_len); n += info_len; }
        block[n++] = counter;

        elips_hmac_sha256(t, prk, ELIPS_HMAC_BYTES, block, n);

        size_t take = out_len - done;
        if (take > sizeof t) take = sizeof t;
        memcpy(out + done, t, take);
        done += take;
        counter++;
    }

    memset(t, 0, sizeof t);
    memset(block, 0, sizeof block);
    return 0;
}
