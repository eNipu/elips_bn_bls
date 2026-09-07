/*
 * SHA-256, FIPS 180-4. See include/elips/sha256.h.
 *
 * Straightforward and unoptimised: it hashes at most a few hundred bytes per
 * hash-to-curve call, against a scalar multiplication that costs milliseconds,
 * so there is nothing here worth tuning. Constant time by construction -- no
 * branch or index depends on the message bytes.
 */
#include "elips/sha256.h"
#include <string.h>

static const uint32_t K[64] = {
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,
    0x923f82a4u,0xab1c5ed5u,0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
    0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,0xe49b69c1u,0xefbe4786u,
    0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,
    0x06ca6351u,0x14292967u,0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
    0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,0xa2bfe8a1u,0xa81a664bu,
    0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,
    0x5b9cca4fu,0x682e6ff3u,0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
    0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

#define ROR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
#define S0(x) (ROR(x, 2)  ^ ROR(x, 13) ^ ROR(x, 22))
#define S1(x) (ROR(x, 6)  ^ ROR(x, 11) ^ ROR(x, 25))
#define s0(x) (ROR(x, 7)  ^ ROR(x, 18) ^ ((x) >> 3))
#define s1(x) (ROR(x, 17) ^ ROR(x, 19) ^ ((x) >> 10))

static void compress(uint32_t h[8], const uint8_t block[64])
{
    uint32_t w[64];
    for (int i = 0; i < 16; i++)
        w[i] = ((uint32_t)block[4*i] << 24) | ((uint32_t)block[4*i+1] << 16) |
               ((uint32_t)block[4*i+2] << 8) | (uint32_t)block[4*i+3];
    for (int i = 16; i < 64; i++)
        w[i] = s1(w[i-2]) + w[i-7] + s0(w[i-15]) + w[i-16];

    uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
    uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];

    for (int i = 0; i < 64; i++) {
        uint32_t t1 = hh + S1(e) + ((e & f) ^ (~e & g)) + K[i] + w[i];
        uint32_t t2 = S0(a) + ((a & b) ^ (a & c) ^ (b & c));
        hh = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d;
    h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}

void elips_sha256_init(elips_sha256_t *s)
{
    static const uint32_t iv[8] = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
    };
    memcpy(s->h, iv, sizeof iv);
    s->len = 0;
    s->buf_len = 0;
}

void elips_sha256_update(elips_sha256_t *s, const void *data, size_t len)
{
    const uint8_t *p = (const uint8_t *)data;
    s->len += len;

    if (s->buf_len) {
        size_t take = ELIPS_SHA256_BLOCK - s->buf_len;
        if (take > len) take = len;
        memcpy(s->buf + s->buf_len, p, take);
        s->buf_len += take; p += take; len -= take;
        if (s->buf_len == ELIPS_SHA256_BLOCK) { compress(s->h, s->buf); s->buf_len = 0; }
    }
    while (len >= ELIPS_SHA256_BLOCK) {
        compress(s->h, p);
        p += ELIPS_SHA256_BLOCK; len -= ELIPS_SHA256_BLOCK;
    }
    if (len) { memcpy(s->buf, p, len); s->buf_len = len; }
}

void elips_sha256_final(uint8_t out[ELIPS_SHA256_DIGEST], elips_sha256_t *s)
{
    uint64_t bits = s->len * 8;
    uint8_t pad[ELIPS_SHA256_BLOCK * 2];
    size_t n = ELIPS_SHA256_BLOCK - (size_t)(s->len % ELIPS_SHA256_BLOCK);
    if (n < 9) n += ELIPS_SHA256_BLOCK;      /* room for 0x80 and the length */

    memset(pad, 0, n);
    pad[0] = 0x80;
    for (int i = 0; i < 8; i++)
        pad[n - 1 - i] = (uint8_t)(bits >> (8 * i));
    elips_sha256_update(s, pad, n);

    for (int i = 0; i < 8; i++) {
        out[4*i]     = (uint8_t)(s->h[i] >> 24);
        out[4*i + 1] = (uint8_t)(s->h[i] >> 16);
        out[4*i + 2] = (uint8_t)(s->h[i] >> 8);
        out[4*i + 3] = (uint8_t)(s->h[i]);
    }
}

void elips_sha256(uint8_t out[ELIPS_SHA256_DIGEST], const void *data, size_t len)
{
    elips_sha256_t s;
    elips_sha256_init(&s);
    elips_sha256_update(&s, data, len);
    elips_sha256_final(out, &s);
}
