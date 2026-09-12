/*
 * elips_bls, without Python.
 *
 * tools/verify/crosscheck_bls_pyecc.py is the stronger test: it compares every
 * operation against an independent implementation on randomized inputs, both
 * directions. This one exists because that test needs py_ecc installed, and a
 * suite that cannot run without a Python package is a suite that stops being
 * run. What is here:
 *
 *   - the RFC 5869 vectors for the HKDF underneath KeyGen
 *   - one end-to-end vector pinned against py_ecc, so a regression in signing
 *     is caught even with no network and no pip
 *   - every rejection path, each of which FAILS if its check is removed from
 *     src/bls/bls.c. That is the rule this repository already applies to its
 *     vector runners, and a verifier that cannot reject is not a verifier.
 */
#include <stdio.h>
#include <string.h>

#include "elips/bls.h"
#include "hkdf.h"

static int failures, checks;

static void ok(int cond, const char *what)
{
    checks++;
    if (!cond) { failures++; printf("  FAIL %s\n", what); }
}

static void is(int got, int want, const char *what)
{
    checks++;
    if (got != want) {
        failures++;
        printf("  FAIL %s: got %d (%s), want %d (%s)\n",
               what, got, elips_bls_strerror(got), want, elips_bls_strerror(want));
    }
}

static int unhex(uint8_t *o, const char *h)
{ int n = 0; while (*h) { unsigned v; sscanf(h, "%2x", &v); o[n++] = (uint8_t)v; h += 2; } return n; }

/* Pinned against py_ecc's G2ProofOfPossession, ikm = 32 bytes of 0x01,
 * message "hello elips". Regenerate with tools/verify/crosscheck_bls_pyecc.py
 * if the scheme ever changes; do not edit by hand. */
static const uint8_t VEC_SK[32] = {
    0x14, 0x4b, 0x27, 0x82, 0x8e, 0x30, 0x5a, 0x2d, 0x67, 0xfc, 0x7f, 0x4e, 0xea, 0x6d, 0xe7, 0x06,
    0xb4, 0x05, 0xcd, 0xd1, 0xab, 0x8a, 0xd2, 0xda, 0xec, 0x04, 0x6c, 0xcd, 0xee, 0xec, 0x8b, 0x79,
};
static const uint8_t VEC_PK[48] = {
    0x95, 0xa2, 0x54, 0x50, 0x1b, 0x77, 0x33, 0x23, 0x9e, 0xd3, 0xce, 0xc4, 0xd5, 0x67, 0x37, 0x97,
    0x7b, 0xd0, 0x9e, 0xde, 0x88, 0x1d, 0x8a, 0x23, 0x45, 0x60, 0xe8, 0x3e, 0x55, 0x25, 0x01, 0x7a,
    0xdd, 0x3b, 0x1d, 0xcc, 0x3e, 0xab, 0xfb, 0x85, 0xe1, 0x2a, 0x41, 0x31, 0xb1, 0x9c, 0x25, 0x3b,
};
static const uint8_t VEC_SIG[96] = {
    0x8b, 0xd3, 0x03, 0x37, 0x71, 0xb8, 0xac, 0x8c, 0x47, 0xea, 0x93, 0x75, 0x2c, 0xb1, 0xc1, 0x74,
    0x83, 0x2c, 0x11, 0x89, 0x03, 0x8a, 0xfe, 0xf4, 0x42, 0x3b, 0x71, 0x93, 0x2c, 0xd1, 0x42, 0x6d,
    0x49, 0xb5, 0x24, 0xd3, 0xc9, 0x4c, 0xcc, 0xc7, 0x26, 0xb8, 0x00, 0x41, 0x34, 0x06, 0x89, 0xe5,
    0x04, 0x26, 0x1d, 0xe4, 0xbb, 0xc9, 0x59, 0xc3, 0x1f, 0x75, 0xfb, 0x60, 0x38, 0xe0, 0xd5, 0x05,
    0x91, 0xc3, 0xb5, 0x2c, 0x18, 0xeb, 0x53, 0x1c, 0x8c, 0xd4, 0xb3, 0x53, 0x1d, 0x49, 0x42, 0x8c,
    0x53, 0xe0, 0x71, 0x74, 0x34, 0x7c, 0x73, 0xbb, 0xea, 0x29, 0xca, 0x48, 0x6a, 0x6f, 0x62, 0x06,
};
static const uint8_t VEC_POP[96] = {
    0x84, 0x6a, 0xa1, 0x2a, 0x44, 0x02, 0xeb, 0x67, 0xcb, 0x92, 0xa4, 0x97, 0xe0, 0x71, 0x6d, 0xb5,
    0x73, 0xc8, 0x17, 0xa4, 0x16, 0x37, 0x83, 0x15, 0x3f, 0x0d, 0xdc, 0xa4, 0x75, 0xf4, 0x87, 0x02,
    0x00, 0x04, 0x9d, 0x8e, 0x9e, 0xd3, 0x50, 0x87, 0xc7, 0x86, 0x05, 0x9c, 0x1f, 0x26, 0xfc, 0x9d,
    0x0d, 0x39, 0xe3, 0x09, 0x8f, 0x1b, 0xae, 0x07, 0x4c, 0x06, 0x2f, 0x84, 0xf2, 0x43, 0x53, 0x21,
    0x06, 0x66, 0xbd, 0x58, 0xc0, 0xd9, 0xbe, 0x3f, 0xf7, 0x6b, 0xa9, 0xdd, 0x9c, 0xe9, 0x05, 0xc5,
    0xb6, 0x02, 0xa1, 0x2e, 0x78, 0xa0, 0x43, 0x50, 0x27, 0x5f, 0xaa, 0xcc, 0xe8, 0xb7, 0x13, 0x7d,
};

/* --- RFC 5869, the SHA-256 cases ----------------------------------------- */
static void test_hkdf(void)
{
    struct { const char *ikm, *salt, *info, *okm; int L; } v[] = {
      { "0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b", "000102030405060708090a0b0c",
        "f0f1f2f3f4f5f6f7f8f9",
        "3cb25f25faacd57a90434f64d0362f2a2d2d0a90cf1a5a4c5db02d56ecc4c5bf34007208d5b887185865", 42 },
      { "0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b", "", "",
        "8da4e775a563c18f715f802a063c5a31b8a11f5c5ee1879ec3454e5f3c738d2d9d201395faa4b61a96c8", 42 },
    };
    for (unsigned i = 0; i < sizeof v / sizeof *v; i++) {
        uint8_t ikm[128], salt[128], info[128], want[128], prk[32], got[128];
        int ni = unhex(ikm, v[i].ikm), ns = unhex(salt, v[i].salt);
        int nf = unhex(info, v[i].info);
        unhex(want, v[i].okm);
        elips_hkdf_extract(prk, salt, (size_t)ns, ikm, (size_t)ni);
        elips_hkdf_expand(got, (size_t)v[i].L, prk, info, (size_t)nf);
        ok(memcmp(got, want, (size_t)v[i].L) == 0, "RFC 5869 HKDF-SHA256 vector");
    }
}

int main(void)
{
    printf("elips_bls\n");
    test_hkdf();

    /* --- the pinned vector ------------------------------------------------ */
    uint8_t ikm[32]; memset(ikm, 0x01, sizeof ikm);
    uint8_t sk[32], pk[48], sig[96], pop[96];
    const char *m = "hello elips";
    const size_t ml = 11;

    is(elips_bls_keygen(sk, ikm, sizeof ikm), ELIPS_BLS_OK, "keygen");
    ok(memcmp(sk, VEC_SK, sizeof sk) == 0, "keygen matches the pinned secret key");
    is(elips_bls_sk_to_pk(pk, sk), ELIPS_BLS_OK, "sk_to_pk");
    ok(memcmp(pk, VEC_PK, sizeof pk) == 0, "public key matches the pinned value");
    is(elips_bls_sign(sig, sk, (const uint8_t *)m, ml), ELIPS_BLS_OK, "sign");
    ok(memcmp(sig, VEC_SIG, sizeof sig) == 0, "signature matches the pinned value");
    is(elips_bls_verify(pk, (const uint8_t *)m, ml, sig), ELIPS_BLS_OK, "verify");
    is(elips_bls_pop_prove(pop, sk), ELIPS_BLS_OK, "pop_prove");
    ok(memcmp(pop, VEC_POP, sizeof pop) == 0, "proof matches the pinned value");
    is(elips_bls_pop_verify(pk, pop), ELIPS_BLS_OK, "pop_verify");

    /* KeyGen is deterministic, which is what makes seed phrases work. */
    uint8_t sk2[32];
    is(elips_bls_keygen(sk2, ikm, sizeof ikm), ELIPS_BLS_OK, "keygen again");
    ok(memcmp(sk, sk2, sizeof sk) == 0, "keygen is deterministic in its input");

    /* --- rejection. Each of these fails if its check is removed. ---------- */
    is(elips_bls_verify(pk, (const uint8_t *)"hello elipt", ml, sig),
       ELIPS_BLS_ERR_VERIFY, "tampered message is rejected");

    uint8_t bad[96];
    memcpy(bad, sig, sizeof bad); bad[95] ^= 1;
    ok(elips_bls_verify(pk, (const uint8_t *)m, ml, bad) != ELIPS_BLS_OK,
       "mangled signature is rejected");

    uint8_t badpk[48];
    memcpy(badpk, pk, sizeof badpk); badpk[47] ^= 1;
    ok(elips_bls_verify(badpk, (const uint8_t *)m, ml, sig) != ELIPS_BLS_OK,
       "wrong public key is rejected");

    /* The compressed encoding of infinity: well formed, and useless as a key.
     * Accepting it would let a signer claim membership with no key at all. */
    uint8_t infpk[48]; memset(infpk, 0, sizeof infpk); infpk[0] = 0xc0;
    is(elips_bls_pk_validate(infpk), ELIPS_BLS_ERR_BAD_KEY,
       "the identity is rejected as a public key");
    is(elips_bls_verify(infpk, (const uint8_t *)m, ml, sig), ELIPS_BLS_ERR_BAD_KEY,
       "verify rejects the identity public key");

    /* Secret keys out of range. Zero has no public key; the order itself and
     * anything above it is not a scalar. */
    uint8_t zero[32]; memset(zero, 0, sizeof zero);
    is(elips_bls_sk_to_pk(pk, zero), ELIPS_BLS_ERR_INVALID, "zero secret key is rejected");
    uint8_t maxsk[32]; memset(maxsk, 0xff, sizeof maxsk);
    is(elips_bls_sk_to_pk(pk, maxsk), ELIPS_BLS_ERR_INVALID,
       "an out-of-range secret key is rejected");
    is(elips_bls_sk_to_pk(pk, sk), ELIPS_BLS_OK, "pk restored");

    /* KeyGen requires at least 32 bytes of input keying material. */
    is(elips_bls_keygen(sk2, ikm, 31), ELIPS_BLS_ERR_INVALID,
       "short input keying material is rejected");

    /* --- aggregation ------------------------------------------------------ */
    enum { N = 4 };
    uint8_t sks[N][32], pks[N][48], sigs[N][96], pops[N][96];
    const char *msgs_s[N] = { "alpha", "beta", "gamma", "delta" };
    const uint8_t *msgs[N];
    size_t lens[N];
    for (int i = 0; i < N; i++) {
        uint8_t seed[32]; memset(seed, (uint8_t)(0x10 + i), sizeof seed);
        is(elips_bls_keygen(sks[i], seed, sizeof seed), ELIPS_BLS_OK, "keygen i");
        is(elips_bls_sk_to_pk(pks[i], sks[i]), ELIPS_BLS_OK, "sk_to_pk i");
        msgs[i] = (const uint8_t *)msgs_s[i];
        lens[i] = strlen(msgs_s[i]);
        is(elips_bls_sign(sigs[i], sks[i], msgs[i], lens[i]), ELIPS_BLS_OK, "sign i");
        is(elips_bls_pop_prove(pops[i], sks[i]), ELIPS_BLS_OK, "pop i");
    }

    uint8_t agg[96];
    is(elips_bls_aggregate(agg, (const uint8_t *)sigs, N), ELIPS_BLS_OK, "aggregate");
    is(elips_bls_aggregate_verify((const uint8_t *)pks, N, msgs, lens, agg),
       ELIPS_BLS_OK, "aggregate_verify");

    /* The aggregate is the same size as one signature. That is the point. */
    ok(sizeof agg == ELIPS_BLS_SIG_BYTES, "the aggregate of N is still 96 bytes");

    /* Swap two messages: the same keys, the same signatures, wrong pairing. */
    const uint8_t *swapped[N]; size_t swlens[N];
    for (int i = 0; i < N; i++) { swapped[i] = msgs[i]; swlens[i] = lens[i]; }
    swapped[0] = msgs[1]; swlens[0] = lens[1];
    swapped[1] = msgs[0]; swlens[1] = lens[0];
    is(elips_bls_aggregate_verify((const uint8_t *)pks, N, swapped, swlens, agg),
       ELIPS_BLS_ERR_VERIFY, "aggregate_verify rejects reordered messages");

    /* Repeated messages without proofs of possession are forgeable. */
    const uint8_t *dup[N]; size_t dlens[N];
    for (int i = 0; i < N; i++) { dup[i] = msgs[0]; dlens[i] = lens[0]; }
    is(elips_bls_aggregate_verify((const uint8_t *)pks, N, dup, dlens, agg),
       ELIPS_BLS_ERR_DUP_MESSAGE, "aggregate_verify refuses a repeated message");

    is(elips_bls_aggregate_verify((const uint8_t *)pks, 0, msgs, lens, agg),
       ELIPS_BLS_ERR_INVALID, "aggregate_verify refuses an empty set");
    is(elips_bls_aggregate(agg, (const uint8_t *)sigs, 0),
       ELIPS_BLS_ERR_INVALID, "aggregate refuses an empty set");

    /* --- the shared-message path, which demands the proofs ---------------- */
    const char *shared = "one message, many signers";
    size_t sl = strlen(shared);
    uint8_t ssigs[N][96], sagg[96];
    for (int i = 0; i < N; i++)
        is(elips_bls_sign(ssigs[i], sks[i], (const uint8_t *)shared, sl),
           ELIPS_BLS_OK, "sign shared");
    is(elips_bls_aggregate(sagg, (const uint8_t *)ssigs, N), ELIPS_BLS_OK, "aggregate shared");
    is(elips_bls_fast_aggregate_verify((const uint8_t *)pks, (const uint8_t *)pops,
                                       N, (const uint8_t *)shared, sl, sagg),
       ELIPS_BLS_OK, "fast_aggregate_verify");

    /* One bad proof stops the whole thing, which is why it takes them. */
    uint8_t badpops[N][96];
    memcpy(badpops, pops, sizeof badpops);
    badpops[2][95] ^= 1;
    ok(elips_bls_fast_aggregate_verify((const uint8_t *)pks, (const uint8_t *)badpops,
                                       N, (const uint8_t *)shared, sl, sagg) != ELIPS_BLS_OK,
       "fast_aggregate_verify rejects a bad proof of possession");

    /* A proof only proves possession of its own key. */
    ok(elips_bls_pop_verify(pks[0], pops[1]) != ELIPS_BLS_OK,
       "a proof does not verify under another key");

    printf("  %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
