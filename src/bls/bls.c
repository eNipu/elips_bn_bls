/*
 * BLS signatures, draft-irtf-cfrg-bls-signature. See include/elips/bls.h.
 *
 * Deliberately a thin layer. Everything below it -- hash-to-curve with the
 * registered RFC 9380 suite, the compressed encodings, subgroup checks, the
 * pairing -- is already pinned against published vectors and against py_ecc,
 * and was verified to agree with py_ecc's BLS byte for byte before this file
 * was written. What is genuinely new here is KeyGen, the aggregation paths,
 * and the input validation.
 *
 * ON SCALAR MULTIPLICATION. Secrets go through ep_mul and ep2_mul, the plain
 * constant-time fixed-window ladders, NOT through ep_mul_glv / ep2_mul_glv,
 * which are roughly 1.5x faster. That is not caution, it is a measurement:
 * under clang the GLV paths read max|t| = 32.97 -> 69.36 in dudect while the
 * plain ladders read 0.94 and 1.40, and issue #30 is open on why. Signing is
 * the one place in this library where a scalar is a long-term secret. When #30
 * closes, this decision should be revisited with a benchmark, not assumed.
 */
#include "elips/bls.h"

#include <string.h>

#include "elips/ec.h"
#include "elips/hash_to_curve.h"
#include "elips/pairing.h"
#include "elips/random.h"
#include "elips/serialize.h"
#include "elips/sha256.h"
#include "elips/sysrand.h"

#include "hkdf.h"

/* The registered ciphersuite strings, and not parameters.
 *
 * ONE SCHEME THROUGHOUT: proof of possession. The draft defines three, and the
 * first cut of this file mixed two of them -- it signed under the basic
 * scheme's NUL_ tag and aggregated under the proof-of-possession rules. The
 * cross-check against py_ecc caught it immediately, which is what it is for.
 * Mixing is not a constant to fix, it is an API that cannot work: a signature
 * made under one tag is not a signature under the other, so the aggregate
 * paths would have silently refused signatures this same file produced.
 *
 * Proof of possession is the one the ecosystem deploys, and it is the only one
 * of the three under which every function in elips/bls.h is meaningful.
 *
 * Note the tags differ by more than a suffix: a proof of possession is signed
 * under BLS_POP_, not BLS_SIG_, so a proof can never be replayed as a
 * signature over a message that happens to equal an encoded public key. */
static const char DST_SIG[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
static const char DST_POP[] = "BLS_POP_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
#define DST_SIG_LEN (sizeof DST_SIG - 1)
#define DST_POP_LEN (sizeof DST_POP - 1)

/* KeyGen's HKDF output length: ceil((1.5 * ceil(log2(r))) / 8) = 48 for the
 * 255-bit order of BLS12-381. Wider than the order on purpose, so reducing it
 * leaves no measurable bias. */
#define KEYGEN_L 48

/* ------------------------------------------------------------ scalars ---- */

/* 32 bytes big-endian -> limbs little-endian. Returns 0 if the value is a
 * valid secret key: strictly less than the order, and not zero. */
static int sk_to_limbs(limb_t *k, const uint8_t sk[ELIPS_BLS_SK_BYTES])
{
    memset(k, 0, ELIPS_ORDER_LIMBS * sizeof(limb_t));
    for (int i = 0; i < ELIPS_BLS_SK_BYTES; i++) {
        int weight = ELIPS_BLS_SK_BYTES - 1 - i;
        if (weight / 8 < ELIPS_ORDER_LIMBS)
            k[weight / 8] |= (limb_t)sk[i] << (8 * (weight % 8));
    }
    if (!elips_scalar_is_reduced(k)) return ELIPS_BLS_ERR_INVALID;

    limb_t acc = 0;
    for (int i = 0; i < ELIPS_ORDER_LIMBS; i++) acc |= k[i];
    if (acc == 0) return ELIPS_BLS_ERR_INVALID;   /* zero has no public key */
    return ELIPS_BLS_OK;
}

static void limbs_to_sk(uint8_t sk[ELIPS_BLS_SK_BYTES], const limb_t *k)
{
    memset(sk, 0, ELIPS_BLS_SK_BYTES);
    for (int i = 0; i < ELIPS_BLS_SK_BYTES; i++) {
        int weight = ELIPS_BLS_SK_BYTES - 1 - i;
        if (weight / 8 < ELIPS_ORDER_LIMBS)
            sk[i] = (uint8_t)(k[weight / 8] >> (8 * (weight % 8)));
    }
}

/* ------------------------------------------------------------- keygen ---- */

int elips_bls_keygen(uint8_t sk[ELIPS_BLS_SK_BYTES],
                     const uint8_t *ikm, size_t ikm_len)
{
    if (sk == NULL || ikm == NULL || ikm_len < 32) return ELIPS_BLS_ERR_INVALID;

    uint8_t salt[ELIPS_SHA256_DIGEST];
    uint8_t prk[ELIPS_HMAC_BYTES];
    uint8_t okm[KEYGEN_L];
    uint8_t ikm0[512];
    const uint8_t info[2] = { (uint8_t)(KEYGEN_L >> 8), (uint8_t)(KEYGEN_L & 0xff) };

    if (ikm_len > sizeof ikm0 - 1) return ELIPS_BLS_ERR_INVALID;
    memcpy(ikm0, ikm, ikm_len);
    ikm0[ikm_len] = 0x00;                       /* IKM || I2OSP(0, 1) */

    /* salt starts as the literal, then is rehashed at the top of every pass.
     * The loop can only repeat if the reduction lands on zero, which is a
     * 2^-255 event; it is written out rather than assumed away because the
     * draft specifies it. */
    static const char salt0[] = "BLS-SIG-KEYGEN-SALT-";
    elips_sha256(salt, salt0, sizeof salt0 - 1);

    int rc = ELIPS_BLS_ERR_INVALID;
    for (int attempt = 0; attempt < 8; attempt++) {
        elips_hkdf_extract(prk, salt, sizeof salt, ikm0, ikm_len + 1);
        if (elips_hkdf_expand(okm, sizeof okm, prk, info, sizeof info) != 0)
            break;

        /* OS2IP(okm) mod r. okm is 48 bytes big-endian, the order is 4 limbs. */
        limb_t wide[(KEYGEN_L + 7) / 8];
        memset(wide, 0, sizeof wide);
        for (int i = 0; i < KEYGEN_L; i++) {
            int weight = KEYGEN_L - 1 - i;
            wide[weight / 8] |= (limb_t)okm[i] << (8 * (weight % 8));
        }
        elips_scalar_reduce_wide(wide, (int)(sizeof wide / sizeof wide[0]));

        limb_t acc = 0;
        for (int i = 0; i < ELIPS_ORDER_LIMBS; i++) acc |= wide[i];
        if (acc != 0) { limbs_to_sk(sk, wide); rc = ELIPS_BLS_OK; }
        memset(wide, 0, sizeof wide);
        if (rc == ELIPS_BLS_OK) break;

        elips_sha256(salt, salt, sizeof salt);
    }

    memset(prk, 0, sizeof prk);
    memset(okm, 0, sizeof okm);
    memset(ikm0, 0, sizeof ikm0);
    return rc;
}

int elips_bls_keygen_random(uint8_t sk[ELIPS_BLS_SK_BYTES])
{
    uint8_t ikm[32];
    if (elips_random_bytes(ikm, sizeof ikm) != 0) return ELIPS_BLS_ERR_RANDOM;
    int rc = elips_bls_keygen(sk, ikm, sizeof ikm);
    memset(ikm, 0, sizeof ikm);
    return rc;
}

/* --------------------------------------------------------------- keys ---- */

int elips_bls_sk_to_pk(uint8_t pk[ELIPS_BLS_PK_BYTES],
                       const uint8_t sk[ELIPS_BLS_SK_BYTES])
{
    if (pk == NULL || sk == NULL) return ELIPS_BLS_ERR_INVALID;
    limb_t k[ELIPS_ORDER_LIMBS];
    int rc = sk_to_limbs(k, sk);
    if (rc != ELIPS_BLS_OK) return rc;

    ep_t G1, P;
    ep_generator(&G1);
    ep_mul(&P, &G1, k, ELIPS_ORDER_BITS);       /* see the note at the top */
    memset(k, 0, sizeof k);

    if (ep_is_infinity(&P)) return ELIPS_BLS_ERR_BAD_KEY;
    ep_write_compressed(pk, &P);
    return ELIPS_BLS_OK;
}

static int read_pk(ep_t *out, const uint8_t pk[ELIPS_BLS_PK_BYTES])
{
    if (ep_read_compressed(out, pk) != ELIPS_SER_OK) return ELIPS_BLS_ERR_BAD_KEY;
    /* KeyValidate: the identity is a valid encoding and a useless key, and
     * accepting it lets a signer claim membership with no key at all. */
    if (ep_is_infinity(out)) return ELIPS_BLS_ERR_BAD_KEY;
    return ELIPS_BLS_OK;
}

static int read_sig(ep2_t *out, const uint8_t sig[ELIPS_BLS_SIG_BYTES])
{
    if (ep2_read_compressed(out, sig) != ELIPS_SER_OK) return ELIPS_BLS_ERR_BAD_SIG;
    return ELIPS_BLS_OK;
}

int elips_bls_pk_validate(const uint8_t pk[ELIPS_BLS_PK_BYTES])
{
    if (pk == NULL) return ELIPS_BLS_ERR_INVALID;
    ep_t P;
    return read_pk(&P, pk);
}

/* --------------------------------------------------- sign and verify ---- */

static int core_sign(uint8_t sig[ELIPS_BLS_SIG_BYTES],
                     const uint8_t sk[ELIPS_BLS_SK_BYTES],
                     const uint8_t *msg, size_t msg_len,
                     const char *dst, size_t dst_len)
{
    limb_t k[ELIPS_ORDER_LIMBS];
    int rc = sk_to_limbs(k, sk);
    if (rc != ELIPS_BLS_OK) return rc;

    ep2_t H, S;
    if (elips_hash_to_g2(&H, msg, msg_len, (const uint8_t *)dst, dst_len) != 0) {
        memset(k, 0, sizeof k);
        return ELIPS_BLS_ERR_INVALID;
    }
    ep2_mul(&S, &H, k, ELIPS_ORDER_BITS);       /* see the note at the top */
    memset(k, 0, sizeof k);

    ep2_write_compressed(sig, &S);
    return ELIPS_BLS_OK;
}

/* e(-G1, sig) * prod_i e(pk_i, H(m_i)) == 1.
 *
 * Written as one multi-pairing rather than two separate pairings so it costs
 * one final exponentiation instead of n+1. Negating the generator rather than
 * the signature keeps the caller's bytes untouched.
 *
 * The final exponentiation returns e^3 on BLS12 (see pairing.h). That does not
 * affect this: the group has prime order r with gcd(3, r) = 1, so a cube is 1
 * exactly when the value is. */
static int pairing_check(const ep_t *pks, const ep2_t *hs, size_t n,
                         const ep2_t *sig)
{
    enum { MAXN = 64 };
    if (n == 0 || n > MAXN - 1) return ELIPS_BLS_ERR_INVALID;

    ep_t  P[MAXN];
    ep2_t Q[MAXN];

    ep_t G1;
    ep_generator(&G1);
    ep_neg(&P[0], &G1);
    ep2_copy(&Q[0], sig);
    for (size_t i = 0; i < n; i++) {
        ep_copy(&P[i + 1], &pks[i]);
        ep2_copy(&Q[i + 1], &hs[i]);
    }

    fp12_t f;
    if (!elips_pairing_multi(f, P, Q, n + 1)) return ELIPS_BLS_ERR_VERIFY;

    fp12_t one;
    fp12_set_one(one);
    return fp12_eq(f, one) ? ELIPS_BLS_OK : ELIPS_BLS_ERR_VERIFY;
}

int elips_bls_sign(uint8_t sig[ELIPS_BLS_SIG_BYTES],
                   const uint8_t sk[ELIPS_BLS_SK_BYTES],
                   const uint8_t *msg, size_t msg_len)
{
    if (sig == NULL || sk == NULL || (msg == NULL && msg_len)) return ELIPS_BLS_ERR_INVALID;
    return core_sign(sig, sk, msg, msg_len, DST_SIG, DST_SIG_LEN);
}

static int core_verify(const uint8_t pk[ELIPS_BLS_PK_BYTES],
                       const uint8_t *msg, size_t msg_len,
                       const uint8_t sig[ELIPS_BLS_SIG_BYTES],
                       const char *dst, size_t dst_len)
{
    if (pk == NULL || sig == NULL || (msg == NULL && msg_len)) return ELIPS_BLS_ERR_INVALID;

    ep_t P;
    int rc = read_pk(&P, pk);
    if (rc != ELIPS_BLS_OK) return rc;

    ep2_t S;
    rc = read_sig(&S, sig);
    if (rc != ELIPS_BLS_OK) return rc;

    ep2_t H;
    if (elips_hash_to_g2(&H, msg, msg_len, (const uint8_t *)dst, dst_len) != 0)
        return ELIPS_BLS_ERR_INVALID;

    return pairing_check(&P, &H, 1, &S);
}

int elips_bls_verify(const uint8_t pk[ELIPS_BLS_PK_BYTES],
                     const uint8_t *msg, size_t msg_len,
                     const uint8_t sig[ELIPS_BLS_SIG_BYTES])
{
    return core_verify(pk, msg, msg_len, sig, DST_SIG, DST_SIG_LEN);
}

/* -------------------------------------------------------- aggregation ---- */

int elips_bls_aggregate(uint8_t out[ELIPS_BLS_SIG_BYTES],
                        const uint8_t *sigs, size_t n)
{
    if (out == NULL || sigs == NULL || n == 0) return ELIPS_BLS_ERR_INVALID;

    ep2_t acc, S;
    int rc = read_sig(&acc, sigs);
    if (rc != ELIPS_BLS_OK) return rc;

    for (size_t i = 1; i < n; i++) {
        rc = read_sig(&S, sigs + i * ELIPS_BLS_SIG_BYTES);
        if (rc != ELIPS_BLS_OK) return rc;
        ep2_add(&acc, &acc, &S);
    }
    ep2_write_compressed(out, &acc);
    return ELIPS_BLS_OK;
}

int elips_bls_aggregate_verify(const uint8_t *pks, size_t n,
                               const uint8_t *const *msgs, const size_t *msg_lens,
                               const uint8_t sig[ELIPS_BLS_SIG_BYTES])
{
    enum { MAXN = 63 };
    if (pks == NULL || msgs == NULL || msg_lens == NULL || sig == NULL)
        return ELIPS_BLS_ERR_INVALID;
    if (n == 0 || n > MAXN) return ELIPS_BLS_ERR_INVALID;

    /* Distinct messages, or the basic scheme is broken. Comparing the messages
     * themselves rather than a hash of them: n is small, and a hash comparison
     * would trade a real check for a probabilistic one to save nothing. */
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            if (msg_lens[i] == msg_lens[j] &&
                memcmp(msgs[i], msgs[j], msg_lens[i]) == 0)
                return ELIPS_BLS_ERR_DUP_MESSAGE;

    ep_t  P[MAXN];
    ep2_t H[MAXN];
    for (size_t i = 0; i < n; i++) {
        int rc = read_pk(&P[i], pks + i * ELIPS_BLS_PK_BYTES);
        if (rc != ELIPS_BLS_OK) return rc;
        if (elips_hash_to_g2(&H[i], msgs[i], msg_lens[i],
                             (const uint8_t *)DST_SIG, DST_SIG_LEN) != 0)
            return ELIPS_BLS_ERR_INVALID;
    }

    ep2_t S;
    int rc = read_sig(&S, sig);
    if (rc != ELIPS_BLS_OK) return rc;

    return pairing_check(P, H, n, &S);
}

int elips_bls_fast_aggregate_verify(const uint8_t *pks,
                                    const uint8_t *pops, size_t n,
                                    const uint8_t *msg, size_t msg_len,
                                    const uint8_t sig[ELIPS_BLS_SIG_BYTES])
{
    if (pks == NULL || pops == NULL || sig == NULL) return ELIPS_BLS_ERR_INVALID;
    if (n == 0) return ELIPS_BLS_ERR_INVALID;

    /* The proofs are checked FIRST and unconditionally. Without them, summing
     * public keys over a shared message is forgeable; see the header. */
    ep_t acc, P;
    for (size_t i = 0; i < n; i++) {
        const uint8_t *pk = pks + i * ELIPS_BLS_PK_BYTES;
        int rc = elips_bls_pop_verify(pk, pops + i * ELIPS_BLS_POP_BYTES);
        if (rc != ELIPS_BLS_OK) return rc;

        rc = read_pk(&P, pk);
        if (rc != ELIPS_BLS_OK) return rc;
        if (i == 0) ep_copy(&acc, &P);
        else        ep_add(&acc, &acc, &P);
    }
    if (ep_is_infinity(&acc)) return ELIPS_BLS_ERR_BAD_KEY;

    ep2_t S, H;
    int rc = read_sig(&S, sig);
    if (rc != ELIPS_BLS_OK) return rc;
    if (elips_hash_to_g2(&H, msg, msg_len,
                         (const uint8_t *)DST_SIG, DST_SIG_LEN) != 0)
        return ELIPS_BLS_ERR_INVALID;

    return pairing_check(&acc, &H, 1, &S);
}

/* ------------------------------------------------ proof of possession ---- */

int elips_bls_pop_prove(uint8_t proof[ELIPS_BLS_POP_BYTES],
                        const uint8_t sk[ELIPS_BLS_SK_BYTES])
{
    if (proof == NULL || sk == NULL) return ELIPS_BLS_ERR_INVALID;
    uint8_t pk[ELIPS_BLS_PK_BYTES];
    int rc = elips_bls_sk_to_pk(pk, sk);
    if (rc != ELIPS_BLS_OK) return rc;
    /* Signed under BLS_POP_, a different tag from the BLS_SIG_ used for
     * messages, so a proof can never be replayed as a signature. */
    return core_sign(proof, sk, pk, sizeof pk, DST_POP, DST_POP_LEN);
}

int elips_bls_pop_verify(const uint8_t pk[ELIPS_BLS_PK_BYTES],
                         const uint8_t proof[ELIPS_BLS_POP_BYTES])
{
    if (pk == NULL || proof == NULL) return ELIPS_BLS_ERR_INVALID;
    return core_verify(pk, pk, ELIPS_BLS_PK_BYTES, proof, DST_POP, DST_POP_LEN);
}

const char *elips_bls_strerror(int code)
{
    switch (code) {
    case ELIPS_BLS_OK:                return "ok";
    case ELIPS_BLS_ERR_INVALID:       return "malformed or out-of-range input";
    case ELIPS_BLS_ERR_BAD_KEY:       return "public key is infinity or not in G1";
    case ELIPS_BLS_ERR_BAD_SIG:       return "signature does not decode or is not in G2";
    case ELIPS_BLS_ERR_VERIFY:        return "signature does not verify";
    case ELIPS_BLS_ERR_RANDOM:        return "the system entropy source failed";
    case ELIPS_BLS_ERR_DUP_MESSAGE:   return "aggregate verify needs distinct messages";
    default:                          return "unknown error";
    }
}
