/*
 * BLS signatures, draft-irtf-cfrg-bls-signature, on BLS12-381.
 *
 * This is the API to reach for first. Everything crosses it as bytes: no curve
 * point, no field element, no limb, no Montgomery form, no domain separation
 * tag. If you are verifying a signature, you should not have to know what a
 * pairing is, and nothing here asks you to.
 *
 *     uint8_t sk[ELIPS_BLS_SK_BYTES], pk[ELIPS_BLS_PK_BYTES];
 *     uint8_t sig[ELIPS_BLS_SIG_BYTES];
 *
 *     elips_bls_keygen(sk, ikm, sizeof ikm);
 *     elips_bls_sk_to_pk(pk, sk);
 *     elips_bls_sign(sig, sk, msg, msg_len);
 *     if (elips_bls_verify(pk, msg, msg_len, sig) == ELIPS_BLS_OK) { ... }
 *
 * VARIANT: minimal-pubkey-size. Public keys are 48-byte compressed G1 points
 * and signatures are 96-byte compressed G2 points. This is the variant
 * Ethereum uses, and the one py_ecc implements, which is what lets the test
 * suite check every operation here against an independent implementation on
 * every build.
 *
 * SCHEME: proof of possession, the one the ecosystem deploys. The draft
 * defines three and they are not interchangeable, so this implements exactly
 * one rather than letting a caller combine them by accident.
 *
 *   messages             BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_
 *   proofs of possession BLS_POP_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_
 *
 * Neither is a parameter. A caller who can pass the wrong domain separation
 * tag eventually will, and a wrong tag produces signatures that look perfectly
 * correct locally and verify against nothing else in the world.
 *
 * NOTHING HERE RETURNS A BOOLEAN. Every function returns ELIPS_BLS_OK, which
 * is zero, or a negative error code. `if (elips_bls_verify(...))` is therefore
 * true on FAILURE, which is the safe way round: the mistake is loud rather
 * than silent. Compare against ELIPS_BLS_OK explicitly.
 *
 * WHAT THIS IS NOT: audited. It implements the draft and agrees with an
 * independent implementation on every build, and that is a different claim
 * from having been reviewed by a cryptographer.
 */
#ifndef ELIPS_BLS_H
#define ELIPS_BLS_H

#include <stddef.h>
#include <stdint.h>

#if !defined(ELIPS_CURVE_BLS12_381)
#  error "elips/bls.h is BLS12-381 only: the other curves have no registered ciphersuite."
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ELIPS_BLS_SK_BYTES   32
#define ELIPS_BLS_PK_BYTES   48
#define ELIPS_BLS_SIG_BYTES  96
#define ELIPS_BLS_POP_BYTES  96

#define ELIPS_BLS_OK                 0
#define ELIPS_BLS_ERR_INVALID      (-1)  /* a malformed or out-of-range input */
#define ELIPS_BLS_ERR_BAD_KEY      (-2)  /* public key is infinity, or not in G1 */
#define ELIPS_BLS_ERR_BAD_SIG      (-3)  /* signature does not decode, or not in G2 */
#define ELIPS_BLS_ERR_VERIFY       (-4)  /* decoded fine, and does not verify */
#define ELIPS_BLS_ERR_RANDOM       (-5)  /* the system entropy source failed */
#define ELIPS_BLS_ERR_DUP_MESSAGE  (-6)  /* AggregateVerify needs distinct messages */

/* --- keys ------------------------------------------------------------------
 *
 * KeyGen derives a secret key from input keying material, per the draft: HKDF
 * over SHA-256, reduced modulo the group order. It is deterministic, so the
 * same ikm always gives the same key, which is what makes seed phrases and
 * test vectors work.
 *
 * ikm must be at least 32 bytes; the draft requires it and shorter input is
 * rejected rather than stretched. Pass elips_bls_keygen_random to have the
 * system entropy source supply it. */
int elips_bls_keygen(uint8_t sk[ELIPS_BLS_SK_BYTES],
                     const uint8_t *ikm, size_t ikm_len);
int elips_bls_keygen_random(uint8_t sk[ELIPS_BLS_SK_BYTES]);

int elips_bls_sk_to_pk(uint8_t pk[ELIPS_BLS_PK_BYTES],
                       const uint8_t sk[ELIPS_BLS_SK_BYTES]);

/* Is this 48-byte string a usable public key? Decodes, rejects infinity, and
 * checks G1 membership. Call it on anything that arrived from the network
 * before you store it; verify calls it for you. */
int elips_bls_pk_validate(const uint8_t pk[ELIPS_BLS_PK_BYTES]);

/* --- sign and verify, one key, one message -------------------------------- */

int elips_bls_sign(uint8_t sig[ELIPS_BLS_SIG_BYTES],
                   const uint8_t sk[ELIPS_BLS_SK_BYTES],
                   const uint8_t *msg, size_t msg_len);

int elips_bls_verify(const uint8_t pk[ELIPS_BLS_PK_BYTES],
                     const uint8_t *msg, size_t msg_len,
                     const uint8_t sig[ELIPS_BLS_SIG_BYTES]);

/* --- aggregation -----------------------------------------------------------
 *
 * The property BLS is worth having: n signatures combine into ONE signature of
 * the same 96 bytes, and it verifies against the n public keys at once.
 *
 * sigs is n consecutive 96-byte signatures, pks is n consecutive 48-byte keys,
 * and so on: plain arrays, no structs. */
int elips_bls_aggregate(uint8_t out[ELIPS_BLS_SIG_BYTES],
                        const uint8_t *sigs, size_t n);

/* n keys, n DISTINCT messages, one signature.
 *
 * The messages must differ. With the basic scheme, repeating a message lets an
 * attacker move a signature between signers, so a repeat is rejected here
 * rather than left as a footnote. Same message for everyone is the other
 * function. */
int elips_bls_aggregate_verify(const uint8_t *pks, size_t n,
                               const uint8_t *const *msgs, const size_t *msg_lens,
                               const uint8_t sig[ELIPS_BLS_SIG_BYTES]);

/* n keys, ONE shared message, one signature.
 *
 * This is the one that needs care, and the signature says so: it demands the
 * proofs of possession and will not run without them.
 *
 * Aggregating keys over a shared message is insecure without them. An attacker
 * who registers pk_evil = [t]G1 - sum(pk_honest) can produce an aggregate
 * signature for a group whose keys they never held: the honest keys cancel.
 * Proof of possession is what makes that impossible, so rather than document
 * the requirement and hope, this function takes the proofs and checks them. */
int elips_bls_fast_aggregate_verify(const uint8_t *pks,
                                    const uint8_t *pops, size_t n,
                                    const uint8_t *msg, size_t msg_len,
                                    const uint8_t sig[ELIPS_BLS_SIG_BYTES]);

/* --- proof of possession --------------------------------------------------
 * A signature over your own public key, under a separate ciphersuite so it can
 * never be confused with a signature over a message. */
int elips_bls_pop_prove(uint8_t proof[ELIPS_BLS_POP_BYTES],
                        const uint8_t sk[ELIPS_BLS_SK_BYTES]);
int elips_bls_pop_verify(const uint8_t pk[ELIPS_BLS_PK_BYTES],
                         const uint8_t proof[ELIPS_BLS_POP_BYTES]);

/* Human-readable form of the codes above. Never returns NULL. */
const char *elips_bls_strerror(int code);

#ifdef __cplusplus
}
#endif

#endif /* ELIPS_BLS_H */
