/*
 * Hash to curve, RFC 9380.
 *
 * Turns an arbitrary byte string into a point of G1 or G2, uniformly and
 * without leaking anything about the string. This is what a BLS signature
 * scheme needs and the reason a pairing library that cannot do it is not yet
 * usable for signatures.
 *
 * Two maps, because the curves need different ones:
 *
 *   BLS12-381             simplified SWU over an 11-isogenous (G1) or
 *                         3-isogenous (G2) curve. This is the suite RFC 9380
 *                         registers, so the output is byte-for-byte what any
 *                         conforming implementation produces. Pinned by the
 *                         specification's own test vectors in
 *                         test/kat/h2c_bls12_381.vec.
 *
 *   BLS12-461, BN-462     Shallue-van de Woestijne. Both curves have A = 0, so
 *                         simplified SWU needs an isogeny that nobody has
 *                         standardised for them; SvdW needs none, and RFC 9380
 *                         uses it for BN254 for exactly this reason. There is
 *                         no registered suite for either curve, so there is
 *                         nothing to interoperate with and these vectors pin
 *                         self-consistency only.
 *
 * DOMAIN SEPARATION. The dst argument is not optional and not decorative. Two
 * protocols that hash with the same tag share a random oracle, and a signature
 * from one can then be replayed into the other. RFC 9380 requires a tag that is
 * unique to the protocol AND to the suite; elips_h2c_suite_g1/g2 give the suite
 * part, and a caller is expected to prefix its own protocol name. A tag longer
 * than 255 bytes is hashed down rather than rejected (RFC 9380 5.3.3), so there
 * is no length to think about.
 *
 * CONSTANT TIME. The message may be a secret -- an oblivious PRF or a
 * password-authenticated key exchange hashes one -- so nothing here branches on
 * the message or on any value derived from it.
 */
#ifndef ELIPS_HASH_TO_CURVE_H
#define ELIPS_HASH_TO_CURVE_H

#include <stddef.h>
#include <stdint.h>
#include "elips/ec.h"
/* Brings in ELIPS_H2C_SUITE_G1 and ELIPS_H2C_SUITE_G2, so a domain separation
 * tag can be built at compile time:
 *
 *     static const char DST[] = "MY-PROTOCOL-V01-CS01-" ELIPS_H2C_SUITE_G1;
 */
#include "elips/h2c_params.h"

/* All four return 0 on success, and -1 only when dst_len is zero.
 *
 * An empty tag is the one input that cannot be served: RFC 9380 requires a
 * non-empty DST, and hashing without one silently removes the domain
 * separation that the argument exists to provide. It is also an easy mistake to
 * make -- an uninitialised length, or sizeof on a pointer -- so it is rejected
 * rather than obeyed. On failure the output is set to the identity, so a caller
 * who ignores the return value gets a useless point rather than an
 * undomain-separated one. Nothing else here can fail. */

/* The random-oracle variants. Two field elements, two maps, one addition, then
 * cofactor clearing. Indifferentiable from a random oracle, which is what
 * signature security proofs assume. */
int elips_hash_to_g1(ep_t *out, const uint8_t *msg, size_t msg_len,
                     const uint8_t *dst, size_t dst_len);
int elips_hash_to_g2(ep2_t *out, const uint8_t *msg, size_t msg_len,
                     const uint8_t *dst, size_t dst_len);

/* The non-uniform variants. One field element instead of two, so roughly half
 * the cost -- and NOT a random oracle: the image is a proper subset of the
 * group and its distribution is not uniform. Safe only where the surrounding
 * proof does not need indifferentiability. When in doubt use the hash_to
 * versions; the saving is one map, not one order of magnitude. */
int elips_encode_to_g1(ep_t *out, const uint8_t *msg, size_t msg_len,
                       const uint8_t *dst, size_t dst_len);
int elips_encode_to_g2(ep2_t *out, const uint8_t *msg, size_t msg_len,
                       const uint8_t *dst, size_t dst_len);

/* The suite identifiers, e.g. "BLS12381G1_XMD:SHA-256_SSWU_RO_". RFC 9380
 * requires these to appear in the domain separation tag. */
const char *elips_h2c_suite_g1(void);
const char *elips_h2c_suite_g2(void);

/* --- the pieces, exposed so they can be tested separately ----------------
 * A failure inside expand_message_xmd and a failure inside the map produce the
 * same symptom -- a wrong point -- so the suite tests both directly. */

/* RFC 9380 5.3.1, with SHA-256. len must be at most 255*32 bytes. Returns 0 on
 * success and -1 if len is out of range or dst_len is zero. */
int elips_expand_message_xmd(uint8_t *out, size_t len,
                             const uint8_t *msg, size_t msg_len,
                             const uint8_t *dst, size_t dst_len);

/* RFC 9380 5.2. Writes count elements; count must be 1 or 2. */
void elips_hash_to_field_fp(fp_t *out, int count,
                            const uint8_t *msg, size_t msg_len,
                            const uint8_t *dst, size_t dst_len);
void elips_hash_to_field_fp2(fp2_t *out, int count,
                             const uint8_t *msg, size_t msg_len,
                             const uint8_t *dst, size_t dst_len);

#endif /* ELIPS_HASH_TO_CURVE_H */
