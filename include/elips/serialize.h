/*
 * Point serialization in the compressed and uncompressed formats.
 *
 * Layout follows the BLS12-381 convention that the ZCash specification
 * introduced and that IETF draft-irtf-cfrg-pairing-friendly-curves adopted, so
 * a BLS12-381 point written here is byte-identical to one written by any
 * conforming implementation. The same code serves BLS12-461 and BN-462, which
 * have no standard of their own; only the width changes.
 *
 * A field element occupies FP_SER_BYTES bytes, big-endian, and the top three
 * bits of the first byte carry flags:
 *
 *   bit 7   compression   the encoding holds x only, and y is recovered
 *   bit 6   infinity      the point is the identity; every other bit is zero
 *   bit 5   sign          set when y is the lexicographically larger root
 *
 * FP_SER_BYTES is chosen to leave those three bits free. BLS12-381 needs 381
 * bits and gets 48 bytes, which is the standard's own size. BN-462 needs 462
 * bits: 58 bytes would leave only two spare bits, so it uses 59. Deriving the
 * width from the requirement rather than from ceil(FP_BITS/8) is what keeps one
 * implementation correct for all three curves.
 *
 * Fp2 coordinates are written imaginary part first (c1 then c0), again matching
 * the BLS12-381 convention.
 *
 * READERS VALIDATE. Every deserializer rejects non-canonical coordinates,
 * inconsistent flags, points off the curve, and points outside the order-r
 * subgroup, and returns a distinct code for each. This is the trust boundary
 * (plan section 6, Phase 5, item 3): nothing else in the library re-checks, so
 * a caller that ignores the return value and uses the point has been handed an
 * attacker-chosen group element. The subgroup check is the expensive part and
 * the one that matters -- small-subgroup attacks are exactly what it stops.
 */
#ifndef ELIPS_SERIALIZE_H
#define ELIPS_SERIALIZE_H

#include <stdint.h>
#include "elips/ec.h"

/* Bytes per field element, with three bits left over for the flags. */
#define FP_SER_BYTES  ((FP_BITS + 3 + 7) / 8)

#define EP_SER_COMPRESSED_BYTES    (FP_SER_BYTES)
#define EP_SER_UNCOMPRESSED_BYTES  (2 * FP_SER_BYTES)
#define EP2_SER_COMPRESSED_BYTES   (2 * FP_SER_BYTES)
#define EP2_SER_UNCOMPRESSED_BYTES (4 * FP_SER_BYTES)

/* Return codes. Distinct so a caller can tell a malformed encoding from a
 * well-formed encoding of a point that is not allowed. */
#define ELIPS_SER_OK                 0
#define ELIPS_SER_ERR_FLAGS        (-1)   /* flag bits contradict the length */
#define ELIPS_SER_ERR_NONCANONICAL (-2)   /* a coordinate is >= p */
#define ELIPS_SER_ERR_NO_SQRT      (-3)   /* x^3 + b is not a square */
#define ELIPS_SER_ERR_NOT_ON_CURVE (-4)
#define ELIPS_SER_ERR_NOT_IN_GROUP (-5)   /* on the curve, outside order r */

/* --- E(Fp) --- */
void ep_write_compressed(uint8_t out[EP_SER_COMPRESSED_BYTES], const ep_t *p);
void ep_write_uncompressed(uint8_t out[EP_SER_UNCOMPRESSED_BYTES], const ep_t *p);
int  ep_read_compressed(ep_t *p, const uint8_t in[EP_SER_COMPRESSED_BYTES]);
int  ep_read_uncompressed(ep_t *p, const uint8_t in[EP_SER_UNCOMPRESSED_BYTES]);

/* --- E'(Fp2) --- */
void ep2_write_compressed(uint8_t out[EP2_SER_COMPRESSED_BYTES], const ep2_t *p);
void ep2_write_uncompressed(uint8_t out[EP2_SER_UNCOMPRESSED_BYTES], const ep2_t *p);
int  ep2_read_compressed(ep2_t *p, const uint8_t in[EP2_SER_COMPRESSED_BYTES]);
int  ep2_read_uncompressed(ep2_t *p, const uint8_t in[EP2_SER_UNCOMPRESSED_BYTES]);

/* Human-readable form of the codes above, for diagnostics and test output. */
const char *elips_ser_strerror(int code);

#endif /* ELIPS_SERIALIZE_H */
