/*
 * Point serialization. See include/elips/serialize.h for the format.
 *
 * Nothing here is on a hot path, so it is written for auditability: every
 * validation step is a named function with one job, and the readers run them in
 * a fixed order so an error code always identifies the first thing that was
 * wrong rather than whichever check happened to fire.
 */
#include "elips/serialize.h"
#include "elips/pairing.h"

#include <string.h>

#define FLAG_COMPRESSED 0x80u
#define FLAG_INFINITY   0x40u
#define FLAG_SIGN       0x20u
#define FLAG_MASK       0xe0u

/* --- field element to and from FP_SER_BYTES big-endian bytes --------------- */

static void fp_write_be(uint8_t *out, const fp_t a)
{
    limb_t plain[FP_LIMBS];
    fp_to_limbs(plain, a);                       /* canonical, already < p */

    memset(out, 0, FP_SER_BYTES);
    for (int i = 0; i < (int)FP_LIMBS; i++) {
        for (int b = 0; b < 8; b++) {
            int pos = FP_SER_BYTES - 1 - (i * 8 + b);
            if (pos >= 0) out[pos] = (uint8_t)(plain[i] >> (8 * b));
        }
    }
}

/* Returns 0 on success, ELIPS_SER_ERR_NONCANONICAL if the encoded integer is
 * not already reduced. Two encodings of one element would otherwise both be
 * accepted, which breaks any protocol that hashes or compares serialized
 * points. */
static int fp_read_be(fp_t r, const uint8_t *in)
{
    limb_t plain[FP_LIMBS];
    memset(plain, 0, sizeof plain);

    for (int i = 0; i < FP_SER_BYTES; i++) {
        int idx = FP_SER_BYTES - 1 - i;          /* byte weight, from the end */
        if (i / 8 < (int)FP_LIMBS)
            plain[i / 8] |= (limb_t)in[idx] << (8 * (i % 8));
        else if (in[idx] != 0)
            return ELIPS_SER_ERR_NONCANONICAL;   /* above the limb array */
    }

    /* plain < p, tested by the borrow out of plain - p. */
    limb_t borrow = 0;
    for (int i = 0; i < (int)FP_LIMBS; i++) {
        limb_t x = plain[i], y = FP_MODULUS[i];
        limb_t d = x - y, b1 = (x < y), d2 = d - borrow, b2 = (d < borrow);
        (void)d2;
        borrow = b1 | b2;
    }
    if (!borrow) return ELIPS_SER_ERR_NONCANONICAL;

    fp_from_limbs(r, plain);
    return 0;
}

/* --- flag handling --------------------------------------------------------- */

/* Reads and clears the flag byte. Rejects any combination the format does not
 * define, so a decoder cannot be steered by bits the writer never sets. */
static int take_flags(uint8_t *scratch, const uint8_t *in, size_t len,
                      int want_compressed, int *is_infinity, int *sign)
{
    memcpy(scratch, in, len);
    uint8_t f = scratch[0];

    int compressed = (f & FLAG_COMPRESSED) != 0;
    *is_infinity   = (f & FLAG_INFINITY)   != 0;
    *sign          = (f & FLAG_SIGN)       != 0;

    if (compressed != want_compressed)          return ELIPS_SER_ERR_FLAGS;
    if (!compressed && *sign)                   return ELIPS_SER_ERR_FLAGS;

    scratch[0] = (uint8_t)(f & (uint8_t)~FLAG_MASK);

    if (*is_infinity) {
        /* The identity has one encoding: flags, then nothing. */
        for (size_t i = 0; i < len; i++)
            if (scratch[i] != 0) return ELIPS_SER_ERR_FLAGS;
        if (*sign) return ELIPS_SER_ERR_FLAGS;
    }
    return 0;
}

static void put_flags(uint8_t *out, int compressed, int infinity, int sign)
{
    uint8_t f = 0;
    if (compressed) f |= FLAG_COMPRESSED;
    if (infinity)   f |= FLAG_INFINITY;
    if (sign)       f |= FLAG_SIGN;
    out[0] = (uint8_t)(out[0] | f);
}

/* --- E(Fp) ----------------------------------------------------------------- */

void ep_write_compressed(uint8_t out[EP_SER_COMPRESSED_BYTES], const ep_t *p)
{
    memset(out, 0, EP_SER_COMPRESSED_BYTES);
    fp_t x, y;
    if (!ep_to_affine(x, y, p)) { put_flags(out, 1, 1, 0); return; }
    fp_write_be(out, x);
    put_flags(out, 1, 0, fp_is_lex_largest(y));
}

void ep_write_uncompressed(uint8_t out[EP_SER_UNCOMPRESSED_BYTES], const ep_t *p)
{
    memset(out, 0, EP_SER_UNCOMPRESSED_BYTES);
    fp_t x, y;
    if (!ep_to_affine(x, y, p)) { put_flags(out, 0, 1, 0); return; }
    fp_write_be(out, x);
    fp_write_be(out + FP_SER_BYTES, y);
    put_flags(out, 0, 0, 0);
}

/* y^2 = x^3 + b, with the root chosen so that its sign bit matches want_sign. */
static int ep_recover_y(fp_t y, const fp_t x, int want_sign)
{
    fp_t rhs, b;
    ep_curve_b(b);
    fp_sqr(rhs, x);
    fp_mul(rhs, rhs, x);
    fp_add(rhs, rhs, b);

    if (!fp_sqrt(y, rhs)) return ELIPS_SER_ERR_NO_SQRT;

    fp_t neg;
    fp_neg(neg, y);
    limb_t flip = (limb_t)0 - (limb_t)(fp_is_lex_largest(y) != want_sign);
    fp_cselect(y, neg, y, flip);
    return 0;
}

/* The order the checks run in is the order of the return codes: shape, then
 * canonicity, then the curve equation, then the subgroup. */
static int ep_finish(ep_t *p, const fp_t x, const fp_t y)
{
    ep_from_affine(p, x, y);
    if (!ep_on_curve(p))    { ep_set_infinity(p); return ELIPS_SER_ERR_NOT_ON_CURVE; }
    if (!ep_in_subgroup(p)) { ep_set_infinity(p); return ELIPS_SER_ERR_NOT_IN_GROUP; }
    return ELIPS_SER_OK;
}

int ep_read_compressed(ep_t *p, const uint8_t in[EP_SER_COMPRESSED_BYTES])
{
    uint8_t buf[EP_SER_COMPRESSED_BYTES];
    int inf, sign, rc;

    ep_set_infinity(p);
    if ((rc = take_flags(buf, in, sizeof buf, 1, &inf, &sign)) != 0) return rc;
    if (inf) return ELIPS_SER_OK;

    fp_t x, y;
    if ((rc = fp_read_be(x, buf)) != 0) return rc;
    if ((rc = ep_recover_y(y, x, sign)) != 0) return rc;
    return ep_finish(p, x, y);
}

int ep_read_uncompressed(ep_t *p, const uint8_t in[EP_SER_UNCOMPRESSED_BYTES])
{
    uint8_t buf[EP_SER_UNCOMPRESSED_BYTES];
    int inf, sign, rc;

    ep_set_infinity(p);
    if ((rc = take_flags(buf, in, sizeof buf, 0, &inf, &sign)) != 0) return rc;
    if (inf) return ELIPS_SER_OK;

    fp_t x, y;
    if ((rc = fp_read_be(x, buf)) != 0)                return rc;
    if ((rc = fp_read_be(y, buf + FP_SER_BYTES)) != 0) return rc;
    return ep_finish(p, x, y);
}

/* --- E'(Fp2) --------------------------------------------------------------- */
/* Imaginary part first, matching the BLS12-381 convention. */

static void fp2_write_be(uint8_t *out, const fp2_t a)
{
    fp_write_be(out, a[1]);
    fp_write_be(out + FP_SER_BYTES, a[0]);
}

static int fp2_read_be(fp2_t r, const uint8_t *in)
{
    int rc;
    if ((rc = fp_read_be(r[1], in)) != 0)                return rc;
    if ((rc = fp_read_be(r[0], in + FP_SER_BYTES)) != 0) return rc;
    return 0;
}

void ep2_write_compressed(uint8_t out[EP2_SER_COMPRESSED_BYTES], const ep2_t *p)
{
    memset(out, 0, EP2_SER_COMPRESSED_BYTES);
    fp2_t x, y;
    if (!ep2_to_affine(x, y, p)) { put_flags(out, 1, 1, 0); return; }
    fp2_write_be(out, x);
    put_flags(out, 1, 0, fp2_is_lex_largest(y));
}

void ep2_write_uncompressed(uint8_t out[EP2_SER_UNCOMPRESSED_BYTES], const ep2_t *p)
{
    memset(out, 0, EP2_SER_UNCOMPRESSED_BYTES);
    fp2_t x, y;
    if (!ep2_to_affine(x, y, p)) { put_flags(out, 0, 1, 0); return; }
    fp2_write_be(out, x);
    fp2_write_be(out + 2 * FP_SER_BYTES, y);
    put_flags(out, 0, 0, 0);
}

static int ep2_recover_y(fp2_t y, const fp2_t x, int want_sign)
{
    fp2_t rhs, b;
    ep2_curve_b(b);
    fp2_sqr(rhs, x);
    fp2_mul(rhs, rhs, x);
    fp2_add(rhs, rhs, b);

    if (!fp2_sqrt(y, rhs)) return ELIPS_SER_ERR_NO_SQRT;

    fp2_t neg;
    fp2_neg(neg, y);
    limb_t flip = (limb_t)0 - (limb_t)(fp2_is_lex_largest(y) != want_sign);
    fp2_cselect(y, neg, y, flip);
    return 0;
}

static int ep2_finish(ep2_t *p, const fp2_t x, const fp2_t y)
{
    ep2_from_affine(p, x, y);
    if (!ep2_on_curve(p))    { ep2_set_infinity(p); return ELIPS_SER_ERR_NOT_ON_CURVE; }
    if (!ep2_in_subgroup(p)) { ep2_set_infinity(p); return ELIPS_SER_ERR_NOT_IN_GROUP; }
    return ELIPS_SER_OK;
}

int ep2_read_compressed(ep2_t *p, const uint8_t in[EP2_SER_COMPRESSED_BYTES])
{
    uint8_t buf[EP2_SER_COMPRESSED_BYTES];
    int inf, sign, rc;

    ep2_set_infinity(p);
    if ((rc = take_flags(buf, in, sizeof buf, 1, &inf, &sign)) != 0) return rc;
    if (inf) return ELIPS_SER_OK;

    fp2_t x, y;
    if ((rc = fp2_read_be(x, buf)) != 0) return rc;
    if ((rc = ep2_recover_y(y, x, sign)) != 0) return rc;
    return ep2_finish(p, x, y);
}

int ep2_read_uncompressed(ep2_t *p, const uint8_t in[EP2_SER_UNCOMPRESSED_BYTES])
{
    uint8_t buf[EP2_SER_UNCOMPRESSED_BYTES];
    int inf, sign, rc;

    ep2_set_infinity(p);
    if ((rc = take_flags(buf, in, sizeof buf, 0, &inf, &sign)) != 0) return rc;
    if (inf) return ELIPS_SER_OK;

    fp2_t x, y;
    if ((rc = fp2_read_be(x, buf)) != 0)                    return rc;
    if ((rc = fp2_read_be(y, buf + 2 * FP_SER_BYTES)) != 0) return rc;
    return ep2_finish(p, x, y);
}

const char *elips_ser_strerror(int code)
{
    switch (code) {
    case ELIPS_SER_OK:                 return "ok";
    case ELIPS_SER_ERR_FLAGS:          return "malformed flag bits";
    case ELIPS_SER_ERR_NONCANONICAL:   return "coordinate not reduced mod p";
    case ELIPS_SER_ERR_NO_SQRT:        return "x has no corresponding y";
    case ELIPS_SER_ERR_NOT_ON_CURVE:   return "point is not on the curve";
    case ELIPS_SER_ERR_NOT_IN_GROUP:   return "point is outside the order-r subgroup";
    default:                           return "unknown error";
    }
}
