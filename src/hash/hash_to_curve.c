/*
 * RFC 9380 hash-to-curve. See include/elips/hash_to_curve.h.
 *
 * Structure follows the specification: expand_message_xmd produces uniform
 * bytes, hash_to_field reduces them into the field, map_to_curve sends a field
 * element to the curve, and clear_cofactor lands it in the prime-order
 * subgroup. The two maps live in h2c_tmpl.h, instantiated once per group.
 */
#include "elips/hash_to_curve.h"
#include "elips/h2c_params.h"
#include "elips/pairing.h"
#include "elips/sha256.h"

#include <gmp.h>
#include <string.h>

/* ------------------------------------------------ RFC 9380 5.3.1, expand_xmd */

/* A tag longer than 255 bytes is replaced by H("H2C-OVERSIZE-DST-" || tag)
 * rather than rejected, so no caller has to care about the limit (5.3.3). */
static size_t normalise_dst(uint8_t out[ELIPS_SHA256_DIGEST],
                            const uint8_t **dst, size_t dst_len)
{
    if (dst_len <= 255) return dst_len;

    elips_sha256_t s;
    elips_sha256_init(&s);
    elips_sha256_update(&s, "H2C-OVERSIZE-DST-", 17);
    elips_sha256_update(&s, *dst, dst_len);
    elips_sha256_final(out, &s);
    *dst = out;
    return ELIPS_SHA256_DIGEST;
}

int elips_expand_message_xmd(uint8_t *out, size_t len,
                             const uint8_t *msg, size_t msg_len,
                             const uint8_t *dst, size_t dst_len)
{
    if (dst_len == 0) return -1;              /* RFC 9380 3.1 */

    uint8_t dst_buf[ELIPS_SHA256_DIGEST];
    dst_len = normalise_dst(dst_buf, &dst, dst_len);

    const size_t b_in_bytes = ELIPS_SHA256_DIGEST;   /* 32 */
    const size_t s_in_bytes = ELIPS_SHA256_BLOCK;    /* 64 */
    size_t ell = (len + b_in_bytes - 1) / b_in_bytes;
    if (ell > 255 || len > 65535) return -1;

    uint8_t dst_prime_len = (uint8_t)dst_len;
    uint8_t zpad[ELIPS_SHA256_BLOCK];
    memset(zpad, 0, sizeof zpad);

    uint8_t l_i_b_str[2] = { (uint8_t)(len >> 8), (uint8_t)len };

    /* b_0 = H(Z_pad || msg || l_i_b_str || 0 || DST_prime) */
    uint8_t b0[ELIPS_SHA256_DIGEST];
    {
        elips_sha256_t s;
        elips_sha256_init(&s);
        elips_sha256_update(&s, zpad, s_in_bytes);
        elips_sha256_update(&s, msg, msg_len);
        elips_sha256_update(&s, l_i_b_str, 2);
        elips_sha256_update(&s, "\0", 1);
        elips_sha256_update(&s, dst, dst_len);
        elips_sha256_update(&s, &dst_prime_len, 1);
        elips_sha256_final(b0, &s);
    }

    uint8_t prev[ELIPS_SHA256_DIGEST];
    for (size_t i = 1; i <= ell; i++) {
        elips_sha256_t s;
        elips_sha256_init(&s);
        if (i == 1) {
            elips_sha256_update(&s, b0, b_in_bytes);
        } else {
            uint8_t x[ELIPS_SHA256_DIGEST];
            for (size_t j = 0; j < b_in_bytes; j++) x[j] = (uint8_t)(b0[j] ^ prev[j]);
            elips_sha256_update(&s, x, b_in_bytes);
        }
        uint8_t idx = (uint8_t)i;
        elips_sha256_update(&s, &idx, 1);
        elips_sha256_update(&s, dst, dst_len);
        elips_sha256_update(&s, &dst_prime_len, 1);
        elips_sha256_final(prev, &s);

        size_t off = (i - 1) * b_in_bytes;
        size_t take = len - off < b_in_bytes ? len - off : b_in_bytes;
        memcpy(out + off, prev, take);
    }
    return 0;
}

/* --------------------------------------------- RFC 9380 5.2, hash_to_field */

/* An ELIPS_H2C_L-byte big-endian value reduced modulo p.
 *
 * mpn_sec_div_r rather than a plain division: the message can be a secret, and
 * the whole point of hashing it to a curve is usually to hide it. The operand
 * sizes are compile-time constants, so only they affect the running time. */
#define WIDE_LIMBS ((ELIPS_H2C_L + 7) / 8)

static void reduce_to_fp(fp_t r, const uint8_t *be)
{
    limb_t wide[WIDE_LIMBS];
    mp_limb_t scratch[512];

    memset(wide, 0, sizeof wide);
    for (int i = 0; i < ELIPS_H2C_L; i++) {
        int weight = ELIPS_H2C_L - 1 - i;          /* byte i is the most significant */
        wide[weight / 8] |= (limb_t)be[i] << (8 * (weight % 8));
    }

    mp_size_t itch = mpn_sec_div_r_itch(WIDE_LIMBS, FP_LIMBS);
    if (itch > (mp_size_t)(sizeof scratch / sizeof scratch[0])) {
        fp_set_zero(r);                             /* unreachable at these sizes */
        return;
    }
    mpn_sec_div_r(wide, WIDE_LIMBS, FP_MODULUS, FP_LIMBS, scratch);
    fp_from_limbs(r, wide);
}

void elips_hash_to_field_fp(fp_t *out, int count,
                            const uint8_t *msg, size_t msg_len,
                            const uint8_t *dst, size_t dst_len)
{
    uint8_t buf[2 * ELIPS_H2C_L];
    if (elips_expand_message_xmd(buf, (size_t)count * ELIPS_H2C_L,
                                 msg, msg_len, dst, dst_len) != 0) {
        for (int i = 0; i < count; i++) fp_set_zero(out[i]);
        return;
    }
    for (int i = 0; i < count; i++)
        reduce_to_fp(out[i], buf + (size_t)i * ELIPS_H2C_L);
}

void elips_hash_to_field_fp2(fp2_t *out, int count,
                             const uint8_t *msg, size_t msg_len,
                             const uint8_t *dst, size_t dst_len)
{
    uint8_t buf[4 * ELIPS_H2C_L];
    if (elips_expand_message_xmd(buf, (size_t)count * 2 * ELIPS_H2C_L,
                                 msg, msg_len, dst, dst_len) != 0) {
        for (int i = 0; i < count; i++) fp2_set_zero(out[i]);
        return;
    }
    for (int i = 0; i < count; i++) {
        reduce_to_fp(out[i][0], buf + (size_t)(2 * i)     * ELIPS_H2C_L);
        reduce_to_fp(out[i][1], buf + (size_t)(2 * i + 1) * ELIPS_H2C_L);
    }
}

/* ------------------------------------------------------- the maps to curve */

#define H2C_ELEM(tbl, i) \
    ((const limb_t *)(tbl) + (size_t)(i) * FP_LIMBS)
#define H2C_PT ep
#define H2C_F  fp
#define H2C_FT fp_t
#define H2C_CN H2C_G1_
#include "hash/h2c_tmpl.h"
#undef H2C_PT
#undef H2C_F
#undef H2C_FT
#undef H2C_CN
#undef H2C_ELEM

#define H2C_ELEM(tbl, i) \
    ((const limb_t (*)[FP_LIMBS])((const limb_t *)(tbl) + (size_t)(i) * 2 * FP_LIMBS))
#define H2C_PT ep2
#define H2C_F  fp2
#define H2C_FT fp2_t
#define H2C_CN H2C_G2_
#include "hash/h2c_tmpl.h"
#undef H2C_PT
#undef H2C_F
#undef H2C_FT
#undef H2C_CN
#undef H2C_ELEM

/* -------------------------------------------------------- the entry points */

int elips_hash_to_g1(ep_t *out, const uint8_t *msg, size_t msg_len,
                     const uint8_t *dst, size_t dst_len)
{
    fp_t u[2];
    ep_t q0, q1;
    ep_set_infinity(out);
    if (dst_len == 0) return -1;
    elips_hash_to_field_fp(u, 2, msg, msg_len, dst, dst_len);
    ep_h2c_map_to_curve(&q0, u[0]);
    ep_h2c_map_to_curve(&q1, u[1]);
    ep_add(&q0, &q0, &q1);
    ep_mul(out, &q0, H2C_HEFF_G1, H2C_HEFF_G1_BITS);
    return 0;
}

int elips_encode_to_g1(ep_t *out, const uint8_t *msg, size_t msg_len,
                       const uint8_t *dst, size_t dst_len)
{
    fp_t u[1];
    ep_t q;
    ep_set_infinity(out);
    if (dst_len == 0) return -1;
    elips_hash_to_field_fp(u, 1, msg, msg_len, dst, dst_len);
    ep_h2c_map_to_curve(&q, u[0]);
    ep_mul(out, &q, H2C_HEFF_G1, H2C_HEFF_G1_BITS);
    return 0;
}

/* clear_cofactor on G2.
 *
 * The plain version multiplies by h_eff, which is 636 bits on BLS12-381 and 769
 * on BLS12-461 -- so the cofactor clearing, not the map, dominates
 * elips_hash_to_g2. Budroni and Pintore (ePrint 2017/419) replace it with
 *
 *     [h_eff]Q = [x^2 - x - 1]Q + [x - 1]psi(Q) + psi^2([2]Q)
 *
 * where x is the mother parameter: 64 bits on BLS12-381 against 636. Two short
 * ladders and three applications of psi.
 *
 * The identity is verified numerically in tools/reference/h2c_ref.py against
 * the h_eff the generator derived, on random points of the TWIST rather than of
 * G2 -- on G2 much weaker relations hold and would hide a wrong chain. Beyond
 * that, a wrong chain here cannot pass silently: it computes a different
 * multiple, and the RFC 9380 vectors stop matching.
 *
 * BN gets its own chain, and it is simpler than the BLS12 one. There the G2
 * cofactor is exactly h2 = p + t - 1 = p + 6x^2, so in base p it is
 * 6x^2 + 1*p: the multiplier by p is one and the remainder is only 231 bits.
 * On the twist psi satisfies psi^2 - [t]psi + [p] = 0, so [p] = [t]psi - psi^2
 * and p disappears entirely:
 *
 *     [h2]Q = [6x^2]Q + [t]psi(Q) - psi^2(Q)
 *           = [6x^2](Q + psi(Q)) + psi(Q) - psi^2(Q)
 *
 * folding the two scalar multiplications into one with t = 6x^2 + 1. A SINGLE
 * 231-bit ladder replaces the 462-bit one, and it computes exactly [h_eff],
 * so no vector moves. ELIPS_6XSQ already exists for the G2 subgroup test, so
 * the chain needs no new constants.
 *
 * Published fast chains for BN G2 generally compute a different MULTIPLE of
 * the cofactor, which would change what hash_to_g2 returns. This one does
 * not, which is why it could be adopted without regenerating anything.
 * tools/reference/h2c_ref.py derives it and checks it against [h_eff] on
 * random points of the twist. */
static void clear_cofactor_g2(ep2_t *r, const ep2_t *q)
{
#ifdef ELIPS_H2C_G2_FAST_CLEAR
    ep2_t base0, base1, t2;

    /* [A]Q + [B]psi(Q) in ONE interleaved ladder rather than two separate
     * ones. The two ladders were never able to share a table -- their base
     * points are Q and psi(Q), which are different points -- so the saving is
     * in the doublings, not the table: one set of A_BITS doublings instead of
     * A_BITS + B_BITS.
     *
     * Measured on BLS12-381, where clear_cofactor is 830 us of hash_to_g2's
     * 1841: the two ladders were 541 + 289 us and the interleaved one is
     * about 625.
     *
     * A negative multiplier is applied to the POINT, since ep2_mul2 takes
     * unsigned scalars and negating a point is free. */
    ep2_copy(&base0, q);
#if H2C_G2_CLEAR_A_NEG
    ep2_neg(&base0, &base0);
#endif
    ep2_psi(&base1, q);
#if H2C_G2_CLEAR_B_NEG
    ep2_neg(&base1, &base1);
#endif

    ep2_mul2(r, &base0, H2C_G2_CLEAR_A, H2C_G2_CLEAR_A_BITS,
                &base1, H2C_G2_CLEAR_B, H2C_G2_CLEAR_B_BITS);

    ep2_dbl(&t2, q);
    ep2_psi(&t2, &t2);
    ep2_psi(&t2, &t2);

    ep2_add(r, r, &t2);
#elif defined(ELIPS_FAMILY_BN)
    ep2_t pq, ppq, sum;

    ep2_psi(&pq, q);          /* psi(Q)    */
    ep2_psi(&ppq, &pq);       /* psi^2(Q)  */
    ep2_neg(&ppq, &ppq);

    ep2_add(&sum, q, &pq);    /* Q + psi(Q) */
    ep2_mul(r, &sum, ELIPS_6XSQ, ELIPS_6XSQ_BITS);
    ep2_add(r, r, &pq);
    ep2_add(r, r, &ppq);
#else
    ep2_mul(r, q, H2C_HEFF_G2, H2C_HEFF_G2_BITS);
#endif
}

int elips_hash_to_g2(ep2_t *out, const uint8_t *msg, size_t msg_len,
                     const uint8_t *dst, size_t dst_len)
{
    fp2_t u[2];
    ep2_t q0, q1;
    ep2_set_infinity(out);
    if (dst_len == 0) return -1;
    elips_hash_to_field_fp2(u, 2, msg, msg_len, dst, dst_len);
    ep2_h2c_map_to_curve(&q0, u[0]);
    ep2_h2c_map_to_curve(&q1, u[1]);
    ep2_add(&q0, &q0, &q1);
    clear_cofactor_g2(out, &q0);
    return 0;
}

int elips_encode_to_g2(ep2_t *out, const uint8_t *msg, size_t msg_len,
                       const uint8_t *dst, size_t dst_len)
{
    fp2_t u[1];
    ep2_t q;
    ep2_set_infinity(out);
    if (dst_len == 0) return -1;
    elips_hash_to_field_fp2(u, 1, msg, msg_len, dst, dst_len);
    ep2_h2c_map_to_curve(&q, u[0]);
    clear_cofactor_g2(out, &q);
    return 0;
}

const char *elips_h2c_suite_g1(void) { return ELIPS_H2C_SUITE_G1; }
const char *elips_h2c_suite_g2(void) { return ELIPS_H2C_SUITE_G2; }
