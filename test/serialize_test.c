/*
 * Point serialization: round trips, the format's fixed points, and above all
 * the rejections.
 *
 * A serializer is easy to test optimistically -- write, read, compare -- and
 * that half of the test would pass against a decoder that validates nothing.
 * The half that matters is the negative one: every malformed encoding an
 * attacker can send has to come back as an error, not as a point. So each
 * failure mode gets its own case, and each asserts the specific error code, not
 * merely "some error", so that a decoder cannot pass by rejecting everything.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <string.h>
#include "elips/serialize.h"
#include "elips/pairing.h"
#include "elips/random.h"

static int fails;
static void ok(int c, const char *w)
{ printf("  [%s] %s\n", c ? "PASS" : "FAIL", w); if (!c) fails++; }

static void okc(int got, int want, const char *w)
{
    int c = (got == want);
    printf("  [%s] %s (got %s)\n", c ? "PASS" : "FAIL", w,
           elips_ser_strerror(got));
    if (!c) fails++;
}

/* --- G1 -------------------------------------------------------------------- */

static void test_ep(void)
{
    uint8_t c[EP_SER_COMPRESSED_BYTES], u[EP_SER_UNCOMPRESSED_BYTES];
    ep_t P, R;
    ep_generator(&P);

    ep_write_compressed(c, &P);
    okc(ep_read_compressed(&R, c), ELIPS_SER_OK, "G1 compressed round trip reads");
    {
        fp_t x1, y1, x2, y2;
        ep_to_affine(x1, y1, &P); ep_to_affine(x2, y2, &R);
        ok(fp_eq(x1, x2) && fp_eq(y1, y2), "G1 compressed round trip is exact");
    }
    ok((c[0] & 0x80u) != 0, "G1 compressed encoding sets the compression bit");

    ep_write_uncompressed(u, &P);
    okc(ep_read_uncompressed(&R, u), ELIPS_SER_OK, "G1 uncompressed round trip reads");
    {
        fp_t x1, y1, x2, y2;
        ep_to_affine(x1, y1, &P); ep_to_affine(x2, y2, &R);
        ok(fp_eq(x1, x2) && fp_eq(y1, y2), "G1 uncompressed round trip is exact");
    }
    ok((u[0] & 0xe0u) == 0, "G1 uncompressed encoding sets no flag bits");

    /* Both roots must be distinguishable, or the sign bit is doing nothing. */
    {
        ep_t N; uint8_t cn[EP_SER_COMPRESSED_BYTES];
        ep_neg(&N, &P);
        ep_write_compressed(cn, &N);
        ok(memcmp(c, cn, sizeof c) != 0, "P and -P have different compressed forms");
        ok(memcmp(c + 1, cn + 1, sizeof c - 1) == 0,
           "P and -P differ only in the sign bit");
        okc(ep_read_compressed(&R, cn), ELIPS_SER_OK, "-P decompresses");
        {
            fp_t x1, y1, x2, y2;
            ep_to_affine(x1, y1, &N); ep_to_affine(x2, y2, &R);
            ok(fp_eq(x1, x2) && fp_eq(y1, y2), "-P round trips to -P, not P");
        }
    }

    /* The identity. */
    {
        ep_t O; ep_set_infinity(&O);
        ep_write_compressed(c, &O);
        ok(c[0] == 0xc0u, "compressed identity is 0xc0 then zeros");
        for (size_t i = 1; i < sizeof c; i++) if (c[i]) { ok(0, "identity tail is zero"); break; }
        okc(ep_read_compressed(&R, c), ELIPS_SER_OK, "compressed identity reads");
        ok(ep_is_infinity(&R), "compressed identity round trips");

        ep_write_uncompressed(u, &O);
        ok(u[0] == 0x40u, "uncompressed identity is 0x40 then zeros");
        okc(ep_read_uncompressed(&R, u), ELIPS_SER_OK, "uncompressed identity reads");
        ok(ep_is_infinity(&R), "uncompressed identity round trips");
    }

    /* --- rejections --- */
    ep_write_compressed(c, &P);

    {   /* compression bit clear in a compressed-length buffer */
        uint8_t bad[EP_SER_COMPRESSED_BYTES];
        memcpy(bad, c, sizeof bad); bad[0] &= (uint8_t)~0x80u;
        okc(ep_read_compressed(&R, bad), ELIPS_SER_ERR_FLAGS,
            "G1 rejects a cleared compression bit");
    }
    {   /* infinity bit set with a non-zero body */
        uint8_t bad[EP_SER_COMPRESSED_BYTES];
        memcpy(bad, c, sizeof bad); bad[0] |= 0x40u;
        okc(ep_read_compressed(&R, bad), ELIPS_SER_ERR_FLAGS,
            "G1 rejects infinity flag over a non-zero body");
    }
    {   /* x = p, the canonical off-by-one an encoder must never emit */
        uint8_t bad[EP_SER_COMPRESSED_BYTES];
        memset(bad, 0, sizeof bad);
        limb_t plain[FP_LIMBS];
        memcpy(plain, FP_MODULUS, sizeof plain);
        for (int i = 0; i < (int)FP_LIMBS; i++)
            for (int b = 0; b < 8; b++) {
                int pos = (int)sizeof bad - 1 - (i * 8 + b);
                if (pos >= 0) bad[pos] = (uint8_t)(plain[i] >> (8 * b));
            }
        bad[0] |= 0x80u;
        okc(ep_read_compressed(&R, bad), ELIPS_SER_ERR_NONCANONICAL,
            "G1 rejects x = p");
    }
    {   /* an x with no matching y. Search deterministically from 1 upward:
         * roughly half of all x are non-residues, so this terminates fast. */
        uint8_t bad[EP_SER_COMPRESSED_BYTES];
        int found = 0;
        for (unsigned v = 1; v < 64 && !found; v++) {
            memset(bad, 0, sizeof bad);
            bad[sizeof bad - 1] = (uint8_t)v;
            bad[0] |= 0x80u;
            if (ep_read_compressed(&R, bad) == ELIPS_SER_ERR_NO_SQRT) found = 1;
        }
        ok(found, "G1 rejects an x with no square root");
    }
    {   /* on the curve is not enough: an uncompressed point off the curve */
        uint8_t bad[EP_SER_UNCOMPRESSED_BYTES];
        ep_write_uncompressed(bad, &P);
        bad[sizeof bad - 1] ^= 1u;         /* perturb y */
        okc(ep_read_uncompressed(&R, bad), ELIPS_SER_ERR_NOT_ON_CURVE,
            "G1 rejects a point off the curve");
    }
    {   /* sign bit is meaningless in the uncompressed form */
        uint8_t bad[EP_SER_UNCOMPRESSED_BYTES];
        ep_write_uncompressed(bad, &P);
        bad[0] |= 0x20u;
        okc(ep_read_uncompressed(&R, bad), ELIPS_SER_ERR_FLAGS,
            "G1 rejects a sign bit on an uncompressed point");
    }
}

/* --- G2 -------------------------------------------------------------------- */

static void test_ep2(void)
{
    uint8_t c[EP2_SER_COMPRESSED_BYTES], u[EP2_SER_UNCOMPRESSED_BYTES];
    ep2_t Q, R;
    ep2_generator(&Q);

    ep2_write_compressed(c, &Q);
    okc(ep2_read_compressed(&R, c), ELIPS_SER_OK, "G2 compressed round trip reads");
    {
        fp2_t x1, y1, x2, y2;
        ep2_to_affine(x1, y1, &Q); ep2_to_affine(x2, y2, &R);
        ok(fp2_eq(x1, x2) && fp2_eq(y1, y2), "G2 compressed round trip is exact");
    }

    ep2_write_uncompressed(u, &Q);
    okc(ep2_read_uncompressed(&R, u), ELIPS_SER_OK, "G2 uncompressed round trip reads");
    {
        fp2_t x1, y1, x2, y2;
        ep2_to_affine(x1, y1, &Q); ep2_to_affine(x2, y2, &R);
        ok(fp2_eq(x1, x2) && fp2_eq(y1, y2), "G2 uncompressed round trip is exact");
    }

    {
        ep2_t N; uint8_t cn[EP2_SER_COMPRESSED_BYTES];
        ep2_neg(&N, &Q);
        ep2_write_compressed(cn, &N);
        ok(memcmp(c + 1, cn + 1, sizeof c - 1) == 0,
           "G2: Q and -Q differ only in the sign bit");
        okc(ep2_read_compressed(&R, cn), ELIPS_SER_OK, "-Q decompresses");
        {
            fp2_t x1, y1, x2, y2;
            ep2_to_affine(x1, y1, &N); ep2_to_affine(x2, y2, &R);
            ok(fp2_eq(x1, x2) && fp2_eq(y1, y2), "-Q round trips to -Q");
        }
    }

    {
        ep2_t O; ep2_set_infinity(&O);
        ep2_write_compressed(c, &O);
        ok(c[0] == 0xc0u, "G2 compressed identity is 0xc0 then zeros");
        okc(ep2_read_compressed(&R, c), ELIPS_SER_OK, "G2 compressed identity reads");
        ok(ep2_is_infinity(&R), "G2 compressed identity round trips");
    }

    ep2_write_compressed(c, &Q);
    {
        uint8_t bad[EP2_SER_COMPRESSED_BYTES];
        memcpy(bad, c, sizeof bad); bad[0] &= (uint8_t)~0x80u;
        okc(ep2_read_compressed(&R, bad), ELIPS_SER_ERR_FLAGS,
            "G2 rejects a cleared compression bit");
    }
    {
        uint8_t bad[EP2_SER_UNCOMPRESSED_BYTES];
        ep2_write_uncompressed(bad, &Q);
        bad[sizeof bad - 1] ^= 1u;
        okc(ep2_read_uncompressed(&R, bad), ELIPS_SER_ERR_NOT_ON_CURVE,
            "G2 rejects a point off the twist");
    }
}

/* --- on-curve but outside the subgroup ------------------------------------- */
/*
 * The check that actually stops small-subgroup attacks, and the one an
 * implementation is most likely to skip because everything still appears to
 * work without it.
 *
 * Finding such a point without a cofactor-clearing map: search x upward for one
 * where x^3 + b is a square, which gives a point of E(Fp). Its order divides
 * #E = r*h, and only a 1/h fraction of them land in the order-r subgroup, so a
 * short search finds one outside it.
 */
static void test_off_subgroup(void)
{
    uint8_t enc[EP_SER_COMPRESSED_BYTES];
    ep_t R;
    int found_any = 0, found_off = 0;

    for (unsigned v = 1; v < 256 && !found_off; v++) {
        memset(enc, 0, sizeof enc);
        enc[sizeof enc - 1] = (uint8_t)v;
        enc[0] |= 0x80u;
        int rc = ep_read_compressed(&R, enc);
        if (rc == ELIPS_SER_OK)                    found_any = 1;
        if (rc == ELIPS_SER_ERR_NOT_IN_GROUP)      found_off = 1;
    }

    if (found_off) {
        ok(1, "G1 rejects an on-curve point outside the order-r subgroup");
    } else if (found_any) {
        /* A prime-order curve, or a cofactor small enough that the search
         * missed. Say so rather than reporting a pass that was never tested. */
        printf("  [SKIP] no off-subgroup point in the search range"
               " (cofactor may be 1)\n");
    } else {
        ok(0, "the search found no valid x at all -- the decoder is broken");
    }
}

/* --- randomness ------------------------------------------------------------ */
/*
 * Not a statistical test: the source is the operating system's and testing its
 * distribution here would test the kernel, not this code. What is tested is the
 * part this code owns -- that it calls the source at all, reduces into range,
 * and does not return the same thing twice.
 */
static void test_random(void)
{
    fp_t a, b;
    ok(fp_rand(a) == 0 && fp_rand(b) == 0, "fp_rand succeeds");
    ok(!fp_eq(a, b), "two field samples differ");

    limb_t k1[ELIPS_ORDER_LIMBS], k2[ELIPS_ORDER_LIMBS];
    ok(elips_random_scalar(k1) == 0 && elips_random_scalar(k2) == 0,
       "elips_random_scalar succeeds");
    ok(memcmp(k1, k2, sizeof k1) != 0, "two scalar samples differ");
    ok(elips_scalar_is_reduced(k1) && elips_scalar_is_reduced(k2),
       "sampled scalars are below the group order");

    {   /* r itself and r+1 must be rejected; r-1 accepted. */
        limb_t k[ELIPS_ORDER_LIMBS];
        memcpy(k, ELIPS_ORDER, sizeof k);
        ok(!elips_scalar_is_reduced(k), "the group order itself is not reduced");
        k[0] += 1;
        ok(!elips_scalar_is_reduced(k), "r+1 is not reduced");
        memcpy(k, ELIPS_ORDER, sizeof k);
        k[0] -= 1;
        ok(elips_scalar_is_reduced(k), "r-1 is reduced");
        memset(k, 0, sizeof k);
        ok(elips_scalar_is_reduced(k), "zero is reduced");
    }

    {   /* A random scalar has to be usable end to end. */
        limb_t k[ELIPS_ORDER_LIMBS];
        ok(elips_random_scalar(k) == 0, "sampled a scalar for the ladder");
        ep_t P, kP;
        ep_generator(&P);
        ep_mul(&kP, &P, k, ELIPS_ORDER_BITS);
        ok(ep_in_subgroup(&kP), "[k]G1 is in the subgroup for random k");

        uint8_t enc[EP_SER_COMPRESSED_BYTES];
        ep_t back;
        ep_write_compressed(enc, &kP);
        okc(ep_read_compressed(&back, enc), ELIPS_SER_OK,
            "a random group element round trips through the encoding");
    }
}

/* --- square roots ---------------------------------------------------------- */
/*
 * The decompressors rest on these, so test them where they are, not only
 * through their callers. A wrong square root that happens to satisfy the curve
 * equation for the generator would otherwise slip through.
 */
static void test_sqrt(void)
{
    int residues = 0, nonresidues = 0, wrong = 0;
    for (int i = 0; i < 64; i++) {
        fp_t a, s, chk;
        if (fp_rand(a) != 0) { ok(0, "fp_rand failed inside the sqrt test"); return; }
        if (fp_sqrt(s, a)) {
            residues++;
            fp_sqr(chk, s);
            if (!fp_eq(chk, a)) wrong++;
        } else {
            nonresidues++;
            /* A non-residue must not come back with a root anyway. */
            fp_t zero; fp_set_zero(zero);
            if (!fp_eq(s, zero)) wrong++;
        }
        /* a^2 is always a residue, whatever a was. */
        fp_t sq, root;
        fp_sqr(sq, a);
        if (!fp_sqrt(root, sq)) wrong++;
        else { fp_sqr(chk, root); if (!fp_eq(chk, sq)) wrong++; }
    }
    ok(wrong == 0, "fp_sqrt: 64 random inputs, every root verified");
    ok(residues > 0 && nonresidues > 0,
       "fp_sqrt saw both residues and non-residues");

    int wrong2 = 0, res2 = 0, non2 = 0;
    for (int i = 0; i < 32; i++) {
        fp2_t a, s, chk, sq, root;
        if (fp_rand(a[0]) != 0 || fp_rand(a[1]) != 0) {
            ok(0, "fp_rand failed inside the fp2 sqrt test"); return;
        }
        if (fp2_sqrt(s, a)) { res2++; fp2_sqr(chk, s); if (!fp2_eq(chk, a)) wrong2++; }
        else                { non2++; }
        fp2_sqr(sq, a);
        if (!fp2_sqrt(root, sq)) wrong2++;
        else { fp2_sqr(chk, root); if (!fp2_eq(chk, sq)) wrong2++; }
    }
    ok(wrong2 == 0, "fp2_sqrt: 32 random inputs, every root verified");
    ok(res2 > 0 && non2 > 0, "fp2_sqrt saw both residues and non-residues");
}

/* --- BLS12-381 interoperability --------------------------------------------
 *
 * The one curve here with a specification, so the one place the format can be
 * checked against something other than itself. These are the published
 * generator encodings; if this test passes, a BLS12-381 point written by this
 * library is the same bytes any conforming implementation would write.
 *
 * It also pins the generators themselves. A library can have a byte-exact
 * format and still fail to interoperate by handing out a different generator of
 * the same subgroup -- which this one did until the parameters were regenerated
 * from the specification.
 */
#ifdef ELIPS_CURVE_BLS12_381
static const char *KAT_G1_COMPRESSED =
    "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905"
    "a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb";
static const char *KAT_G1_UNCOMPRESSED =
    "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905"
    "a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"
    "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af6"
    "00db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
static const char *KAT_G2_COMPRESSED =
    "93e02b6052719f607dacd3a088274f65596bd0d09920b61a"
    "b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e"
    "024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02"
    "b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8";

static int hex_eq(const uint8_t *b, size_t n, const char *h)
{
    if (strlen(h) != 2 * n) return 0;
    for (size_t i = 0; i < n; i++) {
        char buf[3];
        snprintf(buf, sizeof buf, "%02x", b[i]);
        if (buf[0] != h[2 * i] || buf[1] != h[2 * i + 1]) return 0;
    }
    return 1;
}

static void test_kat(void)
{
    ep_t P; ep2_t Q;
    ep_generator(&P); ep2_generator(&Q);

    uint8_t c1[EP_SER_COMPRESSED_BYTES], u1[EP_SER_UNCOMPRESSED_BYTES];
    uint8_t c2[EP2_SER_COMPRESSED_BYTES];

    ep_write_compressed(c1, &P);
    ok(hex_eq(c1, sizeof c1, KAT_G1_COMPRESSED),
       "G1 generator matches the published compressed encoding");
    ep_write_uncompressed(u1, &P);
    ok(hex_eq(u1, sizeof u1, KAT_G1_UNCOMPRESSED),
       "G1 generator matches the published uncompressed encoding");
    ep2_write_compressed(c2, &Q);
    ok(hex_eq(c2, sizeof c2, KAT_G2_COMPRESSED),
       "G2 generator matches the published compressed encoding");
}
#endif

int main(void)
{
    printf("serialization and randomness [%s]\n", ELIPS_CURVE_NAME);
    printf("  sizes: G1 %d/%d, G2 %d/%d bytes (compressed/uncompressed)\n",
           EP_SER_COMPRESSED_BYTES, EP_SER_UNCOMPRESSED_BYTES,
           EP2_SER_COMPRESSED_BYTES, EP2_SER_UNCOMPRESSED_BYTES);

    test_sqrt();
#ifdef ELIPS_CURVE_BLS12_381
    test_kat();
#endif
    test_ep();
    test_ep2();
    test_off_subgroup();
    test_random();

    printf("%s\n", fails ? "FAILED" : "all passed");
    return fails ? 1 : 0;
}
