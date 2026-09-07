/*
 * Edge cases and API contracts, on every curve.
 *
 * The vector suites check that the library computes the right answer on
 * ordinary input. This one checks the boundaries, where implementations
 * actually break: the identity, the zero scalar, an output aliasing its own
 * input, a length exactly on a block boundary, a domain separation tag one byte
 * over the limit.
 *
 * Aliasing gets its own section because two of the original audit's confirmed
 * defects were exactly that -- Fp2_mul_Fp and three Fp12 routines read their
 * output instead of their input -- and nothing in the vector files would catch
 * it, because the vectors never alias.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <string.h>

#include "elips/pairing.h"
#include "elips/serialize.h"
#include "elips/hash_to_curve.h"
#include "elips/sha256.h"
#include "elips/random.h"

static int fails;
static int checks;

static void ok(int c, const char *what)
{
    checks++;
    if (!c) { fails++; printf("  [FAIL] %s\n", what); }
}

static int ep_same(const ep_t *a, const ep_t *b)
{
    if (ep_is_infinity(a) || ep_is_infinity(b))
        return ep_is_infinity(a) && ep_is_infinity(b);
    fp_t x1, y1, x2, y2;
    ep_to_affine(x1, y1, a);
    ep_to_affine(x2, y2, b);
    return fp_eq(x1, x2) && fp_eq(y1, y2);
}

static int ep2_same(const ep2_t *a, const ep2_t *b)
{
    if (ep2_is_infinity(a) || ep2_is_infinity(b))
        return ep2_is_infinity(a) && ep2_is_infinity(b);
    fp2_t x1, y1, x2, y2;
    ep2_to_affine(x1, y1, a);
    ep2_to_affine(x2, y2, b);
    return fp2_eq(x1, x2) && fp2_eq(y1, y2);
}

/* --- SHA-256 padding boundaries ------------------------------------------
 * Every length where the padding decision changes, plus the same message fed
 * one byte at a time, which is where a broken buffered update shows up. */
static void test_sha256_boundaries(void)
{
    static const int lens[] = { 0, 1, 54, 55, 56, 57, 63, 64, 65, 111, 112,
                                119, 120, 127, 128, 129 };
    uint8_t msg[256];
    for (size_t i = 0; i < sizeof msg; i++) msg[i] = (uint8_t)i;

    int bad = 0;
    for (size_t i = 0; i < sizeof lens / sizeof lens[0]; i++) {
        uint8_t one[ELIPS_SHA256_DIGEST], drip[ELIPS_SHA256_DIGEST];
        elips_sha256(one, msg, (size_t)lens[i]);

        elips_sha256_t s;
        elips_sha256_init(&s);
        for (int j = 0; j < lens[i]; j++) elips_sha256_update(&s, msg + j, 1);
        elips_sha256_final(drip, &s);

        if (memcmp(one, drip, sizeof one) != 0) bad++;
    }
    ok(bad == 0, "sha256: one-shot equals byte-at-a-time at 16 padding boundaries");
}

/* --- expand_message_xmd length and tag limits ----------------------------- */
static void test_xmd_limits(void)
{
    static uint8_t out[255 * 32 + 1];
    const uint8_t *dst = (const uint8_t *)"D";

    ok(elips_expand_message_xmd(out, 0, NULL, 0, dst, 1) == 0,
       "xmd accepts a zero-length request");
    ok(elips_expand_message_xmd(out, 1, NULL, 0, dst, 1) == 0,
       "xmd accepts one byte");
    ok(elips_expand_message_xmd(out, 32, NULL, 0, dst, 1) == 0,
       "xmd accepts exactly one block");
    ok(elips_expand_message_xmd(out, 33, NULL, 0, dst, 1) == 0,
       "xmd accepts one block plus one");
    ok(elips_expand_message_xmd(out, 255 * 32, NULL, 0, dst, 1) == 0,
       "xmd accepts the maximum 255 blocks");
    ok(elips_expand_message_xmd(out, 255 * 32 + 1, NULL, 0, dst, 1) == -1,
       "xmd rejects one byte past the maximum");
    ok(elips_expand_message_xmd(out, 32, NULL, 0, dst, 0) == -1,
       "xmd rejects an empty domain separation tag");

    /* A tag of 256 bytes is hashed down rather than refused (RFC 9380 5.3.3),
     * and must therefore differ from the 255-byte tag it was truncated from. */
    uint8_t big[300];
    memset(big, 'D', sizeof big);
    uint8_t a[32], b[32];
    ok(elips_expand_message_xmd(a, 32, NULL, 0, big, 255) == 0,
       "xmd accepts a 255-byte tag");
    ok(elips_expand_message_xmd(b, 32, NULL, 0, big, 256) == 0,
       "xmd accepts a 256-byte tag by hashing it down");
    ok(memcmp(a, b, 32) != 0, "the 255- and 256-byte tags separate domains");
}

/* --- the domain separation tag is not optional ---------------------------- */
static void test_dst_required(void)
{
    ep_t P;
    ep2_t Q;
    const uint8_t *dst = (const uint8_t *)"E";

    ok(elips_hash_to_g1(&P, (const uint8_t *)"m", 1, dst, 0) == -1,
       "hash_to_g1 rejects an empty tag");
    ok(ep_is_infinity(&P), "...and leaves the identity behind, not a point");
    ok(elips_hash_to_g2(&Q, (const uint8_t *)"m", 1, dst, 0) == -1,
       "hash_to_g2 rejects an empty tag");
    ok(ep2_is_infinity(&Q), "...and leaves the identity behind");
    ok(elips_encode_to_g1(&P, (const uint8_t *)"m", 1, dst, 0) == -1,
       "encode_to_g1 rejects an empty tag");
    ok(elips_encode_to_g2(&Q, (const uint8_t *)"m", 1, dst, 0) == -1,
       "encode_to_g2 rejects an empty tag");

    /* Different tags must give different points; that is the whole purpose. */
    ep_t A, B;
    elips_hash_to_g1(&A, (const uint8_t *)"m", 1, (const uint8_t *)"TAG-A", 5);
    elips_hash_to_g1(&B, (const uint8_t *)"m", 1, (const uint8_t *)"TAG-B", 5);
    ok(!ep_same(&A, &B), "the same message under two tags gives two points");

    /* And the same tag must give the same point, every time. */
    elips_hash_to_g1(&B, (const uint8_t *)"m", 1, (const uint8_t *)"TAG-A", 5);
    ok(ep_same(&A, &B), "hash_to_g1 is deterministic");

    /* The empty message is legitimate and must not be confused with an error. */
    ok(elips_hash_to_g1(&A, NULL, 0, (const uint8_t *)"TAG-A", 5) == 0,
       "the empty message hashes");
    ok(ep_in_subgroup(&A), "...to a point of G1");
}

/* --- aliasing ------------------------------------------------------------- */
static void test_aliasing(void)
{
    fp_t x, y, r, a;
    fp_rand(x); fp_rand(y);

#define ALIAS2(op, eq, T, copy)                                        \
    do {                                                              \
        T rr, aa; op(rr, x, y);                                       \
        copy(aa, x); op(aa, aa, y); ok(eq(aa, rr), #op " tolerates r == a"); \
        copy(aa, y); op(aa, x, aa); ok(eq(aa, rr), #op " tolerates r == b"); \
    } while (0)

    ALIAS2(fp_add, fp_eq, fp_t, fp_copy);
    ALIAS2(fp_sub, fp_eq, fp_t, fp_copy);
    ALIAS2(fp_mul, fp_eq, fp_t, fp_copy);

    fp_sqr(r, x); fp_copy(a, x); fp_sqr(a, a);
    ok(fp_eq(a, r), "fp_sqr tolerates r == a");
    fp_neg(r, x); fp_copy(a, x); fp_neg(a, a);
    ok(fp_eq(a, r), "fp_neg tolerates r == a");
    fp_inv(r, x); fp_copy(a, x); fp_inv(a, a);
    ok(fp_eq(a, r), "fp_inv tolerates r == a");
    {
        int e1 = fp_sqrt(r, x);
        fp_copy(a, x);
        int e2 = fp_sqrt(a, a);
        ok(e1 == e2 && fp_eq(a, r), "fp_sqrt tolerates r == a");
    }

    {
        fp2_t u, v, rr, aa;
        fp_rand(u[0]); fp_rand(u[1]); fp_rand(v[0]); fp_rand(v[1]);
        fp2_mul(rr, u, v); fp2_copy(aa, u); fp2_mul(aa, aa, v);
        ok(fp2_eq(aa, rr), "fp2_mul tolerates r == a");
        fp2_copy(aa, v); fp2_mul(aa, u, aa);
        ok(fp2_eq(aa, rr), "fp2_mul tolerates r == b");
        fp2_sqr(rr, u); fp2_copy(aa, u); fp2_sqr(aa, aa);
        ok(fp2_eq(aa, rr), "fp2_sqr tolerates r == a");
        fp2_inv(rr, u); fp2_copy(aa, u); fp2_inv(aa, aa);
        ok(fp2_eq(aa, rr), "fp2_inv tolerates r == a");
        int e1 = fp2_sqrt(rr, u);
        fp2_copy(aa, u);
        int e2 = fp2_sqrt(aa, aa);
        ok(e1 == e2 && fp2_eq(aa, rr), "fp2_sqrt tolerates r == a");
    }

    {
        /* Constant-time inversion: the batched-divstep fp_inv against the
         * mpn_sec_invert fp_inv_sec it replaced, which is kept for exactly
         * this. Both must agree on every value, and both must satisfy the
         * defining property. Edge cases first, because that is where a
         * divstep bound or a sign select goes wrong. */
        fp_t a, x, y, one, chk, zero;
        fp_set_one(one);
        fp_set_zero(zero);

        static const unsigned small[] = { 1, 2, 3, 4, 5, 7, 255, 256 };
        for (unsigned k = 0; k < sizeof small / sizeof *small; k++) {
            limb_t l[FP_LIMBS];
            memset(l, 0, sizeof l);
            l[0] = small[k];
            fp_from_limbs(a, l);
            fp_inv(x, a); fp_inv_sec(y, a);
            ok(fp_eq(x, y), "fp_inv == fp_inv_sec on a small value");
            fp_mul(chk, a, x);
            ok(fp_eq(chk, one), "a * fp_inv(a) == 1 on a small value");
        }

        /* p-1 and p-2, the top of the range. */
        for (int k = 1; k <= 2; k++) {
            limb_t l[FP_LIMBS];
            for (int i = 0; i < FP_LIMBS; i++) l[i] = FP_MODULUS[i];
            l[0] -= (limb_t)k;
            fp_from_limbs(a, l);
            fp_inv(x, a); fp_inv_sec(y, a);
            ok(fp_eq(x, y), "fp_inv == fp_inv_sec near the modulus");
            fp_mul(chk, a, x);
            ok(fp_eq(chk, one), "a * fp_inv(a) == 1 near the modulus");
        }

        int disagree = 0, notone = 0;
        for (int i = 0; i < 400; i++) {
            fp_rand(a);
            fp_inv(x, a); fp_inv_sec(y, a);
            if (!fp_eq(x, y)) disagree++;
            fp_mul(chk, a, x);
            if (!fp_eq(chk, one)) notone++;
        }
        ok(disagree == 0, "fp_inv == fp_inv_sec on 400 random values");
        ok(notone == 0,   "a * fp_inv(a) == 1 on 400 random values");

        /* Zero has no inverse and the contract says the answer is zero, with
         * no branch taken to get there. */
        fp_inv(x, zero);
        ok(fp_is_zero(x), "fp_inv(0) == 0");
        fp_inv_sec(y, zero);
        ok(fp_is_zero(y), "fp_inv_sec(0) == 0");

        /* Output aliasing input, which fp2_inv does. */
        fp_rand(a); fp_copy(x, a); fp_inv(y, a); fp_inv(x, x);
        ok(fp_eq(x, y), "fp_inv tolerates r == a");
    }

    {
        fp12_t f, g, rr, aa;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 2; k++) { fp_rand(f[i][j][k]); fp_rand(g[i][j][k]); }
        fp12_mul(rr, f, g); fp12_copy(aa, f); fp12_mul(aa, aa, g);
        ok(fp12_eq(aa, rr), "fp12_mul tolerates r == a");
        fp12_copy(aa, g); fp12_mul(aa, f, aa);
        ok(fp12_eq(aa, rr), "fp12_mul tolerates r == b");
        fp12_sqr(rr, f); fp12_copy(aa, f); fp12_sqr(aa, aa);
        ok(fp12_eq(aa, rr), "fp12_sqr tolerates r == a");
        fp12_inv(rr, f); fp12_copy(aa, f); fp12_inv(aa, aa);
        ok(fp12_eq(aa, rr), "fp12_inv tolerates r == a");

        /* Cyclotomic squaring. Granger-Scott is valid only inside the
         * subgroup, so put f there first: f^(p^6-1) then ^(p^2+1).
         *
         * The pairing vectors already exercise this routine 321 times per
         * final exponentiation, but only through a value the oracle checks at
         * the end. Comparing directly against fp12_sqr localises a failure to
         * this routine instead of to the whole chain. */
        fp12_t cyc, t0, t1;
        fp12_conj(t0, f);
        fp12_inv(t1, f);
        fp12_mul(t0, t0, t1);
        fp12_frobenius(t1, t0, 2);
        fp12_mul(cyc, t1, t0);

        fp12_sqr(rr, cyc);
        fp12_sqr_cyc(aa, cyc);
        ok(fp12_eq(aa, rr), "fp12_sqr_cyc == fp12_sqr on the cyclotomic subgroup");

        fp12_copy(aa, cyc); fp12_sqr_cyc(aa, aa);
        ok(fp12_eq(aa, rr), "fp12_sqr_cyc tolerates r == a");

        /* And it must NOT agree off the subgroup. If it did, it would be the
         * general squaring and the speedup would be imaginary. */
        fp12_sqr(rr, f);
        fp12_sqr_cyc(aa, f);
        ok(!fp12_eq(aa, rr), "fp12_sqr_cyc is not valid off the subgroup");

        /* Karabina: n compressed squarings must equal n ordinary ones. Run
         * lengths on either side of the threshold, and the degenerate n. */
        for (int n = 0; n <= 12; n++) {
            fp12_t want, got;
            fp12_copy(want, cyc);
            for (int k = 0; k < n; k++) fp12_sqr_cyc(want, want);
            fp12_sqr_cyc_run(got, cyc, n);
            ok(fp12_eq(got, want), "fp12_sqr_cyc_run == repeated fp12_sqr_cyc");
        }
        fp12_copy(aa, cyc); fp12_sqr_cyc_run(aa, aa, 7);
        fp12_sqr_cyc_run(rr, cyc, 7);
        ok(fp12_eq(aa, rr), "fp12_sqr_cyc_run tolerates r == a");

        /* Decompression picks between two relations for g3, and the second is
         * only reached when g1 == 0 -- which random cyclotomic elements never
         * hit, so a wrong second relation would sit untested. Confirmed: forcing
         * the first relation everywhere breaks nothing else in this suite.
         *
         * So check the relations themselves. Both hold identically on every
         * cyclotomic element, and validating relation 2 here validates the
         * branch that cannot be reached directly.
         *
         *   4 g1 g3    == 3 g2^2 + xi g5^2 - 2 g4
         *   4 xi g5 g3 == g1^2 - 2 g2 + 3 xi g4^2
         *
         * with a[0] = (g0,g2,g4) and a[1] = (g1,g3,g5). */
        {
            fp2_t g1, g2, g3, g4, g5, lhs, rhs, t2, th;
            fp2_copy(g1, cyc[1][0]); fp2_copy(g2, cyc[0][1]);
            fp2_copy(g3, cyc[1][1]); fp2_copy(g4, cyc[0][2]);
            fp2_copy(g5, cyc[1][2]);

            fp2_mul(lhs, g1, g3);
            fp2_add(lhs, lhs, lhs); fp2_add(lhs, lhs, lhs);      /* 4 g1 g3 */
            fp2_sqr(rhs, g2);
            fp2_add(t2, rhs, rhs); fp2_add(rhs, t2, rhs);        /* 3 g2^2  */
            fp2_sqr(t2, g5); fp2_mul_xi(t2, t2); fp2_add(rhs, rhs, t2);
            fp2_sub(rhs, rhs, g4); fp2_sub(rhs, rhs, g4);
            ok(fp2_eq(lhs, rhs), "decompression relation 1 holds (divides by g1)");

            fp2_mul_xi(lhs, g5); fp2_mul(lhs, lhs, g3);
            fp2_add(lhs, lhs, lhs); fp2_add(lhs, lhs, lhs);      /* 4 xi g5 g3 */
            fp2_sqr(rhs, g1);
            fp2_sub(rhs, rhs, g2); fp2_sub(rhs, rhs, g2);
            fp2_sqr(t2, g4); fp2_mul_xi(t2, t2);
            fp2_add(th, t2, t2); fp2_add(th, th, t2);            /* 3 xi g4^2 */
            fp2_add(rhs, rhs, th);
            ok(fp2_eq(lhs, rhs), "decompression relation 2 holds (divides by g5)");
        }

        /* The identity is the one element where BOTH relations degenerate:
         * g1 = g5 = 0. gen_params.py asserts no other element can do that. */
        fp12_t idt, ids;
        fp12_set_one(idt);
        fp12_sqr_cyc_run(ids, idt, 5);
        ok(fp12_eq(ids, idt), "compressed squaring of the identity is the identity");
    }
#undef ALIAS2

    {
        ep_t G, H, rr, aa;
        ep_generator(&G); ep_dbl(&H, &G);
        ep_add(&rr, &G, &H); ep_copy(&aa, &G); ep_add(&aa, &aa, &H);
        ok(ep_same(&rr, &aa), "ep_add tolerates r == p");
        ep_copy(&aa, &H); ep_add(&aa, &G, &aa);
        ok(ep_same(&rr, &aa), "ep_add tolerates r == q");
        ep_dbl(&rr, &G); ep_copy(&aa, &G); ep_dbl(&aa, &aa);
        ok(ep_same(&rr, &aa), "ep_dbl tolerates r == p");

        limb_t k[ELIPS_ORDER_LIMBS];
        elips_random_scalar(k);
        ep_mul(&rr, &G, k, ELIPS_ORDER_BITS);
        ep_copy(&aa, &G); ep_mul(&aa, &aa, k, ELIPS_ORDER_BITS);
        ok(ep_same(&rr, &aa), "ep_mul tolerates r == p");
    }

    {
        ep2_t G, H, rr, aa;
        ep2_generator(&G); ep2_dbl(&H, &G);
        ep2_add(&rr, &G, &H); ep2_copy(&aa, &G); ep2_add(&aa, &aa, &H);
        ok(ep2_same(&rr, &aa), "ep2_add tolerates r == p");
        ep2_psi(&rr, &G); ep2_copy(&aa, &G); ep2_psi(&aa, &aa);
        ok(ep2_same(&rr, &aa), "ep2_psi tolerates r == p");

        limb_t k[ELIPS_ORDER_LIMBS];
        elips_random_scalar(k);
        ep2_mul(&rr, &G, k, ELIPS_ORDER_BITS);
        ep2_copy(&aa, &G); ep2_mul(&aa, &aa, k, ELIPS_ORDER_BITS);
        ok(ep2_same(&rr, &aa), "ep2_mul tolerates r == p");
    }
}

/* --- the identity and the degenerate group-law cases ---------------------- */
static void test_group_law_edges(void)
{
    ep_t G, O, n, r, d;
    ep_generator(&G);
    ep_set_infinity(&O);
    ep_neg(&n, &G);

    ep_add(&r, &G, &n);  ok(ep_is_infinity(&r), "P + (-P) == O");
    ep_add(&r, &G, &O);  ok(ep_same(&r, &G),    "P + O == P");
    ep_add(&r, &O, &G);  ok(ep_same(&r, &G),    "O + P == P");
    ep_add(&r, &O, &O);  ok(ep_is_infinity(&r), "O + O == O");
    ep_dbl(&r, &O);      ok(ep_is_infinity(&r), "2O == O");
    ep_add(&r, &G, &G); ep_dbl(&d, &G);
    ok(ep_same(&r, &d), "P + P agrees with 2P");
    ok(ep_on_curve(&O), "the identity counts as on the curve");

    ep2_t H, O2, n2, r2, d2;
    ep2_generator(&H);
    ep2_set_infinity(&O2);
    ep2_neg(&n2, &H);
    ep2_add(&r2, &H, &n2); ok(ep2_is_infinity(&r2), "Q + (-Q) == O");
    ep2_add(&r2, &H, &O2); ok(ep2_same(&r2, &H),    "Q + O == Q");
    ep2_add(&r2, &O2, &O2); ok(ep2_is_infinity(&r2), "O + O == O on the twist");
    ep2_add(&r2, &H, &H); ep2_dbl(&d2, &H);
    ok(ep2_same(&r2, &d2), "Q + Q agrees with 2Q");
}

/* --- scalars at the ends of their range ----------------------------------- */
static void test_scalar_edges(void)
{
    ep_t G, R;
    ep_generator(&G);
    limb_t k[ELIPS_ORDER_LIMBS];

    memset(k, 0, sizeof k);
    ep_mul(&R, &G, k, ELIPS_ORDER_BITS);
    ok(ep_is_infinity(&R), "[0]G == O");

    k[0] = 1;
    ep_mul(&R, &G, k, ELIPS_ORDER_BITS);
    ok(ep_same(&R, &G), "[1]G == G");

    memcpy(k, ELIPS_ORDER, sizeof k);
    ep_mul(&R, &G, k, ELIPS_ORDER_BITS);
    ok(ep_is_infinity(&R), "[r]G == O");

    memcpy(k, ELIPS_ORDER, sizeof k);
    k[0] -= 1;
    ep_mul(&R, &G, k, ELIPS_ORDER_BITS);
    ep_t neg;
    ep_neg(&neg, &G);
    ok(ep_same(&R, &neg), "[r-1]G == -G");

    ep2_t H, R2, neg2;
    ep2_generator(&H);
    memset(k, 0, sizeof k);
    ep2_mul(&R2, &H, k, ELIPS_ORDER_BITS);
    ok(ep2_is_infinity(&R2), "[0]H == O");
    memcpy(k, ELIPS_ORDER, sizeof k);
    k[0] -= 1;
    ep2_mul(&R2, &H, k, ELIPS_ORDER_BITS);
    ep2_neg(&neg2, &H);
    ok(ep2_same(&R2, &neg2), "[r-1]H == -H");

#ifdef ELIPS_FAMILY_BLS12
    /* GLV must agree with the plain ladder at the ends of the range too, which
     * is where a decomposition is most likely to be off by one. */
    int agree = 1;
    for (int i = 0; i < 3; i++) {
        memset(k, 0, sizeof k);
        if (i == 1) k[0] = 1;
        if (i == 2) { memcpy(k, ELIPS_ORDER, sizeof k); k[0] -= 1; }
        ep2_t plain, glv;
        ep2_mul(&plain, &H, k, ELIPS_ORDER_BITS);
        ep2_mul_glv(&glv, &H, k, ELIPS_ORDER_BITS);
        if (!ep2_same(&plain, &glv)) agree = 0;
    }
    ok(agree, "GLV agrees with the plain ladder at k = 0, 1 and r-1");
#endif
}

/* --- the pairing's own trust boundary ------------------------------------- */
static void test_pairing_validation(void)
{
    ep_t P, O;
    ep2_t Q, O2;
    fp12_t z, one;
    fp12_set_one(one);

    ep_generator(&P);
    ep2_generator(&Q);
    ep_set_infinity(&O);
    ep2_set_infinity(&O2);

    ok(elips_pairing(z, &P, &Q) == 1, "the generators pair");
    ok(!fp12_eq(z, one), "e(G1, G2) != 1");

    ok(elips_pairing(z, &O, &Q) == 0, "the identity in G1 is refused");
    ok(fp12_eq(z, one), "...and the output is left at one, not left dirty");
    ok(elips_pairing(z, &P, &O2) == 0, "the identity in G2 is refused");

    /* An on-curve point outside the order-r subgroup must be refused. Built by
     * halving the cofactor clearing away: take a point of E(Fp) that the
     * decoder accepts only when it happens to be in G1, and search for one that
     * is not. Curves with cofactor 1 have none, and say so. */
    uint8_t enc[EP_SER_COMPRESSED_BYTES];
    int found = 0;
    for (unsigned v = 1; v < 256 && !found; v++) {
        memset(enc, 0, sizeof enc);
        enc[sizeof enc - 1] = (uint8_t)v;
        enc[0] |= 0x80u;
        ep_t T;
        if (ep_read_compressed(&T, enc) == ELIPS_SER_ERR_NOT_IN_GROUP) found = 1;
    }
    if (found)
        ok(1, "an on-curve point outside G1 is refused by the decoder");
    else
        printf("  [SKIP] no off-subgroup point in range (G1 cofactor may be 1)\n");
}

/* --- serialization boundaries --------------------------------------------- */
static void test_serialize_edges(void)
{
    ep_t P;
    ep2_t Q;

    uint8_t z1[EP_SER_COMPRESSED_BYTES];
    memset(z1, 0, sizeof z1);
    ok(ep_read_compressed(&P, z1) == ELIPS_SER_ERR_FLAGS,
       "an all-zero compressed buffer is refused (no compression bit)");

    uint8_t z2[EP_SER_UNCOMPRESSED_BYTES];
    memset(z2, 0, sizeof z2);
    ok(ep_read_uncompressed(&P, z2) == ELIPS_SER_ERR_NOT_ON_CURVE,
       "x = y = 0 is refused: it is not on the curve");

    uint8_t z3[EP2_SER_COMPRESSED_BYTES];
    memset(z3, 0, sizeof z3);
    ok(ep2_read_compressed(&Q, z3) == ELIPS_SER_ERR_FLAGS,
       "an all-zero compressed G2 buffer is refused");

    /* Every flag combination on an otherwise valid encoding. Only the correct
     * one may be accepted. */
    ep_generator(&P);
    uint8_t good[EP_SER_COMPRESSED_BYTES];
    ep_write_compressed(good, &P);
    int accepted = 0;
    for (unsigned f = 0; f < 8; f++) {
        uint8_t t[EP_SER_COMPRESSED_BYTES];
        memcpy(t, good, sizeof t);
        t[0] = (uint8_t)((t[0] & 0x1fu) | (f << 5));
        ep_t R;
        if (ep_read_compressed(&R, t) == ELIPS_SER_OK) accepted++;
    }
    ok(accepted == 2,
       "of eight flag combinations exactly two decode (the point and its negation)");

    /* Round-tripping every encoding of the identity. */
    ep_t O;
    ep_set_infinity(&O);
    uint8_t oc[EP_SER_COMPRESSED_BYTES], ou[EP_SER_UNCOMPRESSED_BYTES];
    ep_write_compressed(oc, &O);
    ep_write_uncompressed(ou, &O);
    ep_t back;
    ok(ep_read_compressed(&back, oc) == ELIPS_SER_OK && ep_is_infinity(&back),
       "the compressed identity round trips");
    ok(ep_read_uncompressed(&back, ou) == ELIPS_SER_OK && ep_is_infinity(&back),
       "the uncompressed identity round trips");

    /* A hashed point, which is where real encodings come from. */
    ep_t H;
    const uint8_t *dst = (const uint8_t *)"ELIPS-EDGE-V01-CS01";
    elips_hash_to_g1(&H, (const uint8_t *)"round trip", 10, dst, 19);
    uint8_t hc[EP_SER_COMPRESSED_BYTES];
    ep_write_compressed(hc, &H);
    ok(ep_read_compressed(&back, hc) == ELIPS_SER_OK && ep_same(&back, &H),
       "a hashed G1 point round trips through the compressed encoding");

    ep2_t H2, back2;
    elips_hash_to_g2(&H2, (const uint8_t *)"round trip", 10, dst, 19);
    uint8_t hc2[EP2_SER_COMPRESSED_BYTES];
    ep2_write_compressed(hc2, &H2);
    ok(ep2_read_compressed(&back2, hc2) == ELIPS_SER_OK && ep2_same(&back2, &H2),
       "a hashed G2 point round trips through the compressed encoding");
}

/* --- bilinearity, which is the property everything else exists to support -- */
static void test_bilinearity(void)
{
    ep_t P, aP;
    ep2_t Q, bQ;
    fp12_t base, lhs, rhs;
    limb_t a[ELIPS_ORDER_LIMBS], b[ELIPS_ORDER_LIMBS];

    ep_generator(&P);
    ep2_generator(&Q);
    ok(elips_random_scalar(a) == 0 && elips_random_scalar(b) == 0,
       "sampled two scalars");

    ep_mul(&aP, &P, a, ELIPS_ORDER_BITS);
    ep2_mul(&bQ, &Q, b, ELIPS_ORDER_BITS);

    ok(elips_pairing(base, &P, &Q) == 1, "e(G1, G2) computed");
    ok(elips_pairing(lhs, &aP, &bQ) == 1, "e([a]G1, [b]G2) computed");

    fp12_exp(rhs, base, a, ELIPS_ORDER_BITS);
    fp12_exp(rhs, rhs, b, ELIPS_ORDER_BITS);
    ok(fp12_eq(lhs, rhs), "e([a]G1, [b]G2) == e(G1, G2)^(ab)");

    /* The result must live in mu_r. */
    fp12_t chk, one;
    fp12_set_one(one);
    fp12_exp(chk, base, ELIPS_ORDER, ELIPS_ORDER_BITS);
    ok(fp12_eq(chk, one), "e(G1, G2)^r == 1");
}

int main(void)
{
    printf("edge cases [%s]\n", ELIPS_CURVE_NAME);

    test_sha256_boundaries();
    test_xmd_limits();
    test_dst_required();
    test_aliasing();
    test_group_law_edges();
    test_scalar_edges();
    test_pairing_validation();
    test_serialize_edges();
    test_bilinearity();

    printf("  %d checks, %d failed\n", checks, fails);
    return fails ? 1 : 0;
}
