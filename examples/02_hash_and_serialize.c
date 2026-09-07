/*
 * Example 2 -- turning bytes into points, and points back into bytes.
 *
 * Two things a pairing library needs before it is usable in a protocol:
 *
 *   hash-to-curve   an arbitrary message becomes a point, uniformly, with no
 *                   way to work backwards. RFC 9380.
 *   serialization   a point becomes a fixed-length byte string that another
 *                   implementation can read back.
 *
 * The part people get wrong is the domain separation tag, so this example
 * spends most of its space on it.
 *
 * Run:  ./build/examples/elips_example_hash
 */
#include <stdio.h>
#include <string.h>

#include "elips/hash_to_curve.h"
#include "elips/serialize.h"
#include "elips/pairing.h"

/* The tag has two jobs: separate THIS protocol from every other one, and name
 * the suite so that two implementations agree on what they are computing.
 * ELIPS_H2C_SUITE_G1 supplies the second half; the first half is yours, and it
 * should include a version so you can change the scheme later without
 * colliding with your own old signatures. */
static const char DST_G1[] = "ELIPS-EXAMPLE-V01-CS01-" ELIPS_H2C_SUITE_G1;
static const char DST_G2[] = "ELIPS-EXAMPLE-V01-CS01-" ELIPS_H2C_SUITE_G2;

static void print_hex(const char *label, const uint8_t *b, size_t n)
{
    printf("  %-22s ", label);
    size_t show = n > 16 ? 16 : n;
    for (size_t i = 0; i < show; i++) printf("%02x", b[i]);
    if (n > show) printf("... (%zu bytes)", n);
    printf("\n");
}

int main(void)
{
    printf("== ELiPS hash-to-curve and serialization, curve %s ==\n\n",
           ELIPS_CURVE_NAME);
    printf("suite G1: %s\n", elips_h2c_suite_g1());
    printf("suite G2: %s\n\n", elips_h2c_suite_g2());

    const uint8_t msg[] = "the quick brown fox";
    const size_t  mlen  = sizeof msg - 1;

    /* --- hash to G1 ------------------------------------------------------ */
    ep_t P;
    if (elips_hash_to_g1(&P, msg, mlen,
                         (const uint8_t *)DST_G1, sizeof DST_G1 - 1) != 0) {
        fprintf(stderr, "hash_to_g1 refused the tag\n");
        return 1;
    }
    printf("hash_to_g1(\"%s\"):\n", msg);
    printf("  %-22s %s\n", "on the curve",   ep_on_curve(&P)    ? "yes" : "NO");
    printf("  %-22s %s\n", "in the subgroup", ep_in_subgroup(&P) ? "yes" : "NO");

    uint8_t c1[EP_SER_COMPRESSED_BYTES], u1[EP_SER_UNCOMPRESSED_BYTES];
    ep_write_compressed(c1, &P);
    ep_write_uncompressed(u1, &P);
    print_hex("compressed", c1, sizeof c1);
    print_hex("uncompressed", u1, sizeof u1);

    /* --- and back again -------------------------------------------------- */
    ep_t back;
    int rc = ep_read_compressed(&back, c1);
    printf("  %-22s %s\n", "decodes", elips_ser_strerror(rc));

    /* --- the tag is what separates protocols ----------------------------- */
    printf("\ndomain separation:\n");
    ep_t A, B;
    elips_hash_to_g1(&A, msg, mlen, (const uint8_t *)"PROTOCOL-A", 10);
    elips_hash_to_g1(&B, msg, mlen, (const uint8_t *)"PROTOCOL-B", 10);
    uint8_t ea[EP_SER_COMPRESSED_BYTES], eb[EP_SER_COMPRESSED_BYTES];
    ep_write_compressed(ea, &A);
    ep_write_compressed(eb, &B);
    print_hex("under PROTOCOL-A", ea, sizeof ea);
    print_hex("under PROTOCOL-B", eb, sizeof eb);
    printf("  the same message under two tags gives two unrelated points,\n"
           "  which is why a signature for one protocol cannot be replayed\n"
           "  into another.\n");

    printf("\n  an empty tag is refused:      %s\n",
           elips_hash_to_g1(&A, msg, mlen, (const uint8_t *)"", 0) == -1
               ? "yes" : "NO -- that is a bug");

    /* --- hash to G2 ------------------------------------------------------ */
    ep2_t Q;
    if (elips_hash_to_g2(&Q, msg, mlen,
                         (const uint8_t *)DST_G2, sizeof DST_G2 - 1) != 0) {
        fprintf(stderr, "hash_to_g2 refused the tag\n");
        return 1;
    }
    printf("\nhash_to_g2(\"%s\"):\n", msg);
    printf("  %-22s %s\n", "in the subgroup", ep2_in_subgroup(&Q) ? "yes" : "NO");
    uint8_t c2[EP2_SER_COMPRESSED_BYTES];
    ep2_write_compressed(c2, &Q);
    print_hex("compressed", c2, sizeof c2);

    /* --- the decoder is the trust boundary -------------------------------- */
    printf("\nwhat the decoder refuses, and what it says:\n");
    {
        uint8_t bad[EP_SER_COMPRESSED_BYTES];

        memcpy(bad, c1, sizeof bad);
        bad[0] &= (uint8_t)~0x80u;                     /* compression bit off */
        printf("  %-30s %s\n", "compression bit cleared",
               elips_ser_strerror(ep_read_compressed(&back, bad)));

        memcpy(bad, c1, sizeof bad);
        bad[0] |= 0x40u;                               /* infinity over a body */
        printf("  %-30s %s\n", "infinity flag on a real point",
               elips_ser_strerror(ep_read_compressed(&back, bad)));

        memset(bad, 0xff, sizeof bad);
        bad[0] = (uint8_t)(0x80u | (bad[0] & 0x1fu));  /* x far above p */
        printf("  %-30s %s\n", "coordinate not reduced",
               elips_ser_strerror(ep_read_compressed(&back, bad)));

        uint8_t badu[EP_SER_UNCOMPRESSED_BYTES];
        memcpy(badu, u1, sizeof badu);
        badu[sizeof badu - 1] ^= 1u;                   /* perturb y */
        printf("  %-30s %s\n", "point off the curve",
               elips_ser_strerror(ep_read_uncompressed(&back, badu)));
    }

    printf("\n  Every deserializer checks the flags, that coordinates are\n"
           "  reduced, that the point is on the curve, AND that it is in the\n"
           "  order-r subgroup. The last one is what stops small-subgroup\n"
           "  attacks, and it is the one implementations skip.\n");

#ifdef ELIPS_CURVE_BLS12_381
    printf("\n  On BLS12-381 these encodings are the specification's, so the\n"
           "  bytes above are what any conforming implementation produces.\n");
#else
    printf("\n  Note: %s has no standardised suite, so these encodings are\n"
           "  self-consistent but interoperate with nothing. Use BLS12-381 if\n"
           "  you need to talk to other software.\n", ELIPS_CURVE_NAME);
#endif
    return 0;
}
