/*
 * Example 3 -- a BLS signature, end to end.
 *
 * This is what the library is for. It uses every piece the earlier examples
 * introduced: a CSPRNG for the key, hash-to-curve for the message, the pairing
 * for verification, and the encodings for anything that goes on a wire.
 *
 *   secret key   sk, a uniform scalar in [0, r)
 *   public key   pk = [sk]G2                      96 bytes compressed
 *   signature    sigma = [sk] H(m), H into G1     48 bytes compressed
 *   verify       e(sigma, G2) == e(H(m), pk)
 *
 * The verification identity is just bilinearity read twice:
 *   e([sk]H(m), G2) = e(H(m), G2)^sk = e(H(m), [sk]G2).
 *
 * READ THIS BEFORE COPYING ANY OF IT.
 *
 * This is a demonstration of the library's API, not a signature implementation
 * you should deploy. What a real one needs and this does not have:
 *
 *   - The IETF ciphersuite from draft-irtf-cfrg-bls-signature, including its
 *     exact domain separation tag. The tag below is this example's own, so
 *     these signatures verify against nothing else.
 *   - Proof of possession, or the message augmentation scheme. Aggregating
 *     public keys without one is vulnerable to a rogue-key attack: an attacker
 *     who picks pk_evil = [t]G2 - sum(pk_honest) can forge an aggregate
 *     signature for a group it never had a key in. The aggregation section
 *     below is honest about this.
 *   - Key serialization, storage and zeroization.
 *   - A decision about what happens on a malformed public key from the network.
 *     Here it is checked; make sure yours is.
 *
 * Run:  ./build/examples/elips_example_bls
 */
#include <stdio.h>
#include <string.h>

#include "elips/hash_to_curve.h"
#include "elips/serialize.h"
#include "elips/pairing.h"
#include "elips/random.h"

/* This example's own tag. A real scheme uses the one its specification names. */
static const char DST[] = "ELIPS-EXAMPLE-BLS-SIG-V01-CS01-" ELIPS_H2C_SUITE_G1;
#define DST_LEN (sizeof DST - 1)

typedef struct { limb_t sk[ELIPS_ORDER_LIMBS]; ep2_t pk; } keypair_t;

static int keygen(keypair_t *kp)
{
    if (elips_random_scalar(kp->sk) != 0) return -1;
    /* A zero secret key gives the identity public key, which is useless and
     * which the verifier below would reject anyway. The probability is about
     * 2^-255, but "cannot happen" is not a reason to leave it unhandled. */
    if (!elips_scalar_is_reduced(kp->sk)) return -1;
    ep2_t G2;
    ep2_generator(&G2);
    ep2_mul(&kp->pk, &G2, kp->sk, ELIPS_ORDER_BITS);
    return ep2_is_infinity(&kp->pk) ? -1 : 0;
}

static int sign(ep_t *sig, const keypair_t *kp,
                const uint8_t *msg, size_t msg_len)
{
    ep_t h;
    if (elips_hash_to_g1(&h, msg, msg_len, (const uint8_t *)DST, DST_LEN) != 0)
        return -1;
    ep_mul(sig, &h, kp->sk, ELIPS_ORDER_BITS);
    return 0;
}

/* Returns 1 for a good signature, 0 for anything else. */
static int verify(const ep2_t *pk, const ep_t *sig,
                  const uint8_t *msg, size_t msg_len)
{
    /* Both of these come off a network in real use, so both are checked. The
     * pairing checks them again -- elips_pairing validates its inputs -- but a
     * verifier that relies on that is one refactor away from not doing it. */
    if (ep_is_infinity(sig) || !ep_in_subgroup(sig))  return 0;
    if (ep2_is_infinity(pk) || !ep2_in_subgroup(pk))  return 0;

    ep_t h;
    if (elips_hash_to_g1(&h, msg, msg_len, (const uint8_t *)DST, DST_LEN) != 0)
        return 0;

    ep2_t G2;
    ep2_generator(&G2);

    fp12_t lhs, rhs;
    if (!elips_pairing(lhs, sig, &G2)) return 0;
    if (!elips_pairing(rhs, &h, pk))   return 0;

    /* On BLS12 both sides come back cubed, and (x^3 == y^3) exactly when
     * (x == y) because gcd(3, r) = 1. The comparison is unaffected. */
    return fp12_eq(lhs, rhs);
}

static void show(const char *label, int good)
{
    printf("  %-52s %s\n", label, good ? "accepted" : "rejected");
}

int main(void)
{
    printf("== ELiPS BLS signature example, curve %s ==\n\n", ELIPS_CURVE_NAME);
    printf("tag: %s\n\n", DST);

    keypair_t alice;
    if (keygen(&alice) != 0) { fprintf(stderr, "key generation failed\n"); return 1; }

    const uint8_t msg[]  = "transfer 10 coins to bob";
    const size_t  mlen   = sizeof msg - 1;

    ep_t sig;
    if (sign(&sig, &alice, msg, mlen) != 0) { fprintf(stderr, "signing failed\n"); return 1; }

    /* What actually goes on the wire. */
    uint8_t sig_bytes[EP_SER_COMPRESSED_BYTES];
    uint8_t pk_bytes[EP2_SER_COMPRESSED_BYTES];
    ep_write_compressed(sig_bytes, &sig);
    ep2_write_compressed(pk_bytes, &alice.pk);
    printf("message    \"%s\"\n", msg);
    printf("public key %d bytes\n", EP2_SER_COMPRESSED_BYTES);
    printf("signature  %d bytes\n\n", EP_SER_COMPRESSED_BYTES);

    /* And the receiving side, which starts from bytes and nothing else. */
    ep_t  sig_in;
    ep2_t pk_in;
    int rc1 = ep_read_compressed(&sig_in, sig_bytes);
    int rc2 = ep2_read_compressed(&pk_in, pk_bytes);
    printf("decoding the wire format:\n");
    printf("  signature   %s\n", elips_ser_strerror(rc1));
    printf("  public key  %s\n\n", elips_ser_strerror(rc2));
    if (rc1 != ELIPS_SER_OK || rc2 != ELIPS_SER_OK) return 1;

    printf("verification:\n");
    show("the genuine signature", verify(&pk_in, &sig_in, msg, mlen));

    const uint8_t tampered[] = "transfer 10 coins to eve";
    show("the same signature, one word of the message changed",
         verify(&pk_in, &sig_in, tampered, sizeof tampered - 1));

    keypair_t mallory;
    if (keygen(&mallory) != 0) return 1;
    show("the right message, somebody else's public key",
         verify(&mallory.pk, &sig_in, msg, mlen));

    {   /* A signature that is a valid curve point but the wrong one. */
        ep_t forged;
        ep_neg(&forged, &sig_in);
        show("the signature negated", verify(&pk_in, &forged, msg, mlen));
    }

    {   /* The identity is a point, and it is not a signature. */
        ep_t zero;
        ep_set_infinity(&zero);
        show("the identity offered as a signature", verify(&pk_in, &zero, msg, mlen));
    }

    /* --- aggregation ------------------------------------------------------ */
    printf("\naggregation, three signers on one message:\n");
    keypair_t signers[3];
    ep_t sigs[3], agg_sig;
    ep2_t agg_pk;
    for (int i = 0; i < 3; i++) {
        if (keygen(&signers[i]) != 0) return 1;
        if (sign(&sigs[i], &signers[i], msg, mlen) != 0) return 1;
    }
    ep_copy(&agg_sig, &sigs[0]);
    ep2_copy(&agg_pk, &signers[0].pk);
    for (int i = 1; i < 3; i++) {
        ep_add(&agg_sig, &agg_sig, &sigs[i]);
        ep2_add(&agg_pk, &agg_pk, &signers[i].pk);
    }
    show("three signatures added, three keys added", verify(&agg_pk, &agg_sig, msg, mlen));

    {   /* Drop one signature but keep all three keys. */
        ep_t partial;
        ep_add(&partial, &sigs[0], &sigs[1]);
        show("only two of the three signatures", verify(&agg_pk, &partial, msg, mlen));
    }

    printf("\n  Aggregation collapses three signatures into %d bytes and three\n"
           "  keys into %d. It is also where BLS is easiest to get wrong: with\n"
           "  no proof of possession, an attacker who chooses its key AFTER\n"
           "  seeing the others can produce a valid aggregate for a group it\n"
           "  never joined. Any real deployment needs proof of possession or\n"
           "  message augmentation.\n",
           EP_SER_COMPRESSED_BYTES, EP2_SER_COMPRESSED_BYTES);

    printf("\n  Cost note: verification above computes two separate pairings,\n"
           "  so it runs the final exponentiation twice. A multi-pairing would\n"
           "  share one -- worth roughly 40%% of a verification, and the reason\n"
           "  it is a scheduled item in MODERNIZATION_PLAN.md section 10.6.\n");
    return 0;
}
