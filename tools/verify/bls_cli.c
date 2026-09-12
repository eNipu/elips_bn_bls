/*
 * A hex-in, hex-out shell around include/elips/bls.h, so a script can drive
 * the C implementation and compare it against another one.
 *
 * It exists for tools/verify/crosscheck_bls_pyecc.py. Keeping the comparison
 * in Python means the reference stays py_ecc, an implementation nobody here
 * wrote, rather than a second C routine that could share a mistake with the
 * first.
 *
 * Prints the result as hex and exits 0, or prints the error and exits 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elips/bls.h"

static int unhex(uint8_t *out, size_t cap, const char *h, size_t *n)
{
    size_t len = strlen(h);
    if (len % 2 || len / 2 > cap) return -1;
    for (size_t i = 0; i < len / 2; i++) {
        unsigned v;
        if (sscanf(h + 2 * i, "%2x", &v) != 1) return -1;
        out[i] = (uint8_t)v;
    }
    *n = len / 2;
    return 0;
}

static void phex(const uint8_t *b, size_t n)
{
    for (size_t i = 0; i < n; i++) printf("%02x", b[i]);
    printf("\n");
}

static int fail(int rc)
{
    fprintf(stderr, "%s\n", elips_bls_strerror(rc));
    printf("ERR %d\n", rc);
    return 1;
}

#define BUF 65536

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <op> [args]\n", argv[0]); return 2; }
    const char *op = argv[1];

    static uint8_t a[BUF], b[BUF], c[BUF];
    size_t na, nb, nc;
    uint8_t sk[ELIPS_BLS_SK_BYTES], pk[ELIPS_BLS_PK_BYTES], sig[ELIPS_BLS_SIG_BYTES];

    if (!strcmp(op, "keygen") && argc == 3) {
        if (unhex(a, sizeof a, argv[2], &na)) return 2;
        int rc = elips_bls_keygen(sk, a, na);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        phex(sk, sizeof sk); return 0;
    }
    if (!strcmp(op, "sk_to_pk") && argc == 3) {
        if (unhex(a, sizeof a, argv[2], &na) || na != sizeof sk) return 2;
        int rc = elips_bls_sk_to_pk(pk, a);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        phex(pk, sizeof pk); return 0;
    }
    if (!strcmp(op, "sign") && argc == 4) {
        if (unhex(a, sizeof a, argv[2], &na) || na != sizeof sk) return 2;
        if (unhex(b, sizeof b, argv[3], &nb)) return 2;
        int rc = elips_bls_sign(sig, a, b, nb);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        phex(sig, sizeof sig); return 0;
    }
    if (!strcmp(op, "verify") && argc == 5) {
        if (unhex(a, sizeof a, argv[2], &na)) return 2;
        if (unhex(b, sizeof b, argv[3], &nb)) return 2;
        if (unhex(c, sizeof c, argv[4], &nc)) return 2;
        if (na != sizeof pk || nc != sizeof sig) { printf("ERR -1\n"); return 1; }
        int rc = elips_bls_verify(a, b, nb, c);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        printf("OK\n"); return 0;
    }
    if (!strcmp(op, "aggregate") && argc >= 3) {
        size_t n = (size_t)(argc - 2);
        for (size_t i = 0; i < n; i++) {
            size_t got;
            if (unhex(a + i * ELIPS_BLS_SIG_BYTES, ELIPS_BLS_SIG_BYTES, argv[2 + i], &got)) return 2;
            if (got != ELIPS_BLS_SIG_BYTES) return 2;
        }
        int rc = elips_bls_aggregate(sig, a, n);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        phex(sig, sizeof sig); return 0;
    }
    /* aggregate_verify <sig> <pk1> <msg1> <pk2> <msg2> ... */
    if (!strcmp(op, "aggregate_verify") && argc >= 5 && (argc - 3) % 2 == 0) {
        size_t n = (size_t)(argc - 3) / 2;
        if (unhex(c, sizeof c, argv[2], &nc) || nc != sizeof sig) return 2;
        static uint8_t pks[64 * ELIPS_BLS_PK_BYTES];
        static uint8_t msgbuf[64][4096];
        const uint8_t *msgs[64];
        size_t lens[64];
        if (n > 64) return 2;
        for (size_t i = 0; i < n; i++) {
            size_t got;
            if (unhex(pks + i * ELIPS_BLS_PK_BYTES, ELIPS_BLS_PK_BYTES, argv[3 + 2 * i], &got)) return 2;
            if (got != ELIPS_BLS_PK_BYTES) return 2;
            if (unhex(msgbuf[i], sizeof msgbuf[i], argv[4 + 2 * i], &lens[i])) return 2;
            msgs[i] = msgbuf[i];
        }
        int rc = elips_bls_aggregate_verify(pks, n, msgs, lens, c);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        printf("OK\n"); return 0;
    }
    /* fast_aggregate_verify <sig> <msg> <pk1> <pop1> <pk2> <pop2> ... */
    if (!strcmp(op, "fast_aggregate_verify") && argc >= 6 && (argc - 4) % 2 == 0) {
        size_t n = (size_t)(argc - 4) / 2;
        if (unhex(c, sizeof c, argv[2], &nc) || nc != sizeof sig) return 2;
        if (unhex(b, sizeof b, argv[3], &nb)) return 2;
        static uint8_t pks[64 * ELIPS_BLS_PK_BYTES], pops[64 * ELIPS_BLS_POP_BYTES];
        if (n > 64) return 2;
        for (size_t i = 0; i < n; i++) {
            size_t g1, g2;
            if (unhex(pks + i * ELIPS_BLS_PK_BYTES, ELIPS_BLS_PK_BYTES, argv[4 + 2 * i], &g1)) return 2;
            if (unhex(pops + i * ELIPS_BLS_POP_BYTES, ELIPS_BLS_POP_BYTES, argv[5 + 2 * i], &g2)) return 2;
            if (g1 != ELIPS_BLS_PK_BYTES || g2 != ELIPS_BLS_POP_BYTES) return 2;
        }
        int rc = elips_bls_fast_aggregate_verify(pks, pops, n, b, nb, c);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        printf("OK\n"); return 0;
    }
    if (!strcmp(op, "pop_prove") && argc == 3) {
        if (unhex(a, sizeof a, argv[2], &na) || na != sizeof sk) return 2;
        int rc = elips_bls_pop_prove(sig, a);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        phex(sig, sizeof sig); return 0;
    }
    if (!strcmp(op, "pop_verify") && argc == 4) {
        if (unhex(a, sizeof a, argv[2], &na)) return 2;
        if (unhex(b, sizeof b, argv[3], &nb)) return 2;
        if (na != sizeof pk || nb != ELIPS_BLS_POP_BYTES) { printf("ERR -1\n"); return 1; }
        int rc = elips_bls_pop_verify(a, b);
        if (rc != ELIPS_BLS_OK) return fail(rc);
        printf("OK\n"); return 0;
    }
    fprintf(stderr, "unknown op or wrong argument count: %s\n", op);
    return 2;
}
