/*
 * Hash-to-curve known-answer runner.
 *
 * Drives the library from test/kat/h2c_<curve>.vec, which the Python reference
 * in tools/reference/h2c_ref.py produced. For BLS12-381 those records use
 * RFC 9380's own domain separation tag and are therefore the specification's
 * published vectors: if this passes, a point this library hashes to is the same
 * point every conforming implementation hashes to.
 *
 * Each layer is checked separately -- SHA-256, then expand_message_xmd, then
 * hash_to_field, then the whole map -- because they all fail the same way, with
 * a wrong point, and a single end-to-end check would not say which one broke.
 *
 * Exit codes: 0 all passed, 1 a vector failed, 2 the file could not be read.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "elips/hash_to_curve.h"
#include "elips/h2c_params.h"
#include "elips/pairing.h"
#include "elips/sha256.h"

static long n_pass, n_fail;
static long cur_line;
static mpz_t PRIME;

static void fail(const char *op, const char *what)
{
    n_fail++;
    fprintf(stderr, "FAIL line %ld [%s]: %s\n", cur_line, op, what);
}

/* "-" encodes the empty string, which is a message the suite must handle and
 * which no other spelling survives a whitespace-separated format. */
static size_t unhex(uint8_t *out, size_t cap, const char *h)
{
    if (strcmp(h, "-") == 0) return 0;
    size_t n = strlen(h) / 2;
    if (n > cap) return (size_t)-1;
    for (size_t i = 0; i < n; i++) {
        unsigned v;
        if (sscanf(h + 2 * i, "%2x", &v) != 1) return (size_t)-1;
        out[i] = (uint8_t)v;
    }
    return n;
}

static void fp_to_mpz(mpz_t out, const fp_t a)
{
    limb_t l[FP_LIMBS];
    fp_to_limbs(l, a);
    mpz_import(out, FP_LIMBS, -1, sizeof(limb_t), 0, 0, l);
}

static int cmp_fp(const fp_t got, const char *want_hex, const char *op,
                  const char *what)
{
    mpz_t g, w;
    mpz_init(g); mpz_init_set_str(w, want_hex, 16);
    fp_to_mpz(g, got);
    int ok = mpz_cmp(g, w) == 0;
    if (!ok) {
        char *gs = mpz_get_str(NULL, 16, g);
        fail(op, what);
        fprintf(stderr, "    got  %s\n    want %s\n", gs, want_hex);
        free(gs);
    }
    mpz_clear(g); mpz_clear(w);
    return ok;
}

#define MAX_MSG 4096
#define MAX_OUT 4096

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <vectors.vec>\n", argv[0]); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }

    mpz_init(PRIME);
    printf("hash-to-curve vectors [%s]\n  suites: %s / %s\n",
           ELIPS_CURVE_NAME, elips_h2c_suite_g1(), elips_h2c_suite_g2());

    char *line = NULL; size_t cap = 0; ssize_t len;
    static uint8_t dst[MAX_MSG], msg[MAX_MSG], out[MAX_OUT];

    while ((len = getline(&line, &cap, f)) > 0) {
        cur_line++;
        if (line[0] == '#' || line[0] == '\n') continue;

        char *save = NULL;
        char *op = strtok_r(line, " \t\n", &save);
        if (!op) continue;

        if (strcmp(op, "prime") == 0) {
            char *h = strtok_r(NULL, " \t\n", &save);
            mpz_set_str(PRIME, h, 16);
            /* The vector file and the compiled curve have to be the same one;
             * every value below would otherwise be quietly wrong. */
            mpz_t mine; mpz_init(mine);
            mpz_import(mine, FP_LIMBS, -1, sizeof(limb_t), 0, 0, FP_MODULUS);
            if (mpz_cmp(mine, PRIME) != 0) {
                fprintf(stderr, "vector file is for a different curve than %s\n",
                        ELIPS_CURVE_NAME);
                mpz_clear(mine); free(line); fclose(f); mpz_clear(PRIME);
                return 2;
            }
            mpz_clear(mine);
            continue;
        }

        if (strcmp(op, "sha256") == 0) {
            char *mh = strtok_r(NULL, " \t\n", &save);
            strtok_r(NULL, " \t\n", &save);                  /* "=" */
            char *wh = strtok_r(NULL, " \t\n", &save);
            size_t ml = unhex(msg, MAX_MSG, mh);
            uint8_t dg[ELIPS_SHA256_DIGEST], want[ELIPS_SHA256_DIGEST];
            elips_sha256(dg, msg, ml);
            unhex(want, sizeof want, wh);
            if (memcmp(dg, want, sizeof dg) != 0) fail(op, "digest mismatch");
            else n_pass++;
            continue;
        }

        if (strcmp(op, "xmd") == 0) {
            char *dh = strtok_r(NULL, " \t\n", &save);
            char *mh = strtok_r(NULL, " \t\n", &save);
            long n = strtol(strtok_r(NULL, " \t\n", &save), NULL, 10);
            strtok_r(NULL, " \t\n", &save);                  /* "=" */
            char *wh = strtok_r(NULL, " \t\n", &save);
            size_t dl = unhex(dst, MAX_MSG, dh), ml = unhex(msg, MAX_MSG, mh);
            static uint8_t want[MAX_OUT];
            unhex(want, MAX_OUT, wh);
            if (elips_expand_message_xmd(out, (size_t)n, msg, ml, dst, dl) != 0)
                fail(op, "expand_message_xmd refused a valid length");
            else if (memcmp(out, want, (size_t)n) != 0)
                fail(op, "expanded bytes mismatch");
            else n_pass++;
            continue;
        }

        int is_g2 = strstr(op, "_g2") != NULL;
        char *dh = strtok_r(NULL, " \t\n", &save);
        char *mh = strtok_r(NULL, " \t\n", &save);
        char *eq = strtok_r(NULL, " \t\n", &save);
        if (!dh || !mh || !eq || strcmp(eq, "=")) { fail(op, "malformed record"); continue; }
        size_t dl = unhex(dst, MAX_MSG, dh), ml = unhex(msg, MAX_MSG, mh);
        if (dl == (size_t)-1 || ml == (size_t)-1) { fail(op, "hex too long"); continue; }

        char *w[4];
        int nw = is_g2 ? 4 : 2;
        int short_rec = 0;
        for (int i = 0; i < nw; i++) {
            w[i] = strtok_r(NULL, " \t\n", &save);
            if (!w[i]) short_rec = 1;
        }
        if (short_rec) { fail(op, "short record"); continue; }

        int ok = 1;
        if (strcmp(op, "h2f_g1") == 0) {
            fp_t u[2];
            elips_hash_to_field_fp(u, 2, msg, ml, dst, dl);
            ok &= cmp_fp(u[0], w[0], op, "u[0]");
            ok &= cmp_fp(u[1], w[1], op, "u[1]");
        } else if (strcmp(op, "h2f_g2") == 0) {
            fp2_t u[2];
            elips_hash_to_field_fp2(u, 2, msg, ml, dst, dl);
            ok &= cmp_fp(u[0][0], w[0], op, "u[0].c0");
            ok &= cmp_fp(u[0][1], w[1], op, "u[0].c1");
            ok &= cmp_fp(u[1][0], w[2], op, "u[1].c0");
            ok &= cmp_fp(u[1][1], w[3], op, "u[1].c1");
        } else if (strcmp(op, "h2c_g1") == 0 || strcmp(op, "enc_g1") == 0) {
            ep_t P; fp_t x, y;
            if (op[0] == 'h') elips_hash_to_g1(&P, msg, ml, dst, dl);
            else              elips_encode_to_g1(&P, msg, ml, dst, dl);
            if (!ep_to_affine(x, y, &P)) { fail(op, "result is the identity"); continue; }
            ok &= cmp_fp(x, w[0], op, "x");
            ok &= cmp_fp(y, w[1], op, "y");
            if (!ep_in_subgroup(&P)) { fail(op, "result is outside G1"); ok = 0; }
        } else if (strcmp(op, "h2c_g2") == 0 || strcmp(op, "enc_g2") == 0) {
            ep2_t P; fp2_t x, y;
            if (op[0] == 'h') elips_hash_to_g2(&P, msg, ml, dst, dl);
            else              elips_encode_to_g2(&P, msg, ml, dst, dl);
            if (!ep2_to_affine(x, y, &P)) { fail(op, "result is the identity"); continue; }
            ok &= cmp_fp(x[0], w[0], op, "x.c0");
            ok &= cmp_fp(x[1], w[1], op, "x.c1");
            ok &= cmp_fp(y[0], w[2], op, "y.c0");
            ok &= cmp_fp(y[1], w[3], op, "y.c1");
            if (!ep2_in_subgroup(&P)) { fail(op, "result is outside G2"); ok = 0; }
        } else {
            fprintf(stderr, "unknown op '%s' line %ld\n", op, cur_line);
            free(line); fclose(f); mpz_clear(PRIME);
            return 2;
        }
        if (ok) n_pass++;
    }

    free(line);
    fclose(f);
    mpz_clear(PRIME);

    printf("  %ld passed, %ld failed\n", n_pass, n_fail);
    return n_fail ? 1 : 0;
}
