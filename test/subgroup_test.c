/*
 * Do the fast subgroup tests accept exactly the right points?
 *
 * The direction that matters is rejection. A test that is too strict fails the
 * first time a valid point is checked, so ordinary use finds it. A test that is
 * too permissive fails silently: it accepts attacker-chosen points of small
 * order and every other test in the suite still passes. The old tests were
 * [r]P, which cannot be wrong in that direction; the endomorphism tests can, so
 * this runner exists.
 *
 * test/kat/subgroup_<curve>.vec carries both classes, built by the Python
 * oracle from the real group structure: ing1/ing2 are of order r and must be
 * accepted, outg1/outg2 are on the curve but in the cofactor part and must be
 * rejected.
 *
 * Exit codes: 0 all passed, 1 a check failed, 2 the vector file is unusable.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "elips/pairing.h"

static int fails, checks, accepted, rejected;
static long cur_line;

static void ok(int c, const char *what)
{
    checks++;
    if (!c) { fails++; printf("  [FAIL] %s (line %ld)\n", what, cur_line); }
}

static void ld(fp_t r, mpz_t tmp, const char *hex)
{
    limb_t l[FP_LIMBS];
    memset(l, 0, sizeof l);
    mpz_set_str(tmp, hex, 16);
    mpz_export(l, NULL, -1, sizeof(limb_t), 0, 0, tmp);
    fp_from_limbs(r, l);
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <subgroup.vec>\n", argv[0]); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }

    printf("subgroup vectors [%s]\n", ELIPS_CURVE_NAME);

    mpz_t tmp, prime;
    mpz_inits(tmp, prime, NULL);

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    long records = 0;

#define BAIL(msg) do { fprintf(stderr, msg " at line %ld\n", cur_line); \
                       free(line); fclose(f); \
                       mpz_clears(tmp, prime, NULL); return 2; } while (0)

    while ((len = getline(&line, &cap, f)) > 0) {
        cur_line++;
        if (line[0] == '#' || line[0] == '\n') continue;

        char *save = NULL;
        char *op = strtok_r(line, " \t\n", &save);
        if (!op) continue;

        if (strcmp(op, "prime") == 0) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h) BAIL("malformed prime");
            mpz_set_str(prime, h, 16);
            mpz_import(tmp, FP_LIMBS, -1, sizeof(limb_t), 0, 0, FP_MODULUS);
            if (mpz_cmp(tmp, prime) != 0) {
                fprintf(stderr, "vector file is for a different curve than %s\n",
                        ELIPS_CURVE_NAME);
                free(line); fclose(f); mpz_clears(tmp, prime, NULL); return 2;
            }
            continue;
        }

        int g2   = (strncmp(op, "ing2", 4) == 0) || (strncmp(op, "outg2", 5) == 0);
        int want = (strncmp(op, "in", 2) == 0);
        if (!want && strncmp(op, "out", 3) != 0) BAIL("unknown op");

        char *tok[4];
        int need = g2 ? 4 : 2;
        for (int i = 0; i < need; i++) {
            tok[i] = strtok_r(NULL, " \t\n", &save);
            if (!tok[i]) BAIL("malformed record");
        }

        if (g2) {
            fp2_t qx, qy;
            ld(qx[0], tmp, tok[0]); ld(qx[1], tmp, tok[1]);
            ld(qy[0], tmp, tok[2]); ld(qy[1], tmp, tok[3]);
            ep2_t Q;
            ep2_from_affine(&Q, qx, qy);
            /* The file's own claim about the point must hold first, or the
             * membership answer below is being compared against nothing. */
            ok(ep2_on_curve(&Q), "vector G2 point is on the twist");
            ok(ep2_in_subgroup(&Q) == want,
               want ? "order-r G2 point is accepted"
                    : "cofactor G2 point is REJECTED");
        } else {
            fp_t px, py;
            ld(px, tmp, tok[0]); ld(py, tmp, tok[1]);
            ep_t P;
            ep_from_affine(&P, px, py);
            ok(ep_on_curve(&P), "vector G1 point is on the curve");
            ok(ep_in_subgroup(&P) == want,
               want ? "order-r G1 point is accepted"
                    : "cofactor G1 point is REJECTED");
        }
        if (want) accepted++; else rejected++;
        records++;
    }

    free(line);
    fclose(f);
    mpz_clears(tmp, prime, NULL);

    if (records == 0) {
        fprintf(stderr, "no subgroup records in %s\n", argv[1]);
        return 2;
    }
    /* A file with nothing to reject would make this runner vacuous. Every
     * curve has G2 cofactor points; only BN's G1 cofactor is 1. */
    if (rejected == 0) {
        fprintf(stderr, "vector file carries no points to reject\n");
        return 2;
    }

    printf("  %ld points (%d to accept, %d to reject), %d checks, %d failed\n",
           records, accepted, rejected, checks, fails);
    return fails ? 1 : 0;
}
