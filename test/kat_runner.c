/*
 * Runs the Phase 0 known-answer vectors against the NEW Montgomery tower.
 *
 * Same .vec files that validated the legacy mpz layer, so the new arithmetic is
 * checked against the same independent Python reference rather than against the
 * code it replaces. If both layers pass the same vectors, they agree with the
 * oracle, which is a stronger statement than them agreeing with each other.
 *
 * Elliptic-curve records are handled by the EC section further down once that
 * layer exists; until then they are counted as skipped and reported, so a
 * silently unverified operation cannot hide.
 */
/* getline, strtok_r and ssize_t are POSIX, and the build compiles with
 * -std=c11 (no GNU extensions), which hides them on glibc unless asked for. */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "elips/fpx.h"
#include "elips/ec.h"

static mpz_t PRIME;
static long n_pass, n_fail, n_skip;
static char cur_op[64];
static long cur_line;

static void load_fp(fp_t r, const mpz_t v)
{
    limb_t l[FP_LIMBS];
    memset(l, 0, sizeof l);
    mpz_export(l, NULL, -1, sizeof(limb_t), 0, 0, v);
    fp_from_limbs(r, l);
}
static void store_fp(mpz_t out, const fp_t a)
{
    limb_t l[FP_LIMBS];
    fp_to_limbs(l, a);
    mpz_import(out, FP_LIMBS, -1, sizeof(limb_t), 0, 0, l);
}

static void expect(const mpz_t got, const mpz_t want, int idx)
{
    if (mpz_cmp(got, want) != 0) {
        n_fail++;
        if (n_fail <= 10)
            gmp_fprintf(stderr, "MISMATCH line %ld op %s token %d\n"
                                "  expected %Zx\n  actual   %Zx\n",
                        cur_line, cur_op, idx, want, got);
    }
}

/* Generic driver: read `in_n` field elements, run `fn`, compare `out_n`. */
typedef void (*op_fn)(mpz_t *in, mpz_t *out);
typedef struct { const char *name; int in_n, out_n; op_fn run; } Op;

#define LOADN(dst, src, n)  for (int _i = 0; _i < (n); _i++) load_fp(dst[_i], src[_i])
#define STOREN(dst, src, n) for (int _i = 0; _i < (n); _i++) store_fp(dst[_i], src[_i])

/* --- fp --- */
static void o_fp_add(mpz_t *i, mpz_t *o){ fp_t a,b,r; load_fp(a,i[0]); load_fp(b,i[1]); fp_add(r,a,b); store_fp(o[0],r); }
static void o_fp_sub(mpz_t *i, mpz_t *o){ fp_t a,b,r; load_fp(a,i[0]); load_fp(b,i[1]); fp_sub(r,a,b); store_fp(o[0],r); }
static void o_fp_mul(mpz_t *i, mpz_t *o){ fp_t a,b,r; load_fp(a,i[0]); load_fp(b,i[1]); fp_mul(r,a,b); store_fp(o[0],r); }
static void o_fp_inv(mpz_t *i, mpz_t *o){ fp_t a,r;   load_fp(a,i[0]); fp_inv(r,a);     store_fp(o[0],r); }
static void o_fp_neg(mpz_t *i, mpz_t *o){ fp_t a,r;   load_fp(a,i[0]); fp_neg(r,a);     store_fp(o[0],r); }

/* --- fp2 --- */
#define FP2_BIN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp2_t a,b,r; LOADN(a,i,2); LOADN(b,(i+2),2); fn(r,a,b); STOREN(o,r,2); }
#define FP2_UN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp2_t a,r; LOADN(a,i,2); fn(r,a); STOREN(o,r,2); }
FP2_BIN(o_fp2_add, fp2_add) FP2_BIN(o_fp2_sub, fp2_sub) FP2_BIN(o_fp2_mul, fp2_mul)
FP2_UN(o_fp2_sqr, fp2_sqr)  FP2_UN(o_fp2_inv, fp2_inv)  FP2_UN(o_fp2_xi, fp2_mul_xi)

/* --- fp6 : six fp per element, laid out c0.a c0.b c1.a c1.b c2.a c2.b --- */
static void ld6(fp6_t r, mpz_t *v){ for(int k=0;k<3;k++){ load_fp(r[k][0], v[2*k]); load_fp(r[k][1], v[2*k+1]); } }
static void st6(mpz_t *o, const fp6_t a){ for(int k=0;k<3;k++){ store_fp(o[2*k], a[k][0]); store_fp(o[2*k+1], a[k][1]); } }
#define FP6_BIN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp6_t a,b,r; ld6(a,i); ld6(b,i+6); fn(r,a,b); st6(o,r); }
#define FP6_UN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp6_t a,r; ld6(a,i); fn(r,a); st6(o,r); }
FP6_BIN(o_fp6_add, fp6_add) FP6_BIN(o_fp6_sub, fp6_sub) FP6_BIN(o_fp6_mul, fp6_mul)
FP6_UN(o_fp6_sqr, fp6_sqr)  FP6_UN(o_fp6_inv, fp6_inv)  FP6_UN(o_fp6_v, fp6_mul_v)

/* --- fp12 : d0 then d1, six fp each --- */
static void ld12(fp12_t r, mpz_t *v){ ld6(r[0], v); ld6(r[1], v+6); }
static void st12(mpz_t *o, const fp12_t a){ st6(o, a[0]); st6(o+6, a[1]); }
#define FP12_BIN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp12_t a,b,r; ld12(a,i); ld12(b,i+12); fn(r,a,b); st12(o,r); }
#define FP12_UN(nm, fn) static void nm(mpz_t *i, mpz_t *o){ \
    fp12_t a,r; ld12(a,i); fn(r,a); st12(o,r); }
FP12_BIN(o_fp12_add, fp12_add) FP12_BIN(o_fp12_sub, fp12_sub) FP12_BIN(o_fp12_mul, fp12_mul)
FP12_UN(o_fp12_sqr, fp12_sqr)  FP12_UN(o_fp12_inv, fp12_inv)


/* --- elliptic curve records ------------------------------------------------
 * Vectors store affine points as (infinity flag, x, y). The new layer is
 * Jacobian, so each record converts in, operates, and converts back out. That
 * exercises the conversions as well as the group law. */
static void put_ep(mpz_t *o, const ep_t *P)
{
    fp_t x, y;
    if (!ep_to_affine(x, y, P)) { mpz_set_ui(o[0],1); mpz_set_ui(o[1],0); mpz_set_ui(o[2],0); return; }
    mpz_set_ui(o[0], 0); store_fp(o[1], x); store_fp(o[2], y);
}
static void get_ep(ep_t *P, mpz_t *i)
{
    if (mpz_cmp_ui(i[0], 0) != 0) { ep_set_infinity(P); return; }
    fp_t x, y; load_fp(x, i[1]); load_fp(y, i[2]);
    ep_from_affine(P, x, y);
}
static void put_ep2(mpz_t *o, const ep2_t *P)
{
    fp2_t x, y;
    if (!ep2_to_affine(x, y, P)) { for (int k=0;k<5;k++) mpz_set_ui(o[k],0); mpz_set_ui(o[0],1); return; }
    mpz_set_ui(o[0], 0);
    store_fp(o[1], x[0]); store_fp(o[2], x[1]);
    store_fp(o[3], y[0]); store_fp(o[4], y[1]);
}
static void get_ep2(ep2_t *P, mpz_t *i)
{
    if (mpz_cmp_ui(i[0], 0) != 0) { ep2_set_infinity(P); return; }
    fp2_t x, y;
    load_fp(x[0], i[1]); load_fp(x[1], i[2]);
    load_fp(y[0], i[3]); load_fp(y[1], i[4]);
    ep2_from_affine(P, x, y);
}
static void scalar_limbs(limb_t *k, int *bits, const mpz_t s)
{
    memset(k, 0, sizeof(limb_t) * FP_LIMBS);
    mpz_export(k, NULL, -1, sizeof(limb_t), 0, 0, s);
    *bits = (int)mpz_sizeinbase(s, 2);
}
static void o_efp_dbl(mpz_t *i, mpz_t *o){ ep_t P,R; get_ep(&P,i); ep_dbl(&R,&P); put_ep(o,&R); }
static void o_efp_add(mpz_t *i, mpz_t *o){ ep_t P,Q,R; get_ep(&P,i); get_ep(&Q,i+3); ep_add(&R,&P,&Q); put_ep(o,&R); }
static void o_efp_mul(mpz_t *i, mpz_t *o){
    ep_t P,R; get_ep(&P,i);
    limb_t k[FP_LIMBS]; int bits; scalar_limbs(k,&bits,i[3]);
    ep_mul(&R,&P,k,bits); put_ep(o,&R); }
static void o_efp2_dbl(mpz_t *i, mpz_t *o){ ep2_t P,R; get_ep2(&P,i); ep2_dbl(&R,&P); put_ep2(o,&R); }
static void o_efp2_add(mpz_t *i, mpz_t *o){ ep2_t P,Q,R; get_ep2(&P,i); get_ep2(&Q,i+5); ep2_add(&R,&P,&Q); put_ep2(o,&R); }
static void o_efp2_mul(mpz_t *i, mpz_t *o){
    ep2_t P,R; get_ep2(&P,i);
    limb_t k[FP_LIMBS]; int bits; scalar_limbs(k,&bits,i[5]);
    ep2_mul(&R,&P,k,bits); put_ep2(o,&R); }

static const Op OPS[] = {
    {"fp_add",2,1,o_fp_add}, {"fp_sub",2,1,o_fp_sub}, {"fp_mul",2,1,o_fp_mul},
    {"fp_inv",1,1,o_fp_inv}, {"fp_neg",1,1,o_fp_neg},
    {"fp2_add",4,2,o_fp2_add}, {"fp2_sub",4,2,o_fp2_sub}, {"fp2_mul",4,2,o_fp2_mul},
    {"fp2_sqr",2,2,o_fp2_sqr}, {"fp2_inv",2,2,o_fp2_inv}, {"fp2_mulbasis",2,2,o_fp2_xi},
    {"fp6_add",12,6,o_fp6_add}, {"fp6_sub",12,6,o_fp6_sub}, {"fp6_mul",12,6,o_fp6_mul},
    {"fp6_sqr",6,6,o_fp6_sqr}, {"fp6_inv",6,6,o_fp6_inv}, {"fp6_mulbasis",6,6,o_fp6_v},
    {"fp12_add",24,12,o_fp12_add}, {"fp12_sub",24,12,o_fp12_sub},
    {"fp12_mul",24,12,o_fp12_mul}, {"fp12_sqr",12,12,o_fp12_sqr},
    {"fp12_inv",12,12,o_fp12_inv},
    {"efp_dbl",3,3,o_efp_dbl}, {"efp_add",6,3,o_efp_add}, {"efp_mul",4,3,o_efp_mul},
    {"efp2_dbl",5,5,o_efp2_dbl}, {"efp2_add",10,5,o_efp2_add}, {"efp2_mul",6,5,o_efp2_mul},
    {NULL,0,0,NULL}
};

/* Any op outside the table is a hard error, so a typo in a vector file cannot
 * be mistaken for an unimplemented feature. */
static const char *KNOWN_SKIP[] = { NULL };

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: kat_runner_new <file.vec>\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }

    mpz_init(PRIME);
    mpz_t in[40], out[16], want[16], got;
    for (int i = 0; i < 40; i++) mpz_init(in[i]);
    for (int i = 0; i < 16; i++) { mpz_init(out[i]); mpz_init(want[i]); }
    mpz_init(got);

    char *line = NULL; size_t cap = 0; ssize_t len;
    while ((len = getline(&line, &cap, f)) > 0) {
        cur_line++;
        if (line[0] == '#' || line[0] == '\n') continue;
        char *save = NULL, *tok = strtok_r(line, " \t\n", &save);
        if (!tok) continue;

        if (strcmp(tok, "prime") == 0) {
            char *h = strtok_r(NULL, " \t\n", &save);
            mpz_set_str(PRIME, h, 16);
            /* The vector file's prime must be the one compiled into this
             * binary, or we would be silently testing the wrong curve. */
            mpz_t mine; mpz_init(mine);
            limb_t m[FP_LIMBS]; memcpy(m, FP_MODULUS, sizeof m);
            mpz_import(mine, FP_LIMBS, -1, sizeof(limb_t), 0, 0, m);
            if (mpz_cmp(mine, PRIME) != 0) {
                gmp_fprintf(stderr, "curve mismatch: built for %s (p=%Zx) but "
                                    "vectors are for p=%Zx\n",
                            ELIPS_CURVE_NAME, mine, PRIME);
                return 2;
            }
            mpz_clear(mine);
            continue;
        }

        const Op *op = NULL;
        for (int i = 0; OPS[i].name; i++)
            if (strcmp(tok, OPS[i].name) == 0) { op = &OPS[i]; break; }
        if (!op) {
            int known = 0;
            for (int i = 0; KNOWN_SKIP[i]; i++)
                if (strcmp(tok, KNOWN_SKIP[i]) == 0) { known = 1; break; }
            if (!known) { fprintf(stderr, "unknown op '%s' line %ld\n", tok, cur_line); return 2; }
            n_skip++; continue;
        }
        snprintf(cur_op, sizeof cur_op, "%s", op->name);

        for (int i = 0; i < op->in_n; i++) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h) { fprintf(stderr, "short record line %ld\n", cur_line); return 2; }
            mpz_set_str(in[i], h, 16);
        }
        char *eq = strtok_r(NULL, " \t\n", &save);
        if (!eq || strcmp(eq, "=")) { fprintf(stderr, "missing '=' line %ld\n", cur_line); return 2; }
        for (int i = 0; i < op->out_n; i++) {
            char *h = strtok_r(NULL, " \t\n", &save);
            if (!h) { fprintf(stderr, "short result line %ld\n", cur_line); return 2; }
            mpz_set_str(want[i], h, 16);
        }

        long before = n_fail;
        op->run(in, out);
        for (int i = 0; i < op->out_n; i++) expect(out[i], want[i], i);
        if (n_fail == before) n_pass++;
    }
    free(line); fclose(f);

    /* Release everything. LeakSanitizer is enabled by default alongside
     * AddressSanitizer on Linux and reports leaks in the harness as loudly as
     * leaks in the library, so a test that abandons its mpz_t turns the Asan CI
     * job red for no reason -- and a red job that is known to be noise stops
     * being read. */
    for (int i = 0; i < 40; i++) mpz_clear(in[i]);
    for (int i = 0; i < 16; i++) { mpz_clear(out[i]); mpz_clear(want[i]); }
    mpz_clear(got);
    mpz_clear(PRIME);

    printf("%s [%s]: %ld passed, %ld failed", argv[1], ELIPS_CURVE_NAME, n_pass, n_fail);
    if (n_skip) printf(", %ld skipped (elliptic-curve records)", n_skip);
    printf("\n");
    return n_fail ? 1 : 0;
}
