/*
 * map_to_curve, written once and instantiated for E(Fp) and E'(Fp2).
 *
 * Same reasoning as src/arith/ec_tmpl.h: the two versions are the same algebra
 * over different fields, and two hand-copied versions drift. Here it matters
 * more than usual, because a hash-to-curve map that is subtly wrong still
 * returns points on the curve and in the right subgroup -- it just computes a
 * different function from every other implementation, and nothing but a test
 * vector will tell you.
 *
 * Expects H2C_PT (ep / ep2), H2C_F (fp / fp2), H2C_FT (fp_t / fp2_t) and
 * H2C_CN (the H2C_G1_ / H2C_G2_ constant prefix) from the includer.
 *
 * Constant time throughout. Both branches of every value-dependent choice are
 * computed and selected under a mask, including both square roots, because the
 * message being hashed is not always public: an oblivious PRF or a
 * password-authenticated exchange hashes a secret.
 */

#define H2C_CAT_(a, b) a##b
#define H2C_CAT(a, b)  H2C_CAT_(a, b)
#define M(name)        H2C_CAT(H2C_CN, name)
#define F(name)        H2C_CAT(H2C_F, H2C_CAT(_, name))
#define PT(name)       H2C_CAT(H2C_PT, H2C_CAT(_, name))
#define PTT            H2C_CAT(H2C_PT, _t)
#define FN(name)       H2C_CAT(H2C_CAT(H2C_PT, _h2c_), name)

/* r = mask ? a : b */
static void FN(cmov)(H2C_FT r, const H2C_FT a, const H2C_FT b, limb_t mask)
{
    F(cselect)(r, a, b, mask);
}

static limb_t FN(mask)(int cond) { return (limb_t)0 - (limb_t)(cond & 1); }

/* g(x) = x^3 + A x + B for the curve the map runs on. With A = 0 the middle
 * term drops out, which is why the SvdW branch below never mentions it. */
static void FN(gx)(H2C_FT out, const H2C_FT x, const H2C_FT A, const H2C_FT B)
{
    H2C_FT t;
    F(sqr)(t, x);
    F(mul)(t, t, x);
    F(mul)(out, A, x);
    F(add)(t, t, out);
    F(add)(out, t, B);
}

#if defined(ELIPS_H2C_SSWU)

/* RFC 9380 6.6.2, as its definition rather than as the straight-line listing.
 *
 *   x1 = (-B/A) (1 + 1/(Z^2 u^4 + Z u^2)),  or B/(ZA) when that denominator is 0
 *   x2 = Z u^2 x1
 *   x  = x1 if g(x1) is square, else x2
 *   y  = sqrt(g(x)), with sgn0(y) forced to sgn0(u)
 *
 * The two square roots are both computed and selected. fp_sqrt returns zero and
 * a false flag on a non-residue, so the unused branch is harmless.
 */
#if defined(H2C_SSWU_ONE_ROOT)
/* One exponentiation instead of four, for the map's square root.
 *
 * map_sswu needs sqrt(g(x1)) when g(x1) is square and sqrt(g(x2)) otherwise,
 * and it cannot look first: choosing would branch on the message. The obvious
 * answer, computing both roots and selecting, is what the other branch below
 * does, and it costs two square roots of two exponentiations each.
 *
 * ONE exponentiation gives either. q = p^2 is 9 mod 16 on every curve here, so
 * the 2-Sylow of Fp2* has order 8. Write y = g(x1)^((q+7)/16). Then
 *
 *     t = y^2 / g(x1)
 *
 * is an 8th root of unity, and its index k decides everything:
 *
 *     k even   g(x1) is the square, and (y * zeta^(-k/2))^2      == g(x1)
 *     k odd    g(x2) is the square, and (y * FIX[k] * u^3)^2     == g(x2)
 *
 * The odd case works because g(x2) = (Z u^2)^3 g(x1) holds by construction, so
 * the missing factor is sqrt(Z^3 zeta^-k). That exists precisely because both
 * Z^3 and an odd power of zeta are non-squares, and a product of two
 * non-squares is a square. gen_h2c_params.py derives FIX and re-checks the
 * whole table against the map before emitting it.
 *
 * t is never formed, because that would need an inversion. The loop compares
 * y^2 against g(x1)*zeta^k instead, which is eight multiplications against a
 * divstep, and selects under a mask so k stays hidden. k is compared with a
 * loop counter, which is public; only the match is data.
 *
 * Returns all-ones if x1 was the right abscissa, zero if x2. */
static limb_t FN(sswu_root)(H2C_FT y, const H2C_FT gx1, const H2C_FT u)
{
    H2C_FT cand, ysq, t, fix, u3, one, sel;
    limb_t even = 0;

    F(exp)(cand, gx1, M(SQRT_EXP), M(SQRT_EXP_BITS));
    F(sqr)(ysq, cand);

    F(set_zero)(fix);
    for (int k = 0; k < 8; k++) {
        F(mul)(t, gx1, H2C_ELEM(M(ZETA), k));
        limb_t m = FN(mask)(F(eq)(t, ysq));
        FN(cmov)(fix, H2C_ELEM(M(FIX), k), fix, m);
        if (!(k & 1)) even |= m;
    }

    F(sqr)(u3, u);
    F(mul)(u3, u3, u);
    F(set_one)(one);
    FN(cmov)(sel, one, u3, even);           /* even -> 1, odd -> u^3 */

    F(mul)(y, cand, fix);
    F(mul)(y, y, sel);
    return even;
}
#endif /* H2C_SSWU_ONE_ROOT */

static void FN(map_sswu)(H2C_FT xo, H2C_FT yo, const H2C_FT u)
{
    H2C_FT zu2, t, tinv, x1, x1e, x2, gx1, gx2, y1, y2, ny;

    F(sqr)(zu2, u);
    F(mul)(zu2, zu2, M(Z));

    F(sqr)(t, zu2);
    F(add)(t, t, zu2);                       /* Z^2 u^4 + Z u^2 */

    F(inv)(tinv, t);                         /* inv0: zero maps to zero */
    F(set_one)(x1);
    F(add)(x1, x1, tinv);
    F(mul)(x1, x1, M(MB_OVER_A));
    F(copy)(x1e, M(B_OVER_ZA));
    FN(cmov)(x1, x1e, x1, FN(mask)(F(is_zero)(t)));

    F(mul)(x2, zu2, x1);

    FN(gx)(gx1, x1, M(ISO_A), M(ISO_B));

#if defined(H2C_SSWU_ONE_ROOT)
    limb_t m = FN(sswu_root)(yo, gx1, u);
    FN(cmov)(xo, x1, x2, m);
    (void)gx2; (void)y1; (void)y2;
#else
    FN(gx)(gx2, x2, M(ISO_A), M(ISO_B));

    int is_sq = F(sqrt)(y1, gx1);
    (void)F(sqrt)(y2, gx2);

    limb_t m = FN(mask)(is_sq);
    FN(cmov)(xo, x1, x2, m);
    FN(cmov)(yo, y1, y2, m);
#endif

    F(neg)(ny, yo);
    FN(cmov)(yo, yo, ny, FN(mask)(F(sgn0)(yo) == F(sgn0)(u)));
}

/* x' = x_num(x)/x_den(x), y' = y * y_num(x)/y_den(x), by Horner.
 *
 * The coefficient tables are limb_t[n][FP_LIMBS] for Fp and
 * limb_t[n][2][FP_LIMBS] for Fp2, so indexing one generically needs the
 * includer's H2C_ELEM to rebuild the right pointer type. Casting to the field
 * type directly is not possible: these are array types, not pointer types. */
static void FN(horner)(H2C_FT out, const void *tbl, int n, const H2C_FT x)
{
    H2C_FT acc;
    F(copy)(acc, H2C_ELEM(tbl, n - 1));
    for (int i = n - 2; i >= 0; i--) {
        F(mul)(acc, acc, x);
        F(add)(acc, acc, H2C_ELEM(tbl, i));
    }
    F(copy)(out, acc);
}

#define HORNER(dst, tbl, x) \
    FN(horner)(dst, tbl, (int)(sizeof tbl / sizeof tbl[0]), x)

static void FN(iso_map)(H2C_FT xo, H2C_FT yo, const H2C_FT x, const H2C_FT y)
{
    H2C_FT xn, xd, yn, yd, t;

    HORNER(xn, M(ISO_XNUM), x);
    HORNER(xd, M(ISO_XDEN), x);
    HORNER(yn, M(ISO_YNUM), x);
    HORNER(yd, M(ISO_YDEN), x);

    F(inv)(t, xd);
    F(mul)(xo, xn, t);
    F(inv)(t, yd);
    F(mul)(yn, yn, y);
    F(mul)(yo, yn, t);
}

#undef HORNER

static void FN(map_to_curve)(PTT *out, const H2C_FT u)
{
    H2C_FT x, y, xi, yi;
    FN(map_sswu)(x, y, u);
    FN(iso_map)(xi, yi, x, y);
    PT(from_affine)(out, xi, yi);
}

#else  /* Shallue-van de Woestijne */

/* RFC 9380 6.6.1. Three candidate abscissae, of which the construction
 * guarantees at least one has g(x) square; take the first that does.
 *
 * All three square roots are computed and the answer selected, so the timing
 * says nothing about which candidate won. */
static void FN(map_to_curve)(PTT *out, const H2C_FT u)
{
    H2C_FT one, tv1, tv2, tv3, tv4, x1, x2, x3, gx1, gx2, gx3;
    H2C_FT y1, y2, y3, x, y, ny;

    F(set_one)(one);

    F(sqr)(tv1, u);
    F(mul)(tv1, tv1, M(C1));                 /* u^2 g(Z) */
    F(add)(tv2, one, tv1);
    F(sub)(tv1, one, tv1);

    F(mul)(tv3, tv1, tv2);
    F(inv)(tv3, tv3);                        /* inv0 */

    F(mul)(tv4, u, tv1);
    F(mul)(tv4, tv4, tv3);
    F(mul)(tv4, tv4, M(C3));

    F(sub)(x1, M(C2), tv4);
    F(add)(x2, M(C2), tv4);

    F(sqr)(x3, tv2);
    F(mul)(x3, x3, tv3);
    F(sqr)(x3, x3);
    F(mul)(x3, x3, M(C4));
    F(add)(x3, x3, M(Z));

    H2C_FT A, B;
    F(set_zero)(A);                          /* both curves have A = 0 */
    PT(curve_b)(B);

    FN(gx)(gx1, x1, A, B);
    FN(gx)(gx2, x2, A, B);
    FN(gx)(gx3, x3, A, B);

    int e1 = F(sqrt)(y1, gx1);
    int e2 = F(sqrt)(y2, gx2);
    (void)F(sqrt)(y3, gx3);

    limb_t m1 = FN(mask)(e1);
    limb_t m2 = FN(mask)(e2 & !e1);

    F(copy)(x, x3); F(copy)(y, y3);
    FN(cmov)(x, x2, x, m2); FN(cmov)(y, y2, y, m2);
    FN(cmov)(x, x1, x, m1); FN(cmov)(y, y1, y, m1);

    F(neg)(ny, y);
    FN(cmov)(y, y, ny, FN(mask)(F(sgn0)(y) == F(sgn0)(u)));

    PT(from_affine)(out, x, y);
}

#endif

#undef H2C_CAT_
#undef H2C_CAT
#undef M
#undef F
#undef PT
#undef PTT
#undef FN
