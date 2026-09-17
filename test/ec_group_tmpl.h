/*
 * The body of ec_group_test, parameterised over the group.
 *
 * Instantiated twice, for E(Fp) and E'(Fp2), because since issue #51 both run
 * the same Jacobian formulas and both are reached with points an attacker
 * chose. Expects EC_PT, EC_F, EC_FT, EC_RANDF and EC_EXPECT_OFFSUB.
 */

#define CAT_(a,b) a##b
#define CAT(a,b)  CAT_(a,b)
#define PT(name)  CAT(EC_PT, CAT(_, name))
#define PTT       CAT(EC_PT, _t)
#define F(name)   CAT(EC_F,  CAT(_, name))
#define AFF       CAT(aff_, EC_PT)
#define L(name)   CAT(EC_PT, CAT(_grp_, name))

/* ---- the affine oracle. inf is carried as a flag, not as a coordinate ---- */
typedef struct { EC_FT x, y; int inf; } AFF;

static void L(aff_dbl)(AFF *r, const AFF *p)
{
    if (p->inf || F(is_zero)(p->y)) { r->inf = 1; return; }
    EC_FT l, t, x3;
    F(sqr)(l, p->x);
    F(add)(t, l, l); F(add)(l, t, l);            /* 3x^2 */
    F(add)(t, p->y, p->y);                       /* 2y   */
    F(inv)(t, t);
    F(mul)(l, l, t);                             /* lambda */
    F(sqr)(x3, l);
    F(sub)(x3, x3, p->x); F(sub)(x3, x3, p->x);
    F(sub)(t, p->x, x3);
    F(mul)(t, t, l);
    F(sub)(r->y, t, p->y);
    F(copy)(r->x, x3);
    r->inf = 0;
}

static void L(aff_add)(AFF *r, const AFF *p, const AFF *q)
{
    if (p->inf) { *r = *q; return; }
    if (q->inf) { *r = *p; return; }
    if (F(eq)(p->x, q->x)) {
        if (F(eq)(p->y, q->y)) { L(aff_dbl)(r, p); return; }
        r->inf = 1; return;                      /* q = -p */
    }
    EC_FT l, t, x3;
    F(sub)(l, q->y, p->y);
    F(sub)(t, q->x, p->x);
    F(inv)(t, t);
    F(mul)(l, l, t);
    F(sqr)(x3, l);
    F(sub)(x3, x3, p->x); F(sub)(x3, x3, q->x);
    F(sub)(t, p->x, x3);
    F(mul)(t, t, l);
    F(sub)(r->y, t, p->y);
    F(copy)(r->x, x3);
    r->inf = 0;
}

static void L(to_aff)(AFF *a, const PTT *P)
{
    a->inf = !PT(to_affine)(a->x, a->y, P);
}

static int L(aff_eq)(const AFF *a, const AFF *b)
{
    if (a->inf || b->inf) return a->inf && b->inf;
    return F(eq)(a->x, b->x) && F(eq)(a->y, b->y);
}

/* Almost never in the subgroup where the cofactor is large, which is the
 * point: the group law has to be right off the subgroup too. */
static int L(rand_point)(PTT *P)
{
    EC_FT x, y, t, b;
    EC_RANDF(x);
    F(sqr)(t, x); F(mul)(t, t, x);
    PT(curve_b)(b); F(add)(t, t, b);
    if (!F(sqrt)(y, t)) return 0;
    PT(from_affine)(P, x, y);
    return 1;
}

/* Same point, a different representative: (X:Y:Z) -> (l^2 X : l^3 Y : l Z). */
static void L(rescale)(PTT *P)
{
    EC_FT l, l2, l3;
    do { EC_RANDF(l); } while (F(is_zero)(l));
    F(sqr)(l2, l); F(mul)(l3, l2, l);
    F(mul)(P->x, P->x, l2);
    F(mul)(P->y, P->y, l3);
    F(mul)(P->z, P->z, l);
}

static void L(check_dbl)(const PTT *P, const char *what)
{
    PTT d; AFF a, b, pa;
    PT(dbl)(&d, P);
    L(to_aff)(&a, &d);
    L(to_aff)(&pa, P);
    L(aff_dbl)(&b, &pa);
    ok(L(aff_eq)(&a, &b), what);
    ok(PT(on_curve)(&d), "dbl output is on the curve");
    /* aliasing */
    PTT c; PT(copy)(&c, P); PT(dbl)(&c, &c);
    ok(PT(eq)(&c, &d), "aliased dbl agrees");
}

static void L(check_add)(const PTT *P, const PTT *Q, const char *what)
{
    PTT s; AFF a, b, pa, qa;
    PT(add)(&s, P, Q);
    L(to_aff)(&a, &s);
    L(to_aff)(&pa, P); L(to_aff)(&qa, Q);
    L(aff_add)(&b, &pa, &qa);
    ok(L(aff_eq)(&a, &b), what);
    ok(PT(on_curve)(&s), "add output is on the curve");
    /* commutative, and both aliasing directions */
    PTT t; PT(add)(&t, Q, P);
    ok(PT(eq)(&t, &s), "add is commutative");
    PT(copy)(&t, P); PT(add)(&t, &t, Q);
    ok(PT(eq)(&t, &s), "add aliasing the first argument agrees");
    PT(copy)(&t, Q); PT(add)(&t, P, &t);
    ok(PT(eq)(&t, &s), "add aliasing the second argument agrees");
}

/* Returns 0 on success, 2 if the run proved nothing, 1 if a check failed. */
static int L(run)(long want, const char *group)
{
    long used = 0, offsub = 0;
    int before = fails;
    PTT O; PT(set_infinity)(&O);

    for (long tried = 0; used < want && tried < want * 8; tried++) {
        PTT P, Q, N;
        if (!L(rand_point)(&P)) continue;
        used++;
        if (!PT(in_subgroup)(&P)) offsub++;
        L(rescale)(&P);

        L(check_dbl)(&P, "dbl matches affine");

        if (L(rand_point)(&Q)) {
            L(rescale)(&Q);
            L(check_add)(&P, &Q, "add matches affine");
        }

        /* the degenerate pairs, on every point rather than once */
        L(check_add)(&P, &P, "P + P matches affine (the doubling path)");
        PT(neg)(&N, &P);
        L(check_add)(&P, &N, "P + (-P) is the identity");
        L(check_add)(&P, &O, "P + O is P");
        L(check_add)(&O, &P, "O + P is P");
    }

    L(check_dbl)(&O, "dbl(O) is O");
    L(check_add)(&O, &O, "O + O is O");
    {
        PTT d; PT(dbl)(&d, &O);
        ok(PT(is_infinity)(&d), "dbl(O) reports infinity");
        PT(add)(&d, &O, &O);
        ok(PT(is_infinity)(&d), "O + O reports infinity");
    }
    /* O reached the way a ladder reaches it, from a real point */
    {
        PTT G, acc, d;
        PT(generator)(&G);
        PT(mul)(&acc, &G, ELIPS_ORDER, ELIPS_ORDER_BITS);
        ok(PT(is_infinity)(&acc), "[r]G is O");
        PT(dbl)(&d, &acc);
        ok(PT(is_infinity)(&d), "dbl of a computed O is O");
        PT(add)(&d, &acc, &G);
        ok(PT(eq)(&d, &G), "a computed O plus G is G");
    }

    printf("  %s: %ld points (%ld outside %s), %d checks, %d failed\n",
           group, used, offsub, group, checks, fails);
    if (used == 0) { printf("  %s: no points generated\n", group); return 2; }
    /* A run where every point landed in the subgroup proves nothing about the
     * off-subgroup behaviour it exists to test -- unless the cofactor is 1,
     * where there is no off-subgroup point to find. */
    if (offsub == 0 && EC_EXPECT_OFFSUB) {
        printf("  %s: every point was in the subgroup, test is vacuous\n", group);
        return 2;
    }
    if (offsub == 0)
        printf("  %s: cofactor is 1, so every curve point is in the subgroup\n", group);
    return fails > before ? 1 : 0;
}

#undef CAT_
#undef CAT
#undef PT
#undef PTT
#undef F
#undef AFF
#undef L
