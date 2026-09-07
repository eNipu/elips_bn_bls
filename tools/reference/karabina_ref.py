"""Derive Karabina's compressed squaring, and decide from measurement whether
it is worth using here.

WHY THIS FILE EXISTS EVEN THOUGH THE LIBRARY DOES NOT USE THE RESULT.

Plan section 10.5 scheduled Karabina compression as a 10 to 15 percent win on
the final exponentiation. Once Granger-Scott landed (the section wrongly assumed
it was already there), the remaining margin was measured rather than assumed,
and it does not pay on two of the three curves. Keeping the derivation makes
that a decision anyone can re-run rather than a claim they have to trust, and
re-deciding costs one command if the inputs change.

THE FORMULAS ARE FITTED, NOT RECALLED.

The compressed squaring is a quadratic map, so its coefficients can be solved
for. Sample random cyclotomic elements, square them exactly, and solve the
resulting linear system over Fp2 for the coefficients of every degree-2 monomial.
If the recovered coefficients come out as small integers and small integers
times xi, and the fitted map then reproduces squaring on samples it never saw,
the fit found the formula.

Which coordinates to keep was found the same way. Fitting was attempted for all
fifteen ways of choosing four of the six Fp2 coordinates; keeping (g1, g2, g4,
g5) -- which in the Fp4 view is exactly the two coordinates c1 = (g1,g4) and
c2 = (g2,g5), dropping c0 = (g0,g3) -- is the one that works.

Exit status 0 if every formula fits, has small coefficients, and reproduces
squaring and decompression on held-out samples.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, Fp2, Fp6, Fp12

FIT_N, HOLD_N = 45, 10


def rand_fp12(p):
    r = lambda: Fp2(p, random.randrange(p), random.randrange(p))
    return Fp12(Fp6(r(), r(), r()), Fp6(r(), r(), r()))


def cyclotomic(p, f):
    g = f.conj() * f.inv()
    return g.frob(2) * g


def coords(g):
    """g0 .. g5 by powers of w. d0 carries w^0, w^2, w^4; d1 carries w^1, w^3, w^5."""
    return [g.d0.c0, g.d1.c0, g.d0.c1, g.d1.c1, g.d0.c2, g.d1.c2]


def solve(A, b, p):
    """Gaussian elimination over Fp2, which is a field since x^2+1 is
    irreducible for every p = 3 mod 4 here. Returns None if inconsistent."""
    n, m = len(A), len(A[0])
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    piv, r = [], 0
    for c in range(m):
        s = next((i for i in range(r, n) if not M[i][c].is_zero()), None)
        if s is None:
            continue
        M[r], M[s] = M[s], M[r]
        inv = M[r][c].inv()
        M[r] = [x * inv for x in M[r]]
        for i in range(n):
            if i != r and not M[i][c].is_zero():
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    if any(not M[i][m].is_zero() for i in range(r, n)):
        return None
    sol = [Fp2.zero(p)] * m
    for i, c in enumerate(piv):
        sol[c] = M[i][m]
    return sol


def monomials(v, p):
    out = [Fp2.one(p)] + list(v)
    for i in range(len(v)):
        for j in range(i, len(v)):
            out.append(v[i] * v[j])
    return out


def names_for(labels):
    return (["1"] + list(labels) +
            [f"{labels[i]}*{labels[j]}"
             for i in range(len(labels)) for j in range(i, len(labels))])


def small(c, p):
    """A coefficient is 'small' if it is a little integer, or one times xi.
    Anything else means the fit latched onto noise rather than the formula."""
    xi = Fp2(p, 1, 1)
    for k in range(-16, 17):
        if k and (c - Fp2(p, k % p, 0)).is_zero():
            return str(k)
        if k and (c - xi.mul_int(k)).is_zero():
            return f"{k}xi"
    return None


def render(sol, labels, p):
    nm = names_for(labels)
    terms = []
    for i, c in enumerate(sol):
        if c.is_zero():
            continue
        s = small(c, p)
        if s is None:
            return None
        terms.append(s if nm[i] == "1" else f"{s}*{nm[i]}")
    return " + ".join(terms)


def check_curve(name):
    cv = CURVES[name]
    p = cv.p
    ok = True
    print(f"=== {name} ===")

    fit = [coords(cyclotomic(p, rand_fp12(p))) for _ in range(FIT_N)]
    fitsq = [coords(cyclotomic(p, rand_fp12(p))) for _ in range(0)]
    pairs = []
    while len(pairs) < FIT_N:
        g = cyclotomic(p, rand_fp12(p))
        pairs.append((coords(g), coords(g.sqr())))
    hold = []
    while len(hold) < HOLD_N:
        g = cyclotomic(p, rand_fp12(p))
        hold.append((coords(g), coords(g.sqr())))

    KEEP = (1, 2, 4, 5)
    LBL = [f"g{i}" for i in KEEP]

    # ---- compressed squaring
    sq = {}
    for out in KEEP:
        A = [monomials([c[i] for i in KEEP], p) for c, _ in pairs]
        b = [h[out] for _, h in pairs]
        sol = solve(A, b, p)
        txt = render(sol, LBL, p) if sol else None
        if txt is None:
            ok = False
            print(f"  h{out}: no fit with small coefficients")
        else:
            sq[out] = sol
            print(f"  h{out} = {txt}")

    bad = 0
    for c, h in hold:
        m = monomials([c[i] for i in KEEP], p)
        for out in KEEP:
            acc = Fp2.zero(p)
            for coef, mm in zip(sq[out], m):
                acc = acc + coef * mm
            if not (acc - h[out]).is_zero():
                bad += 1
    print(f"  squaring reproduced on {HOLD_N} held-out samples: "
          f"{'yes' if bad == 0 else f'NO ({bad} wrong)'}")
    ok = ok and bad == 0

    # ---- decompression: recover g3, then g0
    tgt = [c[1] * c[3] * Fp2(p, 4, 0) for c in fit]
    A = [monomials([c[i] for i in KEEP], p) for c in fit]
    sol3 = solve(A, tgt, p)
    txt = render(sol3, LBL, p) if sol3 else None
    if txt is None:
        ok = False
        print("  4*g1*g3: no fit with small coefficients")
    else:
        print(f"  4*g1*g3 = {txt}")

    LBL5 = ["g1", "g2", "g3", "g4", "g5"]
    A5 = [monomials([c[1], c[2], c[3], c[4], c[5]], p) for c in fit]
    sol0 = solve(A5, [c[0] for c in fit], p)
    txt = render(sol0, LBL5, p) if sol0 else None
    if txt is None:
        ok = False
        print("  g0: no fit with small coefficients")
    else:
        print(f"  g0 = {txt}")

    # decompression must reproduce on held-out samples too
    bad = 0
    for c, _ in hold:
        m = monomials([c[i] for i in KEEP], p)
        n3 = Fp2.zero(p)
        for coef, mm in zip(sol3, m):
            n3 = n3 + coef * mm
        g3 = n3 * (c[1] * Fp2(p, 4, 0)).inv()
        if not (g3 - c[3]).is_zero():
            bad += 1
            continue
        m5 = monomials([c[1], c[2], g3, c[4], c[5]], p)
        g0 = Fp2.zero(p)
        for coef, mm in zip(sol0, m5):
            g0 = g0 + coef * mm
        if not (g0 - c[0]).is_zero():
            bad += 1
    print(f"  decompression reproduced on {HOLD_N} held-out samples: "
          f"{'yes' if bad == 0 else f'NO ({bad} wrong)'}")
    print("  NOTE: decompression divides by g1. g1 = 0 is an exceptional case,")
    print("        and it is a branch on a secret-derived value.")
    print()
    return ok and bad == 0


def main():
    random.seed(20260907)
    bad = 0
    for name in ("BLS12-381", "BLS12-461", "BN-462"):
        if not check_curve(name):
            bad += 1
    if bad:
        print(f"FAILED: {bad} curve(s) did not verify")
        return 1
    print("Karabina compression derived and verified on all curves.")
    print()
    print("Whether it is WORTH using is a separate, measured question; see")
    print("MODERNIZATION_PLAN.md section 10.5. On the parameters this library")
    print("ships, the break-even run is 40 to 48 squarings and the longest")
    print("uninterrupted run is 32 (BLS12-381), 33 (BLS12-461) and 87 (BN-462),")
    print("so it pays on one run of one curve and is a loss everywhere else.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
