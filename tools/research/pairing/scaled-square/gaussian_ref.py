"""Research-only denominator-free Gaussian interpolation for Miller squaring.

The kernel returns FOUR times a^2, never the exact raw square. It may replace
general squaring in a reduced-pairing-only path, not a bit-identical Miller API.
Variable-time Python; not production code or a publication-level novelty claim.
"""
import random
import sys
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "tools/reference"))
from elips_ref import CURVES, Fp2, Fp6, Fp12
import normalized_miller_ref as nm


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def twice(a):
    return add(a, a)


def times_u(a):
    return tuple(Fp2(c.p, -c.b % c.p, c.a) for c in a)


def times_s(a):
    return a[1].mul_xi(), a[0]


def sqr4(a):
    x, y = a
    xx, yy = x.sqr(), y.sqr()
    xy2 = (x + y).sqr() - xx - yy
    return xx + yy.mul_xi(), xy2


def gaussian_square(f):
    """Return 4*f^2 with 5 S4 = 15 S2 = 30 Mp, plus linear operations.

    Fp12 = Fp4[w]/(w^3-s), s^2=xi; evaluate at 0, 1, -1, u, infinity.
    u^2=-1 already in Fp2, so multiplication by u is only swap/negate.
    Clear the interpolation denominator 4 rather than performing division.
    """
    a = (f.d0.c0, f.d1.c1)
    b = (f.d1.c0, f.d0.c2)
    c = (f.d0.c1, f.d1.c2)
    p0, pinf = sqr4(a), sqr4(c)
    ac = add(a, c)
    pp, pm = sqr4(add(ac, b)), sqr4(sub(ac, b))
    pu = sqr4(add(sub(a, c), times_u(b)))
    S, D = add(pp, pm), sub(pp, pm)
    fourT = twice(twice(add(p0, pinf)))
    # J = 4*a*b - 4*b*c, as -u*(2*pu+S-4*(p0+pinf)).
    t = times_u(sub(add(twice(pu), S), fourT))
    J = (-t[0], -t[1])
    h0 = add(twice(twice(p0)), times_s(sub(D, J)))
    h1 = add(add(D, J), times_s(twice(twice(pinf))))
    h2 = sub(twice(S), fourT)
    return Fp12(Fp6(h0[0], h2[0], h1[1]),
                Fp6(h1[0], h0[1], h2[1]))


def fft4_square(f):
    """Return 4*f^2 from a four-point Fourier transform plus infinity.

    Evaluate at 1, -1, u, -u, infinity. This has the same 15 S2 cost as
    gaussian_square, but less interpolation work (19 rather than 23 A4).
    """
    a = (f.d0.c0, f.d1.c1)
    b = (f.d1.c0, f.d0.c2)
    c = (f.d0.c1, f.d1.c2)
    ac, amc, ub = add(a, c), sub(a, c), times_u(b)
    pp, pm = sqr4(add(ac, b)), sqr4(sub(ac, b))
    pu, mu = sqr4(add(amc, ub)), sqr4(sub(amc, ub))
    c4 = twice(twice(sqr4(c)))
    S, T, U = add(pp, pm), add(pu, mu), sub(pp, pm)
    t = times_u(sub(pu, mu))
    V = (-t[0], -t[1])
    h0 = add(sub(add(S, T), c4), times_s(sub(U, V)))
    h1 = add(add(U, V), times_s(c4))
    h2 = sub(S, T)
    return Fp12(Fp6(h0[0], h2[0], h1[1]),
                Fp6(h1[0], h0[1], h2[1]))


def scale12(f, k):
    def six(a):
        return Fp6(*(c.mul_int(k) for c in (a.c0, a.c1, a.c2)))
    return Fp12(six(f.d0), six(f.d1))


def replay(terms, square=None):
    f, pos, scale = Fp12.one(nm.P_MOD), 0, 1
    prepared = [(P, table, pow(P.y, -1, nm.P_MOD)) for P, table in terms]
    for digit in nm.DIGITS:
        f = square(f) if square else f.sqr()
        scale = scale * scale * (4 if square else 1) % nm.P_MOD
        for P, table, iy in prepared:
            xy = P.x * iy % nm.P_MOD
            for j in range(pos, pos + 1 + (digit != 0)):
                line = table[j]
                f = nm.mul_line_scaled8(f, line.beta.mul_int(iy),
                                       line.gamma.mul_int(xy))
        pos += 1 + (digit != 0)
    return f, scale


def main():
    rng = random.Random(2026090903)
    for cv in CURVES.values():
        def r2():
            return Fp2(cv.p, rng.randrange(cv.p), rng.randrange(cv.p))
        def r12():
            return Fp12(Fp6(r2(), r2(), r2()), Fp6(r2(), r2(), r2()))
        zero = Fp12.zero(cv.p)
        cases = [zero, Fp12.one(cv.p), -Fp12.one(cv.p)]
        # Basis vectors, equal and opposite coefficients, then full-width cases.
        for i in range(6):
            cs = [Fp2.zero(cv.p) for _ in range(6)]
            cs[i] = Fp2.one(cv.p)
            cases.append(Fp12(Fp6(*cs[:3]), Fp6(*cs[3:])))
        cases += [r12() for _ in range(1000)]
        for f in cases:
            check(gaussian_square(f) == scale12(f.sqr(), 4),
                  f"{cv.name}: scaled square mismatch")
            check(fft4_square(f) == scale12(f.sqr(), 4),
                  f"{cv.name}: Fourier square mismatch")
        counts = {"S2": 0, "M2": 0, "A4": 0}
        original_sqr, original_mul = Fp2.sqr, Fp2.__mul__
        original_add, original_sub = add, sub
        def counted_sqr(a):
            counts["S2"] += 1
            return original_sqr(a)
        def counted_mul(a, b):
            counts["M2"] += 1
            return original_mul(a, b)
        def counted_add(a, b):
            counts["A4"] += 1
            return original_add(a, b)
        def counted_sub(a, b):
            counts["A4"] += 1
            return original_sub(a, b)
        for square, additions in ((gaussian_square, 23), (fft4_square, 19)):
            counts = {"S2": 0, "M2": 0, "A4": 0}
            with patch.object(Fp2, "sqr", counted_sqr), \
                 patch.object(Fp2, "__mul__", counted_mul), \
                 patch.dict(globals(), {"add": counted_add, "sub": counted_sub}):
                square(cases[-1])
            check(counts == {"S2": 15, "M2": 0, "A4": additions}, str(counts))
        E = (cv.p**12 - 1) // cv.r
        check(E % (cv.p-1) == 0, "scale not killed")
        print(f"PASS {cv.name}: {len(cases)} general squares; {counts}", flush=True)

    records = nm._load_vectors(nm.DEFAULT_VECTORS)
    tables = []
    for i, (P, Q, expected) in enumerate(records):
        table = nm.prepare_g2(Q)
        tables.append((P, table))
        old, _ = replay([(P, table)])
        for square in (gaussian_square, fft4_square):
            new, scale = replay([(P, table)], square)
            check(new == scale12(old, scale), f"KAT {i}: scale recurrence")
            check(nm.final_exponentiate(new) == expected, f"KAT {i}: exact pairing")
        print(f"PASS KAT {i}: both exact pairings and raw scale recurrences", flush=True)
    for terms in ([], tables[:2], tables, tables + tables[:1]):
        old, _ = replay(terms)
        new, scale = replay(terms, fft4_square)
        check(new == scale12(old, scale), "multi-pairing scale recurrence")
        check(nm.final_exponentiate(new) == nm.final_exponentiate(old),
              "multi-pairing final value")
        print(f"PASS shared product n={len(terms)}", flush=True)
    P, table = tables[0]
    cancellation, _ = replay([(P, table), (-P, table)], fft4_square)
    check(nm.final_exponentiate(cancellation).is_one(), "pairing cancellation")
    print("PASS e(P,Q)*e(-P,Q)=1", flush=True)
    # Negative control: an incorrect interpolation term must be detected.
    wrong = gaussian_square(Fp12.one(nm.P_MOD)) + Fp12.one(nm.P_MOD)
    check(wrong != scale12(Fp12.one(nm.P_MOD), 4), "negative control")
    print("PASS all Gaussian-squaring reference checks", flush=True)


if __name__ == "__main__":
    main()
