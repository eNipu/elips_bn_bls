"""Derive the GLV decomposition on G1 for both families, and prove the bounds.

A GLV routine that decomposes a scalar slightly wrong does not fail loudly. It
returns a point, and that point is correct for most scalars and wrong near a
digit boundary. So this file does not implement a recalled basis: it reduces
the lattice per curve, proves the digit bound that the ladder length depends
on, and then checks the whole decomposition against real curve points.

THE SETUP

E has an endomorphism phi(x, y) = (beta*x, y) with beta of order 3 in Fp.
phi^3 = 1 and phi != 1, so phi^2 + phi + 1 = 0, and on the order-r subgroup phi
acts as multiplication by a root lambda of that quadratic mod r. Writing

    k == k1 + k2*lambda   (mod r)

turns [k]P into [k1]P + [k2]phi(P), which a joint ladder evaluates in
max(bits(k1), bits(k2)) doublings instead of bits(r).

That is only a win if k1 and k2 are short, and shortness comes from the lattice

    L = { (z1, z2) in Z^2 : z1 + z2*lambda == 0 (mod r) },

which has determinant r and so has vectors of norm about sqrt(r). Given a
reduced basis v1, v2, Babai rounding of (k, 0) lands within half a basis vector
of the lattice, and the remainder is the pair (k1, k2).

WHAT COMES OUT, PER FAMILY

BLS12 needs no lattice reduction. lambda = -x^2 mod r, so (x^2, 1) is already
in L and already of norm about sqrt(r) -- x^4 is r to within the low-degree
terms. Reduction confirms it: the reduced basis is (x^2, 1) and (-1, x^2) up to
sign. So the decomposition is nothing more than writing k in base x^2, which is
one division and no rounding. That is what src/arith/ec.c does there.

BN is the case that needs this file. lambda = 36x^3 + 18x^2 + 6x + 1 is a
degree-3 polynomial in x, about 348 bits against sqrt(r) = 231, so no base-B
split exists and the basis has to be reduced. Reduction gives, in closed form,

    v1 = (-(2x + 1),      6x^2 + 2x)
    v2 = (-(6x^2+4x+1),  -(2x + 1))

Writing a = 2x + 1, c = 6x^2 + 2x and d = 6x^2 + 4x + 1 = a + c, the two
vectors are (-a, c) and (-d, -a), their determinant a^2 + c*d is exactly r, and
the decomposition of k is

    q1 = round(k*a / r)        q2 = round(k*c / r)
    k1 = k - q1*a - q2*d       k2 = q1*c - q2*a

with both q's non-negative for 0 <= k < r. Only three constants are needed and
d is derived from the other two, so src/arith/ec.c stores a and c alone.

THE BOUND

Since (k, 0) = b1*v1 + b2*v2 over the rationals and the q's round each b to
within 1/2, the remainder is -e1*v1 - e2*v2 with |e_i| <= 1/2. Componentwise

    |k1| <= (|v1x| + |v2x|) / 2      |k2| <= (|v1y| + |v2y|) / 2

which for BN is under 3x^2 + 3x + 1, or 231 bits: half the order, so the ladder
halves. That is a proof, not a sample, and the sampling below only confirms it.

Exit status 0 if every basis reduces, every bound is proved and the
decomposition reproduces [k]P on real points.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, EFp


def rdiv(a, b):
    """round(a/b) for b > 0, exact for integers of any size.

    Python's a/b would go through a float and lose every bit past the 53rd,
    which for a 700-bit numerator is all of them.
    """
    return (2 * a + b) // (2 * b)


def gauss_reduce(v1, v2):
    """Lagrange-Gauss reduction of a rank-2 lattice basis."""
    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1]

    if dot(v1, v1) > dot(v2, v2):
        v1, v2 = v2, v1
    while True:
        mu = rdiv(dot(v1, v2), dot(v1, v1))
        v2 = (v2[0] - mu * v1[0], v2[1] - mu * v1[1])
        if dot(v2, v2) >= dot(v1, v1):
            return v1, v2
        v1, v2 = v2, v1


_GEN_CACHE = {}


def g1_generator(cv):
    """Smallest x making x^3 + b a square, lifted and cleared to order r."""
    if cv.name in _GEN_CACHE:
        return _GEN_CACHE[cv.name]
    p, r = cv.p, cv.r
    tr = (cv.X + 1) if cv.family == "bls12" else (6 * cv.X ** 2 + 1)
    h1 = (p + 1 - tr) // r
    for xx in range(1, 4000):
        rhs = (xx ** 3 + cv.b_signed) % p
        if pow(rhs, (p - 1) // 2, p) != 1:
            continue
        cand = EFp(cv, xx, pow(rhs, (p + 1) // 4, p)).mul(h1)
        if not cand.inf and cand.mul(r).inf:
            _GEN_CACHE[cv.name] = cand
            return cand
    raise AssertionError("no G1 generator found")


_EIG_CACHE = {}


def eigenvalue(cv):
    """The (beta, lambda) pair: beta of order 3 in Fp with phi = [lambda] on G1.

    Both cube roots of unity in Fp act as SOME eigenvalue, and the two
    eigenvalues are lambda and -1-lambda. Which one the library wants is fixed
    by what the rest of the code assumes, so it is named per family here rather
    than picked by size:

        BLS12   -x^2 mod r, because src/arith/ec.c splits in base x^2 and
                negates phi(P). The OTHER root is x^2 - 1, which is shorter,
                and choosing it by size would silently break that code.
        BN      36x^3 + 18x^2 + 6x + 1, the root that reduces to the short
                basis below. The other is a full-length 462 bits.

    The named root is then matched to a beta on a real point, so a wrong name
    fails here instead of in the ladder.
    """
    key = cv.name
    if key in _EIG_CACHE:
        return _EIG_CACHE[key]
    p, r, x = cv.p, cv.r, cv.X
    want = (-x * x) % r if cv.family == "bls12" else \
           (36 * x ** 3 + 18 * x ** 2 + 6 * x + 1) % r
    assert pow(want, 3, r) == 1 and want != 1, "lambda must have order exactly 3"
    g = g1_generator(cv)
    for t in range(2, 200):
        b = pow(t, (p - 1) // 3, p)
        if b == 1:
            continue
        if EFp(cv, b * g.x % p, g.y) == g.mul(want):
            assert pow(b, 3, p) == 1, "beta must have order exactly 3"
            _EIG_CACHE[key] = (b, want)
            return b, want
    raise AssertionError("no cube root of unity in Fp acts as [lambda] on G1")


_BASIS_CACHE = {}


def basis(cv):
    """A reduced basis of L, as ((v1x, v1y), (v2x, v2y))."""
    if cv.name in _BASIS_CACHE:
        return _BASIS_CACHE[cv.name]
    _beta, lam = eigenvalue(cv)
    v1, v2 = gauss_reduce((cv.r, 0), (-lam, 1))
    for v in (v1, v2):
        assert (v[0] + v[1] * lam) % cv.r == 0, "basis vector is not in the lattice"
    assert abs(v1[0] * v2[1] - v1[1] * v2[0]) == cv.r, "determinant must be r"
    _BASIS_CACHE[cv.name] = (v1, v2)
    return v1, v2


_BN_CACHE = {}


def bn_constants(cv):
    """(a, c, d) for BN, with the closed forms checked against the reduction."""
    assert cv.family == "bn"
    if cv.name in _BN_CACHE:
        return _BN_CACHE[cv.name]
    x = cv.X
    a, c = 2 * x + 1, 6 * x * x + 2 * x
    d = a + c
    assert d == 6 * x * x + 4 * x + 1, "d must be the third basis entry"
    assert a * a + c * d == cv.r, "the basis determinant must be r"
    # The reduction is what proves these are short; the closed forms only make
    # them cheap to emit. Compare as sets of +-vectors, since Gauss reduction
    # fixes neither the order nor the sign.
    got = {tuple(sorted((abs(v[0]), abs(v[1])))) for v in basis(cv)}
    want = {tuple(sorted((a, c))), tuple(sorted((d, a)))}
    assert got == want, "closed form disagrees with the reduced basis"
    _BN_CACHE[cv.name] = (a, c, d)
    return a, c, d


def decompose(cv, k):
    """(k1, k2) with k == k1 + k2*lambda mod r, both about sqrt(r)."""
    r = cv.r
    if cv.family == "bn":
        a, c, d = bn_constants(cv)
        q1, q2 = rdiv(k * a, r), rdiv(k * c, r)
        return k - q1 * a - q2 * d, q1 * c - q2 * a
    # BLS12: lambda is -x^2, so writing k = e0 + e1*x^2 gives
    # k == e0 + (-e1)*lambda, and the second digit is the negated quotient.
    b = cv.X * cv.X
    return k % b, -(k // b)


def digit_bound(cv):
    """The proved bound on max(|k1|, |k2|), and its bit length."""
    if cv.family == "bn":
        v1, v2 = basis(cv)
        bx = (abs(v1[0]) + abs(v2[0])) // 2 + 1
        by = (abs(v1[1]) + abs(v2[1])) // 2 + 1
        return max(bx, by)
    return cv.X * cv.X


def main():
    random.seed(20260907)
    bad = 0
    for name, cv in sorted(CURVES.items()):
        beta, lam = eigenvalue(cv)
        v1, v2 = basis(cv)
        bound = digit_bound(cv)
        print("%s  family=%s" % (name, cv.family))
        print("  r        %4d bits" % cv.r.bit_length())
        print("  lambda   %4d bits" % lam.bit_length())
        print("  v1 = (%d, %d)   %d, %d bits"
              % (v1[0], v1[1], abs(v1[0]).bit_length(), abs(v1[1]).bit_length()))
        print("  v2 = (%d, %d)   %d, %d bits"
              % (v2[0], v2[1], abs(v2[0]).bit_length(), abs(v2[1]).bit_length()))
        print("  digit bound   %4d bits  (ladder is %.2fx the plain one)"
              % (bound.bit_length(), bound.bit_length() / cv.r.bit_length()))
        if cv.family == "bn":
            a, c, d = bn_constants(cv)
            print("  a = 2x+1        %d bits" % a.bit_length())
            print("  c = 6x^2+2x     %d bits" % c.bit_length())
            print("  d = a + c       %d bits" % d.bit_length())
            print("  widest product  k*c is %d bits" % (cv.r * c).bit_length())

        # The bound is proved above; this only confirms the proof, and would
        # catch a rounding convention that is off by one.
        g = g1_generator(cv)
        worst = 0
        for i in range(400):
            k = (0, 1, cv.r - 1, 2)[i] if i < 4 else random.randrange(cv.r)
            k1, k2 = decompose(cv, k)
            if (k1 + k2 * lam - k) % cv.r != 0:
                print("  FAIL: k = %d does not recombine" % k)
                bad += 1
                break
            worst = max(worst, abs(k1), abs(k2))
        if worst > bound:
            print("  FAIL: sampled digit %d bits exceeds the proved bound %d"
                  % (worst.bit_length(), bound.bit_length()))
            bad += 1
        print("  sampled max   %4d bits" % worst.bit_length())

        # And the decomposition has to move real points, not just integers.
        for _ in range(12):
            pt = g.mul(random.randrange(1, cv.r))
            k = random.randrange(cv.r)
            k1, k2 = decompose(cv, k)
            phi = EFp(cv, beta * pt.x % cv.p, pt.y)
            if phi != pt.mul(lam):
                print("  FAIL: phi is not [lambda] on G1")
                bad += 1
                break
            if pt.mul(k) != pt.mul(k1) + phi.mul(k2):
                print("  FAIL: [k]P != [k1]P + [k2]phi(P)")
                bad += 1
                break
        else:
            print("  recombination on curve points OK")
        print()
    print("glv_ref: %s" % ("FAILED" if bad else "OK"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
