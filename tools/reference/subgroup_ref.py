"""Derive the fast subgroup membership tests, and prove they are exact.

A subgroup test that is too permissive does not fail loudly. It accepts
attacker-chosen points of small order, which is the hole the test exists to
close. So this file does not implement a recalled formula: it derives the test
per curve, proves the tests reject everything outside the subgroup, and then
checks the proof numerically against points that really are outside.

THE ARGUMENT

Both tests have the shape "E(P) == [m]P" for an endomorphism E and a short
multiplier m. The set of points satisfying it is ker(E - m), a subgroup. The
easy direction -- that G contains no counterexample -- is just the statement
that E acts as [m] on G. The hard direction is that nothing else satisfies it.

Every endomorphism used here satisfies a known quadratic over Z:

  G1, BLS12   phi(x, y) = (beta*x, y) with beta of order 3 in Fp.
              phi^3 = 1 and phi != 1, so phi^2 + phi + 1 = 0.

  G2, both    psi, the untwist-Frobenius-twist. Its characteristic polynomial
              is the Frobenius one, psi^2 - [t] psi + [p] = 0.

Suppose T has prime order l and satisfies the test. Then m is an eigenvalue of
the endomorphism acting on E[l], so m is a root of that quadratic mod l:

  G1   m^2 + m + 1  == 0  (mod l)
  G2   m^2 - t*m + p == 0  (mod l)

So if no prime dividing the cofactor divides the corresponding integer, no
point of prime order in the cofactor part can satisfy the test, the kernel
meets the cofactor part trivially, and the test is exact. Detecting a shared
prime needs no factorisation: it is a gcd.

WHAT COMES OUT

G1 on BLS12 is exact for free, with no condition on the seed. With m = -x^2,

    m^2 + m + 1 = x^4 - x^2 + 1 = r

exactly, and r is prime and does not divide the cofactor. That is why the
literature reports no exceptional BLS12 seeds; here it is a one-line identity
rather than a citation.

G1 on BN needs no test at all. #E(Fp) = r, so the cofactor is 1 and being on
the curve already puts a point in G1.

G2 is where a curve could fail, and where the published proof had to be
corrected once (ePrint 2022/352 fixing 2021/1130). The multiplier is p mod r,
which is x on BLS12 and 6x^2 on BN, and m^2 - t*m + p is not forced to be
coprime to the cofactor by any identity. It has to be checked per curve, and
this file checks it.

Exit status 0 if every test is derived, proved exact and confirmed numerically.
"""
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, Fp2, EFp, EFp2
from selftest import fp2_sqrt
from h2c_ref import group_orders as twist_group_orders


# ------------------------------------------------------------ group orders --

def trace(cv):
    """Trace of Frobenius, from the family's own parametrisation."""
    return (cv.X + 1) if cv.family == "bls12" else (6 * cv.X * cv.X + 1)


def group_orders(cv):
    """(#E(Fp), h1, #E'(Fp2), h2).

    The twist order comes from h2c_ref.group_orders, which does not pick
    between the two candidate traces algebraically: it builds a real point on
    the twist and keeps the order that annihilates it. Reused rather than
    reimplemented, so there is one place where this can be wrong.
    """
    r = cv.r
    nE, nE2 = twist_group_orders(cv)
    assert nE % r == 0,  "r must divide #E(Fp)"
    assert nE2 % r == 0, "r must divide #E'(Fp2)"
    # G is the unique subgroup of order r only if r appears to the first power.
    assert nE % (r * r) != 0,  "r^2 must not divide #E(Fp)"
    assert nE2 % (r * r) != 0, "r^2 must not divide #E'(Fp2)"
    return nE, nE // r, nE2, nE2 // r


# ---------------------------------------------------------- the G1 test ----

def g1_beta_and_lambda(cv):
    """The cube root of unity that makes phi act as [-x^2], and that multiplier.

    Both primitive cube roots of unity give a valid endomorphism; they act as
    the two roots of m^2+m+1 mod r, which are x^2-1 and -x^2. The pairing
    between them is found by testing, not assumed.
    """
    p, r, X = cv.p, cv.r, cv.X
    m = -X * X
    assert (m * m + m + 1) % r == 0, "-x^2 must be a primitive cube root mod r"
    assert m * m + m + 1 == r, "the identity the exactness proof rests on"

    _, h1, _, _ = group_orders(cv)
    P = None
    while P is None or P.inf:
        P = random_point_fp(cv).mul(h1)
    assert P.mul(r).inf, "P must have order r"

    for beta in cube_roots_of_unity(p):
        if EFp(cv, beta * P.x % p, P.y) == P.mul(m % r):
            return beta, m
    raise AssertionError("no cube root of unity makes phi act as [-x^2]")


def cube_roots_of_unity(p):
    assert (p - 1) % 3 == 0
    out = []
    while len(out) < 2:
        g = random.randrange(2, p)
        b = pow(g, (p - 1) // 3, p)
        if b != 1 and b not in out:
            out.append(b)
    return out


# ---------------------------------------------------------- the G2 test ----

def g2_multiplier(cv):
    """psi acts on G2 as [p mod r]. Return the short representative."""
    p, r, X = cv.p, cv.r, cv.X
    m = X if cv.family == "bls12" else 6 * X * X
    assert (p - m) % r == 0, "psi's multiplier must be p mod r"
    return m


def psi_constants(cv):
    """The same u2, u3 that gen_params.py bakes into PSI_X and PSI_Y."""
    p = cv.p
    xi = Fp2(p, 1, 1)
    g = fp2_pow(xi, (p - 1) // 6)
    u2 = (g * g).inv()
    u3 = (g * g * g).inv()
    return u2, u3


def fp2_pow(a, e):
    r = Fp2.one(a.p)
    while e:
        if e & 1:
            r = r * a
        a = a.sqr()
        e >>= 1
    return r


def psi(cv, Q, u2, u3):
    if Q.inf:
        return EFp2(cv)
    return EFp2(cv, Q.x.conj() * u2, Q.y.conj() * u3)


# ------------------------------------------------------- random points -----

def random_point_fp(cv):
    p = cv.p
    assert p % 4 == 3
    while True:
        x = random.randrange(p)
        y2 = (x * x * x + cv.b_signed) % p
        y = pow(y2, (p + 1) // 4, p)
        if (y * y - y2) % p == 0:
            return EFp(cv, x, y)


def random_point_fp2(cv):
    """A random point on the twist.

    Uses the project's own fp2_sqrt. An earlier version of this function raised
    to (p^2+1)/4, which is simply wrong: p = 3 mod 4 makes p^2 = 1 mod 8, so
    that exponent is not a square-root exponent over Fp2 and the sampler looped
    forever. The library already had the right routine.
    """
    p = cv.p
    b = EFp2.b_twist(cv)
    while True:
        x = Fp2(p, random.randrange(p), random.randrange(p))
        y = fp2_sqrt(x.sqr() * x + b, p)
        if y is not None:
            return EFp2(cv, x, y)


# ------------------------------------------------------------------ main ---

def check_curve(name):
    cv = CURVES[name]
    p, r, X = cv.p, cv.r, cv.X
    t = trace(cv)
    nE, h1, nE2, h2 = group_orders(cv)
    ok = True

    print(f"=== {name} ({cv.family}) ===")
    print(f"  h1 = {h1} ({h1.bit_length()} bits)")
    print(f"  h2 = {h2.bit_length()} bits")

    # ---- G1
    if cv.family == "bn":
        assert h1 == 1, "BN G1 cofactor must be 1"
        print("  G1: cofactor is 1, so on-curve implies in G1. No test needed.")
    else:
        beta, m1 = g1_beta_and_lambda(cv)
        witness = m1 * m1 + m1 + 1
        g = math.gcd(h1, witness)
        print(f"  G1: phi(P) == [-x^2]P     with beta pinned by construction")
        print(f"      m^2+m+1 == r : {witness == r}   gcd(h1, m^2+m+1) = {g}")
        if g != 1:
            ok = False
            print("      NOT EXACT: a cofactor prime is a cube root of unity")
        else:
            print("      exact: no cofactor prime can satisfy the relation")

        # numerical: G1 points pass, cofactor points fail
        passed_in = failed_out = 0
        for _ in range(12):
            R = random_point_fp(cv)
            P = R.mul(h1)                       # in G1
            if not P.inf:
                lhs = EFp(cv, beta * P.x % p, P.y)
                if lhs == P.mul(m1 % r):
                    passed_in += 1
                else:
                    ok = False
            T = R.mul(r)                        # in the cofactor part
            if not T.inf:
                lhs = EFp(cv, beta * T.x % p, T.y)
                if lhs != T.mul(m1 % r):
                    failed_out += 1
                else:
                    ok = False
                    print("      REJECTED POINT ACCEPTED: a non-G1 point passed")
        print(f"      numerically: {passed_in} G1 points accepted,"
              f" {failed_out} cofactor points rejected")

    # ---- G2
    m2 = g2_multiplier(cv)
    witness2 = m2 * m2 - t * m2 + p
    g2 = math.gcd(h2, witness2)
    lbl = "x" if cv.family == "bls12" else "6x^2"
    print(f"  G2: psi(Q) == [{lbl}]Q     ({m2.bit_length()} bits vs r's {r.bit_length()})")
    print(f"      gcd(h2, m^2-t*m+p) = {g2}")
    if g2 != 1:
        ok = False
        print("      NOT EXACT: a cofactor prime satisfies the Frobenius relation")
    else:
        print("      exact: no cofactor prime can satisfy the relation")

    u2, u3 = psi_constants(cv)
    passed_in = failed_out = 0
    for _ in range(2):
        R = random_point_fp2(cv)
        Q = R.mul(h2)                           # in G2
        if not Q.inf:
            if psi(cv, Q, u2, u3) == Q.mul(m2 % r):
                passed_in += 1
            else:
                ok = False
                print("      psi does NOT act as the claimed multiplier on G2")
        T = R.mul(r)                            # in the cofactor part
        if not T.inf:
            if psi(cv, T, u2, u3) != T.mul(m2 % r):
                failed_out += 1
            else:
                ok = False
                print("      REJECTED POINT ACCEPTED: a non-G2 point passed")
    print(f"      numerically: {passed_in} G2 points accepted,"
          f" {failed_out} cofactor points rejected")
    print()
    return ok


def main():
    random.seed(20260907)
    bad = 0
    for name in ("BLS12-381", "BLS12-461", "BN-462"):
        if not check_curve(name):
            bad += 1
    if bad:
        print(f"FAILED: {bad} curve(s) did not verify")
        return 1
    print("Every subgroup test is exact and confirmed on points outside the subgroup.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
