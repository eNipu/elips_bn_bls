"""Symbolically check the library's final exponentiation chains.

Every operation in a hard-part chain acts on a cyclotomic element as a fixed
exponentiation, so the whole chain computes m -> m^E for one integer E.
Tracking the exponent instead of the field element gives E exactly and
instantly, with no discrete logarithm needed.

  fp12_mul        exponents add
  fp12_sqr_cyc    exponent doubles
  fp12_exp_param  exponent multiplies by the signed mother parameter x
  fp12_frobenius  exponent multiplies by p^k
  fp12_conj       exponent negates

Conjugation is negation only because the easy part already put the element in
the cyclotomic subgroup, where the p^6 Frobenius equals inversion. Both chains
rely on that, and so does this trace.

The correct hard exponent is lambda = (p^4 - p^2 + 1)/r. The BLS12 chain
computes 3*lambda by design (see issue #16 and the comment in the C source);
the BN chain computes lambda exactly. This script asserts both, so a change to
either chain that alters its exponent fails here rather than only in the
numeric pairing vectors.

Exit status 0 if every chain traces to its claimed exponent, 1 otherwise.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES


def trace_bls12(p, x):
    """src/pairing/miller.c : pairing_final_exp_fast, ELIPS_FAMILY_BLS12.

    Claims 3*lambda = (x-1)^2 (x+p) (x^2+p^2-1) + 3.
    """
    m = 1                       # exponent of the easy-part output
    mi = -m                     # fp12_conj(mi, m)

    a = x * m                   # fp12_exp_param(a, m)
    a = a + mi                  # fp12_mul(a, a, mi)          -> m^(x-1)

    b = x * a                   # fp12_exp_param(b, a)
    b = b + (-a)                # fp12_mul(b, b, conj(a))     -> m^((x-1)^2)

    c = x * b                   # fp12_exp_param(c, b)
    c = c + p * b               # fp12_mul(c, c, frob(b,1))   -> b^(x+p)

    d = x * c                   # fp12_exp_param(d, c)
    d = x * d                   # fp12_exp_param(d, d)        -> c^(x^2)
    d = d + p**2 * c            # fp12_mul(d, d, frob(c,2))
    d = d + (-c)                # fp12_mul(d, d, conj(c))     -> c^(x^2+p^2-1)

    e = 2 * m + m               # fp12_sqr_cyc then fp12_mul  -> m^3
    return d + e                # fp12_mul(r, d, e)


def trace_bn(p, x):
    """src/pairing/miller.c : pairing_final_exp_fast, BN branch.

    Claims lambda = u3^(-36-36p) * u2^(-30-18p+6p^2) * u1^(-18-12p)
                    * m^(-2+p+p^2+p^3), with uk = m^(x^k).
    """
    m = 1
    u1 = x * m                  # fp12_exp_param(u1, m)
    u2 = x * u1                 # fp12_exp_param(u2, u1)
    u3 = x * u2                 # fp12_exp_param(u3, u2)

    # small_powers() builds a6, a12, a18, a30, a36; the chain then conjugates
    # and Frobenius-twists them, which is exactly the signed sum below.
    acc = (-36 - 36 * p) * u3
    acc += (-30 - 18 * p + 6 * p**2) * u2
    acc += (-18 - 12 * p) * u1
    acc += (-2 + p + p**2 + p**3) * m
    return acc


def main():
    bad = 0
    for name in ("BN-462", "BLS12-461", "BLS12-381"):
        cv = CURVES[name]
        p, r, x = cv.p, cv.r, cv.X

        assert (p**4 - p**2 + 1) % r == 0, "r must divide p^4-p^2+1"
        lam = (p**4 - p**2 + 1) // r

        if cv.family == "bn":
            E, want, claim = trace_bn(p, x), lam, "lambda"
        else:
            E, want, claim = trace_bls12(p, x), 3 * lam, "3*lambda"

        print(f"\n=== {name} ===")
        print(f"  lambda = (p^4-p^2+1)/r : {lam.bit_length()} bits")
        print(f"  chain traces to        : {E.bit_length()} bits")
        print(f"  chain claims           : {claim}")

        if E == want:
            print("  OK: exponent is exactly what the chain claims")
        else:
            bad += 1
            print("  MISMATCH")
            q, rem = divmod(E, lam)
            if rem == 0:
                print(f"    chain = lambda x {q}, not the claimed {claim}")
            else:
                k = (E % r) * pow(lam % r, -1, r) % r
                print(f"    chain is not an integer multiple of lambda;"
                      f" modulo r it is lambda x {k}")

        # Whatever multiple it is, it must stay a pairing: a factor sharing a
        # divisor with r would send some inputs to 1.
        k = (E // lam) if E % lam == 0 else ((E % r) * pow(lam % r, -1, r) % r)
        g = math.gcd(k % r, r)
        print(f"  gcd(multiplier, r) = {g}"
              f"  ({'non-degenerate' if g == 1 else 'DEGENERATE'})")
        if g != 1:
            bad += 1

    print()
    if bad:
        print(f"FAILED: {bad} problem(s)")
        return 1
    print("All final exponentiation chains trace to their claimed exponents.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
