"""Symbolically trace the library's final exponentiation chains.

Every operation in the chain (square, multiply, Frobenius, conjugate) acts on a
cyclotomic element as a fixed exponentiation, so the whole chain computes
f -> f^E' for one integer E'. Tracking the exponent instead of the field element
gives E' exactly and instantly, with no discrete logarithm needed.

After the easy part the element lies in the cyclotomic subgroup, where the p^6
Frobenius equals inversion. The trace therefore uses -1 for frobenius_p6 in the
hard part, matching what the code relies on.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES


def trace_bn(p, X):
    """src/bn_final_exp.c : bn_final_exp_optimal, hard part only."""
    A  = 1                       # exponent of the (already easy-part-reduced) input
    t0 = -A                      # frobenius_p6 then... see below
    # t0 = A^X ; then frobenius_p6 -> A^(-X)
    t0 = -X * A
    t1 = 2 * t0                  # Fp12_sqr(t1, t0)          -> -2X
    t0b = 2 * t1                 # Fp12_sqr(t0, t1)          -> -4X
    t0 = t1 + t0b                # Fp12_mul(t0, t1, t0)      -> -6X
    t2 = X * t0                  # pow_motherparam           -> -6X^2
    t2 = -t2                     # frobenius_p6              -> +6X^2
    t3 = (X + 1) * t2            # pow_motherparam with X_binary[0]=0 -> ^(X+1)
    t3 = 2 * t3                  # Fp12_sqr
    t0 = -t0                     # frobenius_p6              -> +6X
    t3 = t3 + t0                 # Fp12_mul(t3, t3, t0)
    t2 = -t2                     # frobenius_p6              -> -6X^2
    t2 = t3 + t2                 # Fp12_mul(t2, t3, t2)
    t3 = t3 + A                  # Fp12_mul(t3, t3, A)
    t1 = t1 + t2                 # Fp12_mul(t1, t1, t2)
    t0 = -A                      # frobenius_p6(A)           -> A^-1
    t0 = t1 + t0                 # Fp12_mul(t0, t1, t0)
    t0 = p**3 * t0               # frobenius_p3
    t0 = t0 + t3
    t1 = p * t1                  # frobenius_p1
    t0 = t0 + t1
    t2 = p**2 * t2               # frobenius_p2
    t0 = t0 + t2
    return t0


def trace_bls12(p, X):
    """src/bls12_finalexp.c : bls12_finalexp_optimal, hard part only."""
    A  = 1
    t0 = 2 * A                          # Fp12_sqr(t0, A)
    t1 = X * t0                         # pow_motherparam
    t2 = (X // 2) * t1                  # pow with the halved mother parameter
    t3 = -A                             # frobenius_p6(A) -> A^-1
    t1 = t3 + t1                        # Fp12_mul(t1, t3, t1)
    t1 = -t1                            # frobenius_p6
    t1 = t1 + t2                        # Fp12_mul(t1, t1, t2)
    t2 = X * t1                         # pow_motherparam
    t3 = X * t2                         # pow_motherparam
    t1 = -t1                            # frobenius_p6
    t3 = t1 + t3                        # Fp12_mul(t3, t1, t3)
    t1 = -t1                            # frobenius_p6
    t1 = p**3 * t1                      # frobenius_p3
    t2 = p**2 * t2                      # frobenius_p2
    t1 = t1 + t2                        # Fp12_mul(t1, t1, t2)
    t2 = X * t3                         # pow_motherparam
    t2 = t2 + t0                        # Fp12_mul(t2, t2, t0)
    t2 = t2 + A                         # Fp12_mul(t2, t2, A)
    t1 = t1 + t2                        # Fp12_mul(t1, t1, t2)
    t2 = p * t3                         # frobenius_p1
    return t1 + t2


for name in ("BN-462", "BLS12-461", "BLS12-381"):
    cv = CURVES[name]
    p, r, X = cv.p, cv.r, cv.X
    hard_true = (p**4 - p**2 + 1) // r
    assert (p**4 - p**2 + 1) % r == 0, "r must divide p^4-p^2+1"
    E = trace_bn(p, X) if cv.family == "bn" else trace_bls12(p, X)

    print(f"\n=== {name} ===")
    print(f"  correct hard exponent (p^4-p^2+1)/r : {hard_true.bit_length()} bits")
    print(f"  library chain computes              : {E.bit_length()} bits")
    q, rem = divmod(E, hard_true)
    if rem == 0:
        print(f"  library = correct x {q}   -> pairing^{q}")
        import math
        print(f"  gcd({q}, r) = {math.gcd(q % r, r)}  "
              f"({'still a valid pairing' if math.gcd(q % r, r)==1 else 'DEGENERATE'})")
    else:
        # not an exact integer multiple; compare modulo r, which is what matters
        k = (E % r) * pow(hard_true % r, -1, r) % r
        print(f"  NOT an exact integer multiple of the correct exponent")
        print(f"  but modulo r: library = correct x {k}")
        if k < 10**6:
            print(f"  -> library computes pairing^{k}")
