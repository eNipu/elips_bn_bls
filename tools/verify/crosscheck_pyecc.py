#!/usr/bin/env python3
"""Check the BLS12-381 pairing against an implementation written by someone else.

WHY THIS EXISTS

Every other test in this repository compares the C against tools/reference/,
and both were written by the same hand. A shared misconception passes both in
silence. That is not hypothetical here: issue #16 was a final exponentiation
raising to the wrong exponent, and it survived for years because the result was
still bilinear, still non-degenerate, and still agreed with itself.

py_ecc is the Ethereum Foundation's BLS12-381 implementation. It shares no code,
no author and no representation with this library: its Fp12 is a direct
degree-12 extension, not a tower of Fp2 and Fp6.

THE BASIS CHANGE

The two representations are reconcilable, which is the only reason this
comparison is possible. ELiPS builds Fp12 as Fp2[w]/(w^6 - xi) with xi = 1 + u
and u^2 = -1, so

    u = w^6 - 1   =>   (w^6 - 1)^2 = -1   =>   w^12 - 2w^6 + 2 = 0

which is exactly py_ecc's modulus polynomial. An element of ELiPS is
sum_i g_i w^i with g_i = c0 + c1 u in Fp2, and substituting u = w^6 - 1 gives

    g_i w^i = (c0 - c1) w^i + c1 w^(i+6)

so the twelve tower coefficients map onto py_ecc's twelve w-basis coefficients.
tools/verify/identities.gp proves the modulus polynomial claim symbolically.

WHAT IS COMPARED, AND WHAT IS NOT

The EXACT pairing, from pairing_final_exp_plain, not elips_pairing (which on
BLS12 returns e^3 by design). ELiPS's value comes out as py_ecc's INVERSE: a
pure orientation convention, and both are valid non-degenerate pairings. That
inversion is expected and is asserted, not tolerated.

Only BLS12-381. py_ecc implements no other curve, so BLS12-461 and BN-462 get
no anchor from this file. tools/verify/curve_params.gp is what covers them.

Run:  python3 tools/verify/crosscheck_pyecc.py <path-to-dump-binary>
"""
import subprocess
import sys


def load_elips(binary):
    """Run the dump helper and parse its output."""
    out = subprocess.run([binary], capture_output=True, text=True, check=True).stdout
    pt, coeff = {}, {}
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "e":
            coeff[(int(t[1]), int(t[2]), int(t[3]))] = int(t[4])
        else:
            pt[t[0]] = int(t[1])
    return pt, coeff


def to_w_basis(coeff, p):
    """ELiPS tower coefficients -> py_ecc's basis in powers of w.

    a[0] holds (g0, g2, g4) and a[1] holds (g1, g3, g5), by powers of w.
    """
    slot = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3, (0, 2): 4, (1, 2): 5}
    out = [0] * 12
    for (d, c, b), v in coeff.items():
        i = slot[(d, c)]
        if b == 0:
            out[i] = (out[i] + v) % p
        else:                       # the u part, and u = w^6 - 1
            out[i] = (out[i] - v) % p
            out[i + 6] = (out[i + 6] + v) % p
    return out


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__.strip().splitlines()[-1])
    try:
        from py_ecc.optimized_bls12_381 import G1, G2, pairing, field_modulus, normalize
        from py_ecc.fields import optimized_bls12_381_FQ12 as FQ12
    except ImportError:
        sys.exit("py_ecc is not installed. See tools/verify/README.md")

    p = field_modulus
    pt, coeff = load_elips(sys.argv[1])
    fails = 0

    def check(name, cond):
        nonlocal fails
        print(f"  {'ok   ' if cond else 'FAIL '} {name}")
        if not cond:
            fails += 1

    print("BLS12-381 against py_ecc (independent implementation)\n")

    # If the two sides disagree on the inputs, agreeing on the output would
    # mean nothing.
    g1 = normalize(G1)
    g2 = normalize(G2)
    check("same G1 generator", int(g1[0]) == pt["P.x"] and int(g1[1]) == pt["P.y"])
    check("same G2 generator",
          g2[0].coeffs[0] == pt["Q.x0"] and g2[0].coeffs[1] == pt["Q.x1"]
          and g2[1].coeffs[0] == pt["Q.y0"] and g2[1].coeffs[1] == pt["Q.y1"])

    mine = to_w_basis(coeff, p)
    theirs = pairing(G2, G1)
    inverse = [int(v) % p for v in (FQ12.one() / theirs).coeffs]
    direct = [int(v) % p for v in theirs.coeffs]

    check("exact pairing equals py_ecc's, all 12 coefficients (up to orientation)",
          mine == inverse)
    check("the orientation really is inversion, not equality", mine != direct)

    # A cross-check that cannot fail proves nothing. e^3 is what issue #16
    # actually computed, and it is still bilinear and still non-degenerate.
    wrong = {
        "e^3, the issue #16 bug": theirs ** 3,
        "e^-3": (FQ12.one() / theirs) ** 3,
        "e^2": theirs ** 2,
    }
    for name, val in wrong.items():
        check(f"rejects {name}", mine != [int(v) % p for v in val.coeffs])

    print(f"\n{fails} failure(s)" if fails else "\nBLS12-381 pairing agrees with py_ecc")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
