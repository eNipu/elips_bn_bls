"""Executable examples for RESEARCH_NOTE.md, using only the standard library.

Research/teaching code, not a production pairing interface. All inputs are
public and deterministic. --skip-pairing runs only the short arithmetic examples.
"""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools/reference"))
sys.path.insert(0, str(HERE / "scaled-square"))

from elips_ref import Fp2, Fp6, Fp12
import normalized_miller_ref as nm
from gaussian_ref import fft4_square, replay, scale12


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def toy_examples():
    """Hand-sized arithmetic in F13[z]/(z^3-2), not a cryptographic curve."""
    p, u, eta = 13, 5, 2
    check(u*u % p == p-1, "u is not a square root of -1")
    check(eta not in {pow(t, 3, p) for t in range(p)}, "cubic is reducible")

    a, b, c = 1, 2, 3
    B, C = 4, 5
    m0, m3 = a*B % p, c*C % p
    mp = (a+b+c)*(B+C) % p
    mm = (a-b+c)*(B-C) % p
    H = [2*(m0+eta*m3) % p, (mp-mm-2*m3) % p, (mp+mm-2*m0) % p]
    # Independent schoolbook convolution and reduction.
    product = [a*B, a*C+b*B, b*C+c*B, c*C]
    expected = [2*(product[0]+eta*product[3]) % p,
                2*product[1] % p, 2*product[2] % p]
    check([m0, m3, mp, mm] == [4, 2, 2, 11], "line example intermediates")
    check(H == expected == [3, 0, 5], "line interpolation example")
    print("PASS toy line product: 2*T(z)*(4+5z) = 3+5z^2")

    P, M = (a+b+c)**2 % p, (a-b+c)**2 % p
    I, J = (a+u*b-c)**2 % p, (a-u*b-c)**2 % p
    K = 4*c*c % p
    S, T, U, V = (P+M) % p, (I+J) % p, (P-M) % p, -u*(I-J) % p
    h = [(S+T-K+eta*(U-V)) % p, (U+V+eta*K) % p, (S-T) % p]
    square = [a*a, 2*a*b, b*b+2*a*c, 2*b*c, c*c]
    expected = [4*(square[0]+eta*square[3]) % p,
                4*(square[1]+eta*square[4]) % p, 4*square[2] % p]
    check([P, M, I, J, K] == [10, 4, 12, 1, 10], "Fourier evaluations")
    check([S, T, U, V] == [1, 0, 6, 10], "Fourier combinations")
    check(h == expected == [9, 10, 1], "Fourier example")
    print("PASS toy Fourier square: 4*(1+2z+3z^2)^2 = 9+10z+z^2")
    # A deliberately corrupted coordinate must fail the same comparison.
    corrupted = h[:]
    corrupted[0] = (corrupted[0]+1) % p
    check(corrupted != expected, "corruption control")


def tower_examples():
    p = nm.P_MOD
    for d in (1, 2, 4, 6):
        check(nm.FINAL_EXP % (p**d-1) == 0, "subfield divisibility")

    def scalar(n):
        return Fp2(p, n % p, 0)

    f = Fp12(Fp6(scalar(1), scalar(2), scalar(3)),
             Fp6(scalar(4), scalar(5), scalar(6)))
    zero = Fp2.zero(p)
    B, C = Fp2(p, 7, 8), Fp2(p, 9, 10)
    line = Fp12(Fp6(Fp2.one(p), zero, zero), Fp6(zero, B, C))
    check(nm.mul_line_scaled8(f, B, C) == scale12(f*line, 2),
          "BLS12-381 8M2 line identity")
    check(fft4_square(f) == scale12(f.sqr(), 4), "BLS12-381 Fourier identity")
    combined = nm.mul_line_scaled8(fft4_square(f), B, C)
    check(combined == scale12(f.sqr()*line, 8), "combined one-line update")
    # This is intentionally NOT equality of raw values.
    check(fft4_square(Fp12.one(p)) != Fp12.one(p), "raw-value warning")
    print("PASS BLS12-381: scaled line = 2*f*L; Fourier square = 4*f^2")
    print("PASS combined one-line update = 8*f^2*L")


def pairing_example():
    # The loader returns public, canonical committed test vectors.
    P, Q, expected = nm._load_vectors(nm.DEFAULT_VECTORS)[0]
    table = nm.prepare_g2(Q)
    baseline = nm.miller_product_scaled8([(P, table)])  # validates P
    candidate, scale = replay([(P, table)], fft4_square)
    check(candidate == scale12(baseline, scale), "whole-loop scale relation")
    e = nm.final_exponentiate(candidate)
    check(e == expected, "independent known-answer vector")
    check(not e.is_one(), "nontrivial known-answer control")
    cubed = nm.final_exponentiate(candidate, cubed=True)
    check(cubed == e.pow(3), "e versus e^3 convention")
    cancel, _ = replay([(P, table), (-P, table)], fft4_square)
    check(nm.final_exponentiate(cancel).is_one(), "pairing cancellation")
    # Replay is a trusted-input experiment, not a validated public API.
    print("PASS complete pairing: raw scale relation, exact KAT, and e^3 convention")
    print("PASS cancellation: e(P,Q)*e(-P,Q)=1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pairing", action="store_true")
    args = parser.parse_args()
    toy_examples()
    tower_examples()
    if not args.skip_pairing:
        pairing_example()
    print("All requested tutorial examples passed.")


if __name__ == "__main__":
    main()
