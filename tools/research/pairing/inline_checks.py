"""Reproducible versions of earlier shell-only algebra checks.

Requires SymPy for exact polynomial arithmetic. Failed heuristic searches and
cost estimates are not impossibility proofs; see CORRECTIONS.md.
"""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "reference"))
import sympy as s
from elips_ref import CURVES
from gen_pairing_vectors import generators, FILES
from pairing_ref import untwist, fp_into_fp12, line


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    for name, cv in CURVES.items():
        E, remainder = divmod(cv.p**12-1, cv.r)
        check(remainder == 0, name)
        for d in (1, 2, 4, 6):
            check(E % (cv.p**d-1) == 0, f"{name}: subfield {d}")
        print("PASS subfield scaling:", name)

    # Preserve the earlier bounded polynomial-space search, not a circuit bound.
    x = s.symbols("x")
    r = x**4-x*x+1
    p = s.expand((x-1)**2*r/3+x)
    lam = s.cancel((p**4-p*p+1)/r)
    check(s.expand(3*lam-(x-1)**2*(x+p)*(x*x+p*p-1)-3) == 0,
          "BLS12 hard exponent identity")
    for cap in (3, 4):
        basis = [s.expand(x**i*p**j) for j in range(4) for i in range(cap+1)]
        degree = max(s.degree(b, x) for b in basis)
        multiplier_degree = degree-s.degree(lam, x)
        polys = basis+[-s.expand(lam*x**i) for i in range(max(0, multiplier_degree+1))]
        matrix = s.Matrix([[b.coeff(x, k) for b in polys] for k in range(degree+1)])
        print(f"Bounded polynomial search x-degree <= {cap}: nullity {len(matrix.nullspace())}")

    # The norm/torus cubic relation derived in the earlier experiment.
    a, b, c, v, z, xi = s.symbols("a b c v z xi")
    t = a+b*v+c*v*v
    sig = a+z*b*v+z*z*c*v*v
    tau = a+z*z*b*v+z*c*v*v
    expr = s.rem(s.expand(-z*z*t*tau-z*t*sig-sig*tau-v), z*z+z+1, z)
    expr = s.rem(expr, v**3-xi, v)
    check(s.expand(expr-v*(3*a*b-3*c*c*xi-1)) == 0, "cubic norm relation")
    print("PASS cubic norm relation")

    # Both new interpolation identities as exact polynomials, including reduction.
    eta = s.symbols("eta")
    pp, pm = (a+b+c)**2, (a-b+c)**2
    pu, mu = (a+s.I*b-c)**2, (a-s.I*b-c)**2
    S, D, total = pp+pm, pp-pm, a*a+c*c
    J = -s.I*(2*pu+S-4*total)
    gaussian = [4*a*a+eta*(D-J), D+J+4*eta*c*c, 2*S-4*total]
    T, V = pu+mu, -s.I*(pu-mu)
    fourier = [S+T-4*c*c+eta*(D-V), D+V+4*eta*c*c, S-T]
    expected = s.rem(4*(a+b*v+c*v*v)**2, v**3-eta, v)
    for name, values in (("Gaussian", gaussian), ("Fourier", fourier)):
        check(all(s.expand(values[k]-expected.coeff(v, k)) == 0 for k in range(3)),
              name)
        print("PASS exact interpolation:", name)

    # Preserve the numerical support-growth check, without inferring a rank bound.
    cv = CURVES["BLS12-381"]
    G1, G2 = generators(cv, FILES[cv.name][0])
    P, q = G1.mul(3), G2.mul(5)
    Q = untwist(cv, q.x, q.y)
    xp, yp = fp_into_fp12(P.x, cv.p), fp_into_fp12(P.y, cv.p)
    L0 = line(Q, Q, xp, yp)
    Q2 = Q.dbl()
    L1 = line(Q2, Q2, xp, yp)
    def support(f):
        return sum(not t.is_zero() for half in (f.d0, f.d1)
                   for t in (half.c0, half.c1, half.c2))
    found = [support(f) for f in (L0, L0.sqr(), L0*L1, L0.sqr()*L1)]
    check(found == [3, 5, 5, 6], "line support control")
    print("PASS line supports L0, L0^2, L0*L1, L0^2*L1:", found)
    n, ternary = abs(cv.X), []
    while n:
        ternary.append(n % 3)
        n //= 3
    print("Seed: binary weight", abs(cv.X).bit_count(), "ternary length",
          len(ternary), "ternary weight", sum(d != 0 for d in ternary))
    print("These checks establish their explicit identities, not novelty or optimality.")


if __name__ == "__main__":
    main()
