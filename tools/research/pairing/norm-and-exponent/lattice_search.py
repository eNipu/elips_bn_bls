"""Preserved reconstruction of the earlier inline LLL/Babai experiment.

This is a HEURISTIC coset search, not an optimality or nonexistence proof.
The default runs quick unit tests. --search enables the slow exact-Fraction
research computation. Failed basis-row/target-embedding approaches are retained
as constructors for inspection; they are not valid negative tests.
"""
import argparse
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "reference"))
from elips_ref import CURVES


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def nearest_integer(q):
    return (2*q.numerator + q.denominator) // (2*q.denominator)


def gram_schmidt(basis):
    orthogonal, norms, mu = [], [], []
    for row in basis:
        v = [Fraction(x) for x in row]
        coefficients = []
        for star, norm in zip(orthogonal, norms):
            check(norm != 0, "dependent basis")
            m = sum(x*y for x, y in zip(row, star)) / norm
            coefficients.append(m)
            v = [x-m*y for x, y in zip(v, star)]
        orthogonal.append(v)
        norms.append(sum(x*x for x in v))
        mu.append(coefficients)
    return orthogonal, norms, mu


def lll(basis, delta=Fraction(3, 4)):
    """Exact arithmetic, intentionally simple and slow, as in the experiment."""
    b = [row[:] for row in basis]
    _, norms, mu = gram_schmidt(b)
    k, swaps = 1, 0
    while k < len(b):
        for j in range(k-1, -1, -1):
            q = nearest_integer(mu[k][j])
            if q:
                b[k] = [x-q*y for x, y in zip(b[k], b[j])]
                _, norms, mu = gram_schmidt(b)
        if norms[k] >= (delta-mu[k][k-1]**2)*norms[k-1]:
            k += 1
        else:
            b[k], b[k-1] = b[k-1], b[k]
            swaps += 1
            _, norms, mu = gram_schmidt(b)
            k = max(k-1, 1)
    return b, swaps


def babai(basis, target):
    star, norms, _ = gram_schmidt(basis)
    residual = [Fraction(x) for x in target]
    lattice_vector = [0]*len(target)
    for i in range(len(basis)-1, -1, -1):
        q = nearest_integer(sum(x*y for x, y in zip(residual, star[i])) / norms[i])
        residual = [x-q*y for x, y in zip(residual, basis[i])]
        lattice_vector = [x+q*y for x, y in zip(lattice_vector, basis[i])]
    return lattice_vector


def kernel_basis(monomials, modulus):
    pivot = next(i for i, m in enumerate(monomials) if gcd(m, modulus) == 1)
    inverse = pow(monomials[pivot], -1, modulus)
    rows = []
    for i, m in enumerate(monomials):
        if i != pivot:
            row = [0]*len(monomials)
            row[pivot] = (-m*inverse) % modulus
            row[i] = 1
            rows.append(row)
    row = [0]*len(monomials)
    row[pivot] = modulus
    rows.append(row)
    return rows, pivot, inverse


def old_embedding(monomials, modulus, target=None, W=1 << 20, K=1):
    """Retain the abandoned basis-row/target-embedding construction.

    Merely inspecting LLL basis rows cannot certify that a representation is
    absent. A target row also permits multiples of the target. The earlier
    attempts varied W and K but did not establish a closest-vector bound.
    """
    m = len(monomials)
    width = m + 1 + (target is not None)
    rows = []
    for i, value in enumerate(monomials):
        row = [0]*width
        row[i], row[m] = W, value
        rows.append(row)
    row = [0]*width
    row[m] = modulus
    rows.append(row)
    if target is not None:
        row = [0]*width
        row[m], row[m+1] = target, K
        rows.append(row)
    return rows


def search(max_x, max_p):
    cv = CURVES["BLS12-381"]
    N = cv.p**4-cv.p**2+1
    target = 3*(N//cv.r)
    x = abs(cv.X)
    # Match the original order: nonconstant monomials first, then 1.
    names = [(i, j) for j in range(max_p+1) for i in range(max_x+1) if i or j]
    names.append((0, 0))
    monomials = [pow(x, i, N)*pow(cv.p, j, N) % N for i, j in names]
    rows, pivot, inverse = kernel_basis(monomials, N)
    print(f"Reducing {len(rows)} rows; exact Fraction arithmetic can be slow.", flush=True)
    reduced, swaps = lll(rows)
    particular = [0]*len(monomials)
    particular[pivot] = target*inverse % N
    correction = babai(reduced, particular)
    coefficients = [a-b for a, b in zip(particular, correction)]
    check(sum(a*b for a, b in zip(coefficients, monomials)) % N == target,
          "coset invariant")
    print("LLL swaps:", swaps)
    print("Candidate max|coefficient|:", max(abs(c) for c in coefficients))
    print("Terms:", [(name, c) for name, c in zip(names, coefficients) if c])
    print("This is one heuristic candidate, NOT a shortest-vector certificate.")
    if max_x == 5 and max_p == 3:
        check(max(abs(c) for c in coefficients) <= 3, "degree-five positive control")


def selftest():
    reduced, _ = lll([[4, 1], [3, 2]])
    determinant = reduced[0][0]*reduced[1][1]-reduced[0][1]*reduced[1][0]
    check(abs(determinant) == 5, "lattice determinant")
    check(any(sum(x*x for x in row) == 2 for row in reduced), "short-vector control")
    _, norms, mu = gram_schmidt(reduced)
    for i in range(1, len(reduced)):
        check(all(abs(c) <= Fraction(1, 2) for c in mu[i]), "size reduction")
        check(norms[i] >= (Fraction(3, 4)-mu[i][i-1]**2)*norms[i-1], "Lovasz")
    mons, N, target = [13, 57, 1], 997, 197
    rows, pivot, inverse = kernel_basis(mons, N)
    check(all(sum(a*b for a, b in zip(row, mons)) % N == 0 for row in rows),
          "kernel basis")
    reduced, _ = lll(rows)
    particular = [0, 0, 0]
    particular[pivot] = target*inverse % N
    c = [a-b for a, b in zip(particular, babai(reduced, particular))]
    check(sum(a*b for a, b in zip(c, mons)) % N == target, "small coset")
    check(old_embedding([13, 57], 997, 197, W=256)[-1] == [0, 0, 197, 1],
          "archived embedding")
    cv = CURVES["BLS12-381"]
    x, p = cv.X, cv.p
    N = p**4-p**2+1
    check((x-1)**2*(x+p)*(x*x+p*p-1)+3 == 3*(N//cv.r), "known exact chain")
    print("PASS LLL/Babai unit tests and exact known-chain control")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--max-x", type=int, choices=(3, 4, 5), default=5)
    parser.add_argument("--max-p", type=int, choices=(3, 4, 6), default=3)
    args = parser.parse_args()
    selftest()
    if args.search:
        search(args.max_x, args.max_p)
