"""Local experiment: prescribed norms of a complete Miller function.

This is NOT a compressed Miller implementation or a demonstrated speedup.
It tests the invariant needed by the proposed state representation, with
vertical denominators retained. No production files are changed.
"""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "reference"))

from elips_ref import CURVES, Fp12
from gen_pairing_vectors import generators, FILES
from pairing_ref import untwist, fp_into_fp12, line


def full_miller(cv, Q, xp, yp, m):
    """Monic, denominator-retaining Miller function for a positive integer."""
    if m < 1:
        raise ValueError("positive m only")
    f, T = Fp12.one(cv.p), Q
    for bit in bin(m)[3:]:
        U = T.dbl()
        f = f.sqr() * line(T, T, xp, yp) * (xp - U.x).inv()
        T = U
        if bit == "1":
            U = T.add(Q)
            f = f * line(T, Q, xp, yp) * (xp - U.x).inv()
            T = U
    return f, T


def norm6(a):
    return a * a.conj()


def norm4(a):
    b = a.frob(4)
    return a * b * b.frob(4)


def run(name):
    cv = CURVES[name]
    G1, G2 = generators(cv, FILES[name][0])
    P, q = G1.mul(7), G2.mul(11)
    Q = untwist(cv, q.x, q.y)
    xp, yp = fp_into_fp12(P.x, cv.p), fp_into_fp12(P.y, cv.p)
    vQ, hQ = xp - Q.x, yp - Q.y
    kQ = xp.pow(3) - Q.x.pow(3)
    assert vQ.frob(6) == vQ
    assert hQ.frob(4) == hQ
    assert kQ.frob(2) == kQ
    for m in (2, 3, 5):
        f, T = full_miller(cv, Q, xp, yp, m)
        vT, hT = xp - T.x, yp - T.y
        kT = xp.pow(3) - T.x.pow(3)
        A = vQ.pow(m) * vT.inv()
        B = hQ.pow(m) * hT.inv()
        C = kQ.pow(m) * kT.inv()
        assert norm6(f) == A
        assert norm4(f) == B
        assert norm6(B) == C
        # Inclusion-exclusion of relative norms: maps into Phi_12 subgroup.
        g = f.pow(6) * C * (A.pow(3) * B.sqr()).inv()
        assert norm6(g).is_one()
        assert norm4(g).is_one()
        # Same projector, expressed by the usual easy exponent.
        if m == 2:
            easy = (cv.p**6 - 1) * (cv.p**2 + 1)
            assert g == f.pow(easy * (cv.p**2 - 2))
        print(name, "m =", m, "PASS prescribed norms and cyclotomic projection",
              flush=True)


if __name__ == "__main__":
    for name in ("BLS12-381", "BN-462"):
        run(name)
    print("No compressed recurrence or runtime improvement is asserted.")
