"""Reference optimal ate pairing (issue #13).

Written as plain textbook Miller's algorithm carried out entirely in Fp12 on
untwisted points. That is far slower than the library's sparse-multiplication
approach, and deliberately so: an oracle should share no machinery with the code
it checks. No twist shortcuts, no precomputed Frobenius constants, no sparse
line representation, no addition chain for the final exponent.

Untwisting, derived from the library's own map in src/bn_twist.c:
    x_E12 = (x'/xi) * v^2,  y_E12 = (y'/xi) * v*w
with v = w^2 and xi = v^3 = w^6 this is
    x = x' * w^-2,  y = y' * w^-3
so the twist is D-type and E'(Fp2): y'^2 = x'^3 + b*xi.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, Fp2, Fp6, Fp12


def _w(p):
    """The element w of Fp12 = Fp6[w]/(w^2 - v)."""
    return Fp12(Fp6.zero(p), Fp6.one(p))


def fp2_into_fp12(a):
    """Embed an Fp2 element as the constant term of Fp12."""
    p = a.p
    return Fp12(Fp6(a, Fp2.zero(p), Fp2.zero(p)), Fp6.zero(p))


def fp_into_fp12(x, p):
    return fp2_into_fp12(Fp2(p, x % p, 0))


class Pt12:
    """Affine point on E(Fp12): y^2 = x^3 + b_signed. None means infinity."""
    __slots__ = ("x", "y", "inf", "p")

    def __init__(self, x=None, y=None, p=None):
        self.x, self.y = x, y
        self.inf = x is None
        self.p = p if p is not None else (x.p if x is not None else None)

    def __neg__(self):
        return Pt12(p=self.p) if self.inf else Pt12(self.x, -self.y, self.p)

    def dbl(self):
        if self.inf:
            return Pt12(p=self.p)
        lam = (self.x.sqr() + self.x.sqr() + self.x.sqr()) * (self.y + self.y).inv()
        xr = lam.sqr() - self.x - self.x
        return Pt12(xr, lam * (self.x - xr) - self.y, self.p)

    def add(self, o):
        if self.inf: return o
        if o.inf:    return self
        if (self.x - o.x).d0.is_zero() and (self.x - o.x).d1.is_zero():
            same_y = (self.y - o.y).d0.is_zero() and (self.y - o.y).d1.is_zero()
            return self.dbl() if same_y else Pt12(p=self.p)
        lam = (o.y - self.y) * (o.x - self.x).inv()
        xr = lam.sqr() - self.x - o.x
        return Pt12(xr, lam * (self.x - xr) - self.y, self.p)

    def frob(self, k=1):
        """(x, y) -> (x^(p^k), y^(p^k)); maps E(Fp12) to itself since b is in Fp."""
        if self.inf:
            return self
        return Pt12(self.x.frob(k), self.y.frob(k), self.p)


def untwist(cv, xq, yq):
    """Lift a twist point (x', y') in E'(Fp2) to E(Fp12)."""
    p = cv.p
    w = _w(p)
    winv = w.inv()
    w2i = winv.sqr()
    w3i = w2i * winv
    return Pt12(fp2_into_fp12(xq) * w2i, fp2_into_fp12(yq) * w3i, p)


def line(A, B, xp, yp):
    """Line through A and B (tangent if equal), evaluated at P = (xp, yp).

    xp and yp are Fp12 embeddings of Fp coordinates. Returns an Fp12 element.
    """
    p = A.p
    if A.inf or B.inf:
        return Fp12.one(p)
    dx = A.x - B.x
    if dx.d0.is_zero() and dx.d1.is_zero():
        dy = A.y - B.y
        if not (dy.d0.is_zero() and dy.d1.is_zero()):
            return xp - A.x                      # vertical line
        if A.y.d0.is_zero() and A.y.d1.is_zero():
            return xp - A.x
        lam = (A.x.sqr() + A.x.sqr() + A.x.sqr()) * (A.y + A.y).inv()
    else:
        lam = (B.y - A.y) * (B.x - A.x).inv()
    return (yp - A.y) - lam * (xp - A.x)


def miller(cv, Q, xp, yp, digits):
    """f_{m,Q}(P) for m given as {bit index: +-1}, plus the final T."""
    p = cv.p
    top = max(digits)
    f = Fp12.one(p)
    negQ = -Q
    T = Q if digits[top] > 0 else negQ          # seed with the leading digit
    for i in range(top - 1, -1, -1):
        f = f.sqr() * line(T, T, xp, yp)
        T = T.dbl()
        d = digits.get(i, 0)
        if d == 1:
            f = f * line(T, Q, xp, yp)
            T = T.add(Q)
        elif d == -1:
            f = f * line(T, negQ, xp, yp)
            T = T.add(negQ)
    return f, T


def optimal_ate(cv, px, py, qx, qy):
    """Optimal ate pairing. P = (px, py) in E(Fp); Q = (qx, qy) on the twist."""
    p = cv.p
    xp, yp = fp_into_fp12(px, p), fp_into_fp12(py, p)
    Q = untwist(cv, qx, qy)

    f, T = miller(cv, Q, xp, yp, cv.loop_digits)

    if cv.family == "bn":
        # BN needs two correction lines, matching src/bn_miller_optate.c
        Q1 = Q.frob(1)
        Q2 = Q.frob(2)
        f = f * line(T, Q1, xp, yp)
        T = T.add(Q1)
        f = f * line(T, -Q2, xp, yp)
        T = T.add(-Q2)

    # Final exponentiation, computed directly rather than by addition chain.
    e = (p**12 - 1) // cv.r
    return f.pow(e)


if __name__ == "__main__":
    import json
    src = sys.argv[1]
    data = json.load(open(src))
    cv = CURVES[data["curve"]]
    p = cv.p
    got = optimal_ate(cv, int(data["px"], 16), int(data["py"], 16),
                      Fp2(p, int(data["qx0"], 16), int(data["qx1"], 16)),
                      Fp2(p, int(data["qy0"], 16), int(data["qy1"], 16)))
    from elips_ref import fp12_to_list
    ref = fp12_to_list(got)
    lib = [int(v, 16) % p for v in data["f"]]
    ok = ref == lib
    print(f"curve         : {data['curve']}")
    print(f"reference == library : {ok}")
    if not ok:
        for i, (a, b) in enumerate(zip(ref, lib)):
            if a != b:
                print(f"  coord {i:2d}  ref={hex(a)[:40]}...  lib={hex(b)[:40]}...")
        # Is the library value a fixed power of ours? That would say the
        # algorithms agree but a convention (inverse, conjugate) differs.
        from elips_ref import fp12_from_list
        L = fp12_from_list(p, lib)
        for nm, cand in (("reference^-1", got.inv()),
                         ("conj(reference)", got.conj()),
                         ("reference^p", got.frob(1))):
            if fp12_to_list(cand) == lib:
                print(f"  NOTE: library value equals {nm}")
    print(f"in mu_r       : {got.pow(cv.r).is_one()}")
    print(f"non-degenerate: {not got.is_one()}")
    sys.exit(0 if ok else 1)
