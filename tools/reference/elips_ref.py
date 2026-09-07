"""
Independent reference implementation of the ELiPS field tower and curve
arithmetic, written from the defining equations rather than from the C source.

Purpose: serve as the Phase 0 oracle. If this file and src/ agree on a value,
two independently written implementations agree, which is evidence the value is
correct. If they disagree, one of them is wrong and the disagreement is the
finding.

Tower (must match include/elips/fpx.h and the *_mul_xi
routines in src/):

    Fp2  = Fp[u]  / (u^2 + 1)          element a0 + a1*u
    Fp6  = Fp2[v] / (v^3 - (1+u))      element c0 + c1*v + c2*v^2
    Fp12 = Fp6[w] / (w^2 - v)          element d0 + d1*w

so xi = 1 + u is the cubic non-residue, and v = w^2, xi = v^3 = w^6.

Curves:
    BN-462     y^2 = x^3 - 4     (src/bn_efp.c:EFp_rational_point_bn subtracts b)
    BLS12-461  y^2 = x^3 + 4     (src/bn_efp.c:EFp_rational_point_bls12 adds b)
    BLS12-381  y^2 = x^3 + 4     (standard)
"""

from dataclasses import dataclass


# ---------------------------------------------------------------- parameters

class Curve:
    def __init__(self, name, X, p, r, b, sign_b, loop_digits, family):
        self.name = name
        self.X = X
        self.p = p
        self.r = r
        self.b = b              # magnitude of the curve constant
        self.sign_b = sign_b    # +1 for y^2=x^3+b, -1 for y^2=x^3-b
        self.loop_digits = loop_digits   # {index: +-1} for the Miller loop
        self.family = family    # "bn" or "bls12"

    @property
    def b_signed(self):
        return (self.sign_b * self.b) % self.p


def _bn462():
    x = 2**114 + 2**101 - 2**14 - 1
    p = 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1
    r = 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1
    # src/curve_settings.c generate_bn_mother_parameter(): X_binary_opt = 6x+2
    loop = {116: 1, 115: 1, 103: 1, 102: 1, 16: -1, 15: -1, 2: -1}
    assert sum(s * 2**i for i, s in loop.items()) == 6*x + 2
    return Curve("BN-462", x, p, r, 4, -1, loop, "bn")


def _bls12_461():
    X = -2**77 + 2**50 + 2**33
    p = ((X-1)**2 * (X**4 - X**2 + 1)) // 3 + X
    r = X**4 - X**2 + 1
    loop = {77: -1, 50: 1, 33: 1}            # src/curve_settings.c bls12_generate_X
    assert sum(s * 2**i for i, s in loop.items()) == X
    return Curve("BLS12-461", X, p, r, 4, +1, loop, "bls12")


def _bls12_381():
    X = -0xd201000000010000
    p = ((X-1)**2 * (X**4 - X**2 + 1)) // 3 + X
    r = X**4 - X**2 + 1
    loop = _naf(X)
    assert sum(s * 2**i for i, s in loop.items()) == X
    return Curve("BLS12-381", X, p, r, 4, +1, loop, "bls12")


def _naf(n):
    """Non-adjacent form as {bit index: +-1}."""
    neg = n < 0
    n = abs(n)
    out, i = {}, 0
    while n:
        if n & 1:
            z = 2 - (n % 4)
            out[i] = -z if neg else z
            n -= z
        n >>= 1
        i += 1
    return out


CURVES = {c.name: c for c in (_bn462(), _bls12_461(), _bls12_381())}


# ------------------------------------------------------------------- Fp / Fp2
# Every element carries its own modulus so several curves can coexist in one
# process, which the C library cannot currently do (see plan section 4.4).

@dataclass(frozen=True)
class Fp2:
    p: int
    a: int      # coefficient of 1
    b: int      # coefficient of u

    @staticmethod
    def zero(p):  return Fp2(p, 0, 0)
    @staticmethod
    def one(p):   return Fp2(p, 1, 0)

    def __add__(s, o): return Fp2(s.p, (s.a+o.a) % s.p, (s.b+o.b) % s.p)
    def __sub__(s, o): return Fp2(s.p, (s.a-o.a) % s.p, (s.b-o.b) % s.p)
    def __neg__(s):    return Fp2(s.p, (-s.a) % s.p, (-s.b) % s.p)

    def __mul__(s, o):
        # (a+bu)(c+du) = (ac - bd) + (ad + bc)u   since u^2 = -1
        return Fp2(s.p, (s.a*o.a - s.b*o.b) % s.p, (s.a*o.b + s.b*o.a) % s.p)

    def sqr(s):
        # (a+bu)^2 = (a+b)(a-b) + 2ab*u
        return Fp2(s.p, ((s.a+s.b)*(s.a-s.b)) % s.p, (2*s.a*s.b) % s.p)

    def mul_int(s, k):  return Fp2(s.p, (s.a*k) % s.p, (s.b*k) % s.p)

    def inv(s):
        # conjugate / norm; norm = a^2 + b^2 lives in Fp
        n = (s.a*s.a + s.b*s.b) % s.p
        ninv = pow(n, s.p-2, s.p)
        return Fp2(s.p, (s.a*ninv) % s.p, (-s.b*ninv) % s.p)

    def mul_xi(s):
        # multiply by xi = 1 + u ; matches src/bn_fp2.c Fp2_mul_basis
        return Fp2(s.p, (s.a - s.b) % s.p, (s.a + s.b) % s.p)

    def is_zero(s): return s.a % s.p == 0 and s.b % s.p == 0
    def conj(s):    return Fp2(s.p, s.a % s.p, (-s.b) % s.p)
    def frob(s):    return s.conj()          # x -> x^p on Fp2 is conjugation


# --------------------------------------------------------------------- Fp6

@dataclass(frozen=True)
class Fp6:
    c0: Fp2
    c1: Fp2
    c2: Fp2

    @staticmethod
    def zero(p): return Fp6(Fp2.zero(p), Fp2.zero(p), Fp2.zero(p))
    @staticmethod
    def one(p):  return Fp6(Fp2.one(p),  Fp2.zero(p), Fp2.zero(p))

    @property
    def p(self): return self.c0.p

    def __add__(s, o): return Fp6(s.c0+o.c0, s.c1+o.c1, s.c2+o.c2)
    def __sub__(s, o): return Fp6(s.c0-o.c0, s.c1-o.c1, s.c2-o.c2)
    def __neg__(s):    return Fp6(-s.c0, -s.c1, -s.c2)

    def __mul__(s, o):
        # Karatsuba over Fp2 with v^3 = xi
        t0 = s.c0*o.c0
        t1 = s.c1*o.c1
        t2 = s.c2*o.c2
        e0 = ((s.c1+s.c2)*(o.c1+o.c2) - t1 - t2).mul_xi() + t0
        e1 = ((s.c0+s.c1)*(o.c0+o.c1) - t0 - t1) + t2.mul_xi()
        e2 = ((s.c0+s.c2)*(o.c0+o.c2) - t0 - t2) + t1
        return Fp6(e0, e1, e2)

    def sqr(s): return s * s

    def mul_v(s):
        # multiply by v ; matches src/bn_fp6.c Fp6_mul_basis
        return Fp6(s.c2.mul_xi(), s.c0, s.c1)

    def inv(s):
        # standard Fp6 inversion over Fp2
        A = s.c0.sqr() - (s.c1 * s.c2).mul_xi()
        B = s.c2.sqr().mul_xi() - (s.c0 * s.c1)
        C = s.c1.sqr() - (s.c0 * s.c2)
        F = ((s.c2 * B + s.c1 * C).mul_xi() + s.c0 * A).inv()
        return Fp6(A*F, B*F, C*F)

    def is_zero(s): return s.c0.is_zero() and s.c1.is_zero() and s.c2.is_zero()


# --------------------------------------------------------------------- Fp12

@dataclass(frozen=True)
class Fp12:
    d0: Fp6
    d1: Fp6

    @staticmethod
    def zero(p): return Fp12(Fp6.zero(p), Fp6.zero(p))
    @staticmethod
    def one(p):  return Fp12(Fp6.one(p),  Fp6.zero(p))

    @property
    def p(self): return self.d0.p

    def __add__(s, o): return Fp12(s.d0+o.d0, s.d1+o.d1)
    def __sub__(s, o): return Fp12(s.d0-o.d0, s.d1-o.d1)
    def __neg__(s):    return Fp12(-s.d0, -s.d1)

    def __mul__(s, o):
        # Karatsuba over Fp6 with w^2 = v
        t0 = s.d0 * o.d0
        t1 = s.d1 * o.d1
        return Fp12(t0 + t1.mul_v(),
                    (s.d0 + s.d1) * (o.d0 + o.d1) - t0 - t1)

    def sqr(s): return s * s

    def inv(s):
        # (d0 + d1 w)^-1 = (d0 - d1 w) / (d0^2 - v*d1^2)
        f = (s.d0.sqr() - s.d1.sqr().mul_v()).inv()
        return Fp12(s.d0 * f, (-s.d1) * f)

    def conj(s):  return Fp12(s.d0, -s.d1)      # the p^6 Frobenius

    def pow(s, e):
        if e < 0:
            return s.inv().pow(-e)
        r, base = Fp12.one(s.p), s
        while e:
            if e & 1:
                r = r * base
            base = base.sqr()
            e >>= 1
        return r

    def is_one(s):
        return (s.d0.c0.a % s.p == 1 and s.d0.c0.b % s.p == 0
                and s.d0.c1.is_zero() and s.d0.c2.is_zero() and s.d1.is_zero())

    def frob(s, k=1):
        """x -> x^(p^k), computed by direct exponentiation.

        Deliberately NOT using precomputed Frobenius constants: the point of an
        oracle is to avoid sharing machinery with the implementation under test.
        Slow, and that is acceptable here.
        """
        return s.pow(pow(s.p, k))


# ------------------------------------------------------- flat coordinate view
# The C library stores Fp12 as x0=(c0,c1,c2), x1=(c3,c4,c5) with each ci an Fp2
# pair. Flatten to 12 integers in exactly that order so vectors line up.

def fp12_to_list(f):
    out = []
    for six in (f.d0, f.d1):
        for c in (six.c0, six.c1, six.c2):
            out += [c.a % f.p, c.b % f.p]
    return out


def fp12_from_list(p, vals):
    assert len(vals) == 12
    def six(o):
        return Fp6(Fp2(p, vals[o], vals[o+1]),
                   Fp2(p, vals[o+2], vals[o+3]),
                   Fp2(p, vals[o+4], vals[o+5]))
    return Fp12(six(0), six(6))


# ------------------------------------------------------------ curve over Fp

class EFp:
    """Affine point on y^2 = x^3 + b_signed over Fp. None represents infinity."""

    def __init__(self, curve, x=None, y=None):
        self.cv = curve
        self.x = x
        self.y = y
        self.inf = x is None

    def is_on_curve(self):
        if self.inf:
            return True
        p = self.cv.p
        return (self.y*self.y - self.x*self.x*self.x - self.cv.b_signed) % p == 0

    def __neg__(self):
        return EFp(self.cv) if self.inf else EFp(self.cv, self.x, (-self.y) % self.cv.p)

    def dbl(self):
        p = self.cv.p
        if self.inf or self.y % p == 0:
            return EFp(self.cv)
        lam = (3*self.x*self.x) * pow(2*self.y, p-2, p) % p
        xr = (lam*lam - 2*self.x) % p
        return EFp(self.cv, xr, (lam*(self.x - xr) - self.y) % p)

    def __add__(self, o):
        p = self.cv.p
        if self.inf: return o
        if o.inf:    return self
        if (self.x - o.x) % p == 0:
            return self.dbl() if (self.y - o.y) % p == 0 else EFp(self.cv)
        lam = (o.y - self.y) * pow(o.x - self.x, p-2, p) % p
        xr = (lam*lam - self.x - o.x) % p
        return EFp(self.cv, xr, (lam*(self.x - xr) - self.y) % p)

    def mul(self, k):
        if k < 0:
            return (-self).mul(-k)
        r, base = EFp(self.cv), self
        while k:
            if k & 1:
                r = r + base
            base = base.dbl()
            k >>= 1
        return r

    def __eq__(self, o):
        if self.inf or o.inf:
            return self.inf and o.inf
        p = self.cv.p
        return (self.x - o.x) % p == 0 and (self.y - o.y) % p == 0


# ----------------------------------------------------------- curve over Fp2

class EFp2:
    """Affine point on the sextic twist y^2 = x^3 + b_signed*xi over Fp2.

    The twist equation follows from the library's own untwisting map in
    src/bn_twist.c: x_E12 = (x'/xi)*v^2 and y_E12 = (y'/xi)*v*w. With v = w^2
    and xi = w^6 that is x = x'*w^-2, y = y'*w^-3, and substituting into
    y^2 = x^3 + b gives y'^2 = x'^3 + b*xi. This is a D-type twist.
    """

    def __init__(self, curve, x=None, y=None):
        self.cv = curve
        self.x = x
        self.y = y
        self.inf = x is None

    @staticmethod
    def b_twist(curve):
        p = curve.p
        return Fp2(p, curve.b_signed, 0).mul_xi()

    def is_on_curve(self):
        if self.inf:
            return True
        return (self.y.sqr() - self.x.sqr()*self.x - EFp2.b_twist(self.cv)).is_zero()

    def __neg__(self):
        return EFp2(self.cv) if self.inf else EFp2(self.cv, self.x, -self.y)

    def dbl(self):
        if self.inf or self.y.is_zero():
            return EFp2(self.cv)
        lam = self.x.sqr().mul_int(3) * (self.y.mul_int(2)).inv()
        xr = lam.sqr() - self.x.mul_int(2)
        return EFp2(self.cv, xr, lam*(self.x - xr) - self.y)

    def __add__(self, o):
        if self.inf: return o
        if o.inf:    return self
        if (self.x - o.x).is_zero():
            return self.dbl() if (self.y - o.y).is_zero() else EFp2(self.cv)
        lam = (o.y - self.y) * (o.x - self.x).inv()
        xr = lam.sqr() - self.x - o.x
        return EFp2(self.cv, xr, lam*(self.x - xr) - self.y)

    def mul(self, k):
        if k < 0:
            return (-self).mul(-k)
        r, base = EFp2(self.cv), self
        while k:
            if k & 1:
                r = r + base
            base = base.dbl()
            k >>= 1
        return r

    def __eq__(self, o):
        if self.inf or o.inf:
            return self.inf and o.inf
        return (self.x - o.x).is_zero() and (self.y - o.y).is_zero()
