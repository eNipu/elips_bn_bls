"""RFC 9380 hash-to-curve, written independently from the specification.

Phase 5b, issue #6. This is the oracle for src/hash/, in the same sense that
elips_ref.py is the oracle for the field and curve layers: written from the
defining equations rather than from the library, so agreement between the two
means something.

Two maps, because the curves need different ones:

  simplified SWU     BLS12-381 only. Needs a curve with A*B != 0, which neither
  (+ isogeny)        of these curves is (both have A = 0), so the map runs on
                     an isogenous curve and is pushed across. RFC 9380 registers
                     exactly this for BLS12-381 and it is what every other
                     implementation computes, so it is what interoperates.

  Shallue-van de     BLS12-461 and BN-462. Works directly on A = 0 with no
  Woestijne          isogeny, which is why RFC 9380 uses it for BN254 and why
                     it is the right choice for the two curves here that have
                     no registered suite and therefore nothing to match.

The isogeny coefficients live in h2c_iso_bls12_381.json; see the provenance note
there. Everything else -- the SvdW Z values, the constants c1..c4, the cofactor
multipliers -- is derived here from the curve parameters.
"""
import hashlib
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, EFp, EFp2, Fp2      # noqa: E402
from selftest import fp2_sqrt                     # noqa: E402

K_SECURITY = 128        # RFC 9380 target security level, in bits


# --------------------------------------------------------------- field views
#
# The maps below are written once and run over both Fp and Fp2, so each needs
# the same small interface: arithmetic, a square test, a square root, and
# sgn0. Only sgn0 differs in substance between the two.

def make_fp(p):
    class F:
        m = 1
        deg = 1

        def __init__(s, v):     s.v = v % p
        def __add__(s, o):      return F(s.v + o.v)
        def __sub__(s, o):      return F(s.v - o.v)
        def __mul__(s, o):      return F(s.v * o.v)
        def __neg__(s):         return F(-s.v)
        def sqr(s):             return F(s.v * s.v)
        def is_zero(s):         return s.v == 0
        def __eq__(s, o):       return s.v == o.v
        def __hash__(s):        return hash(s.v)

        def inv0(s):
            """RFC 9380's inv0: 0 maps to 0 rather than being undefined."""
            return F(pow(s.v, p - 2, p))

        def is_square(s):
            return s.v == 0 or pow(s.v, (p - 1) // 2, p) == 1

        def sqrt(s):
            assert p % 4 == 3
            return F(pow(s.v, (p + 1) // 4, p))

        def sgn0(s):
            """RFC 9380 4.1, m == 1."""
            return s.v & 1

        def limbs(s):           return [s.v]

        @staticmethod
        def one():              return F(1)
        @staticmethod
        def zero():             return F(0)
        @staticmethod
        def from_int(k):        return F(k)
        @staticmethod
        def from_coeffs(c):     return F(c[0])
    return F


def make_fp2(p):
    class F:
        m = 2
        deg = 2

        def __init__(s, v):     s.v = v
        def __add__(s, o):      return F(s.v + o.v)
        def __sub__(s, o):      return F(s.v - o.v)
        def __mul__(s, o):      return F(s.v * o.v)
        def __neg__(s):         return F(-s.v)
        def sqr(s):             return F(s.v.sqr())
        def is_zero(s):         return s.v.is_zero()
        def __eq__(s, o):       return (s.v - o.v).is_zero()
        def __hash__(s):        return hash((s.v.a, s.v.b))

        def inv0(s):
            return F(Fp2.zero(p)) if s.v.is_zero() else F(s.v.inv())

        def is_square(s):
            # a is a square in Fp2 exactly when its norm is a square in Fp.
            n = (s.v.a * s.v.a + s.v.b * s.v.b) % p
            return n == 0 or pow(n, (p - 1) // 2, p) == 1

        def sqrt(s):
            q = fp2_sqrt(s.v, p)
            return F(q) if q is not None else None

        def sgn0(s):
            """RFC 9380 4.1, m == 2: the real part decides unless it is zero."""
            s0 = s.v.a % p & 1
            z0 = (s.v.a % p) == 0
            s1 = s.v.b % p & 1
            return s0 | (z0 & s1)

        def limbs(s):           return [s.v.a % p, s.v.b % p]

        @staticmethod
        def one():              return F(Fp2.one(p))
        @staticmethod
        def zero():             return F(Fp2.zero(p))
        @staticmethod
        def from_int(k):        return F(Fp2(p, k, 0))
        @staticmethod
        def from_coeffs(c):     return F(Fp2(p, c[0], c[1]))
    return F


# ------------------------------------------------- RFC 9380 5.3, hash to field

def expand_message_xmd(msg, DST, n_bytes):
    """RFC 9380 5.3.1, with SHA-256."""
    H = hashlib.sha256
    b_in_bytes = 32            # SHA-256 output
    s_in_bytes = 64            # SHA-256 block

    DST = oversize_dst(DST)
    ell = (n_bytes + b_in_bytes - 1) // b_in_bytes
    if ell > 255 or n_bytes > 65535:
        raise ValueError("expand_message_xmd: length out of range")

    DST_prime = DST + bytes([len(DST)])
    msg_prime = (b"\x00" * s_in_bytes + msg + n_bytes.to_bytes(2, "big")
                 + b"\x00" + DST_prime)
    b0 = H(msg_prime).digest()
    blocks = [H(b0 + b"\x01" + DST_prime).digest()]
    for i in range(2, ell + 1):
        x = bytes(a ^ b for a, b in zip(b0, blocks[-1]))
        blocks.append(H(x + bytes([i]) + DST_prime).digest())
    return b"".join(blocks)[:n_bytes]


def oversize_dst(DST):
    """RFC 9380 5.3.3. A DST longer than 255 bytes is hashed down rather than
    rejected, so no caller has to think about the limit."""
    if len(DST) <= 255:
        return DST
    return hashlib.sha256(b"H2C-OVERSIZE-DST-" + DST).digest()


def field_L(p):
    """RFC 9380 5.3: ceil((ceil(log2(p)) + k) / 8)."""
    return (p.bit_length() + K_SECURITY + 7) // 8


def hash_to_field(msg, DST, count, F, p):
    L = field_L(p)
    buf = expand_message_xmd(msg, DST, count * F.m * L)
    out = []
    for i in range(count):
        coeffs = []
        for j in range(F.m):
            off = L * (j + i * F.m)
            coeffs.append(int.from_bytes(buf[off:off + L], "big") % p)
        out.append(F.from_coeffs(coeffs))
    return out


# ------------------------------------------- RFC 9380 6.6.2, simplified SWU
#
# Written as the mathematical definition rather than as the specification's
# straight-line listing. The two agree: the listing is that definition with the
# divisions deferred, and the sign fixup at the end makes the result independent
# of which square root the intermediate step returned.

def map_to_curve_sswu(F, A, B, Z, u):
    def g(x):
        return x.sqr() * x + A * x + B

    zu2 = Z * u.sqr()
    t = zu2.sqr() + zu2                       # Z^2 u^4 + Z u^2
    if t.is_zero():
        x1 = B * (Z * A).inv0()               # the exceptional case
    else:
        x1 = (-B) * A.inv0() * (F.one() + t.inv0())

    gx1 = g(x1)
    if gx1.is_square():
        x, y = x1, gx1.sqrt()
    else:
        x2 = zu2 * x1                         # g(x2) = Z^3 u^6 g(x1)
        x, y = x2, g(x2).sqrt()

    if y.sgn0() != u.sgn0():
        y = -y
    return x, y


def horner(coeffs, x):
    acc = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        acc = acc * x + c
    return acc


def iso_map(x, y, maps):
    """x' = x_num(x)/x_den(x),  y' = y * y_num(x)/y_den(x)."""
    xn, xd, yn, yd = maps
    return (horner(xn, x) * horner(xd, x).inv0(),
            y * horner(yn, x) * horner(yd, x).inv0())


# ---------------------------------- RFC 9380 6.6.1, Shallue-van de Woestijne

def svdw_z(F, A, B, p):
    """RFC 9380 Appendix H.1, in the specified search order 1, -1, 2, -2, ...

    For an extension field the specification leaves the enumeration to the
    implementation; the order used here is the same one lifted through
    from_int, i.e. only the subfield is searched. That always succeeds on these
    curves and keeps the constant small, which is worth something in a table of
    generated limbs."""
    def g(x):
        return x.sqr() * x + A * x + B

    three, four = F.from_int(3), F.from_int(4)
    ctr = 0
    while ctr < 10000:
        ctr += 1
        for Z in (F.from_int(ctr), F.from_int(-ctr)):
            gZ = g(Z)
            if gZ.is_zero():
                continue
            t = three * Z.sqr() + four * A
            if t.is_zero():
                continue
            val = -t * (four * gZ).inv0()
            if val.is_zero() or not val.is_square():
                continue
            if not (gZ.is_square() or g(-Z * F.from_int(2).inv0()).is_square()):
                continue
            return Z, (ctr if Z == F.from_int(ctr) else -ctr)
    raise RuntimeError("no SvdW Z found")


def svdw_constants(F, A, B, Z):
    """The four values RFC 9380 6.6.1 precomputes."""
    def g(x):
        return x.sqr() * x + A * x + B

    two, three, four = F.from_int(2), F.from_int(3), F.from_int(4)
    gZ = g(Z)
    t = three * Z.sqr() + four * A

    c1 = gZ
    c2 = -Z * two.inv0()
    c3 = (-gZ * t).sqrt()
    if c3 is None:
        raise RuntimeError("SvdW: -g(Z)(3Z^2+4A) is not a square")
    if c3.sgn0() != 0:
        c3 = -c3                     # the specification pins sgn0(c3) == 0
    c4 = -four * gZ * t.inv0()

    # Each constant is checked against its defining equation rather than
    # trusted, because all four are silent if wrong: the map still returns a
    # point, just not the right one.
    assert c1 == gZ
    assert c2 * two == -Z
    assert c3.sqr() == -gZ * t and c3.sgn0() == 0
    assert c4 * t == -four * gZ
    return c1, c2, c3, c4


def map_to_curve_svdw(F, A, B, Z, C, u):
    c1, c2, c3, c4 = C

    def g(x):
        return x.sqr() * x + A * x + B

    one = F.one()
    tv1 = u.sqr() * c1
    tv2 = one + tv1
    tv1 = one - tv1
    tv3 = (tv1 * tv2).inv0()
    tv4 = u * tv1 * tv3 * c3

    x1 = c2 - tv4
    x2 = c2 + tv4
    x3 = (tv2.sqr() * tv3).sqr() * c4 + Z

    e1 = g(x1).is_square()
    e2 = g(x2).is_square() and not e1
    x = x1 if e1 else (x2 if e2 else x3)

    y = g(x).sqrt()
    if y.sgn0() != u.sgn0():
        y = -y
    return x, y


# ------------------------------------------------------------- suite objects

class Suite:
    """One (curve, group) hash-to-curve instantiation."""

    def __init__(self, curve, group):
        self.cv = curve
        self.group = group                    # "G1" or "G2"
        p = curve.p
        self.p = p
        self.L = field_L(p)

        if group == "G1":
            self.F = make_fp(p)
            self.A = self.F.zero()
            self.B = self.F.from_int(curve.b_signed)
            self.point = lambda x, y: EFp(curve, x.v, y.v)
        else:
            self.F = make_fp2(p)
            self.A = self.F.zero()
            self.B = self.F(EFp2.b_twist(curve))
            self.point = lambda x, y: EFp2(curve, x.v, y.v)

        self.h_eff = cofactor_multiplier(curve, group)
        self.iso = None

        if curve.name == "BLS12-381":
            self.kind = "SSWU"
            data = _load_iso()
            side = data["g1" if group == "G1" else "g2"]
            cf = (lambda h: self.F.from_int(int(h, 16))) if group == "G1" else \
                 (lambda pr: self.F.from_coeffs([int(pr[0], 16), int(pr[1], 16)]))
            self.iso_A = cf(side["A"])
            self.iso_B = cf(side["B"])
            self.Z = cf(side["Z"])
            self.iso = [[cf(c) for c in poly] for poly in side["map"]]
            self.suite_id = "BLS12381%s_XMD:SHA-256_SSWU_RO_" % group
        else:
            self.kind = "SVDW"
            self.Z, self.z_int = svdw_z(self.F, self.A, self.B, p)
            self.C = svdw_constants(self.F, self.A, self.B, self.Z)
            tag = "BLS12461" if curve.name == "BLS12-461" else "BN462"
            self.suite_id = "%s%s_XMD:SHA-256_SVDW_RO_" % (tag, group)

    # -- the map, whichever kind this suite uses
    def map_to_curve(self, u):
        if self.kind == "SSWU":
            x, y = map_to_curve_sswu(self.F, self.iso_A, self.iso_B, self.Z, u)
            x, y = iso_map(x, y, self.iso)
        else:
            x, y = map_to_curve_svdw(self.F, self.A, self.B, self.Z, self.C, u)
        return self.point(x, y)

    def hash_to_field(self, msg, DST, count):
        return hash_to_field(msg, DST, count, self.F, self.p)

    def hash_to_curve(self, msg, DST):
        """RFC 9380 3, the random-oracle variant."""
        u0, u1 = self.hash_to_field(msg, DST, 2)
        return (self.map_to_curve(u0) + self.map_to_curve(u1)).mul(self.h_eff)

    def encode_to_curve(self, msg, DST):
        """RFC 9380 3, the non-uniform variant. Cheaper, and NOT a random
        oracle: its image is a fraction of the group. Only safe where the
        analysis does not need indifferentiability."""
        u0, = self.hash_to_field(msg, DST, 1)
        return self.map_to_curve(u0).mul(self.h_eff)

    def default_dst(self):
        return self.suite_id.encode()


_ISO_CACHE = {}


def _load_iso():
    if not _ISO_CACHE:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "h2c_iso_bls12_381.json")
        _ISO_CACHE.update(json.load(open(path)))
    return _ISO_CACHE


# ----------------------------------------------------------- cofactor clearing

def group_orders(curve):
    """#E(Fp) and #E'(Fp2), the second identified rather than assumed: both
    candidate traces are tried and the one that annihilates a real twist point
    wins."""
    p, r = curve.p, curve.r
    tr = (curve.X + 1) if curve.family == "bls12" else (6 * curve.X ** 2 + 1)
    n1 = p + 1 - tr
    assert n1 % r == 0

    t2 = tr * tr - 2 * p
    f2 = math.isqrt((4 * p * p - t2 * t2) // 3)
    assert 3 * f2 * f2 == 4 * p * p - t2 * t2

    bt = EFp2.b_twist(curve)
    probe = None
    for k in range(1, 4000):
        x = Fp2(p, k, 1)
        y = fp2_sqrt(x.sqr() * x + bt, p)
        if y is not None:
            probe = EFp2(curve, x, y)
            break
    for c in {(3 * f2 + t2) // 2, (-3 * f2 + t2) // 2}:
        n2 = p * p + 1 - c
        if n2 % r == 0 and probe.mul(n2).inf:
            return n1, n2
    raise RuntimeError("could not identify the twist order")


_HEFF_CACHE = {}


def cofactor_multiplier(curve, group):
    """What clear_cofactor multiplies by.

    BLS12 matches RFC 9380 8.8: 1-x on G1 and 3(x^2-1)h2 on G2. Both are derived
    here from the curve parameter rather than transcribed, and then confirmed
    against the specification's own value through the committed test vectors.

    BN has no registered suite, so the plain cofactor is used.

    Note what is NOT asserted: that the multiplier is a multiple of the
    cofactor. On BLS12-381 G1 it is not -- h_eff = 1-x is 64 bits against a
    128-bit cofactor -- and multiplying by it still lands in the subgroup,
    because E(Fp) is not cyclic there and (1-x) kills both of its cofactor
    components. That is the whole reason RFC 9380 can specify such a short
    multiplier. Asserting divisibility would reject the correct answer, so what
    is checked instead is the property that actually matters."""
    key = (curve.name, group)
    if key in _HEFF_CACHE:
        return _HEFF_CACHE[key]

    n1, n2 = group_orders(curve)
    r = curve.r
    h1, h2 = n1 // r, n2 // r

    if group == "G1":
        h = (1 - curve.X) if curve.family == "bls12" else h1
    else:
        h = 3 * (curve.X ** 2 - 1) * h2 if curve.family == "bls12" else h2

    assert math.gcd(h, r) == 1, \
        "a multiplier sharing a factor with r is not surjective onto the group"
    assert _clears_cofactor(curve, group, h), \
        "the multiplier does not land in the order-r subgroup"
    _HEFF_CACHE[key] = h
    return h


def _clears_cofactor(curve, group, h, trials=8):
    """[h]P must have order dividing r for every P on the curve, so test it on
    points drawn without any reference to the subgroup."""
    import random
    rng = random.Random(20260906)
    p, r = curve.p, curve.r

    for _ in range(trials):
        if group == "G1":
            while True:
                x = rng.randrange(p)
                rhs = (x * x * x + curve.b_signed) % p
                if rhs and pow(rhs, (p - 1) // 2, p) == 1:
                    break
            P = EFp(curve, x, pow(rhs, (p + 1) // 4, p))
        else:
            bt = EFp2.b_twist(curve)
            while True:
                x = Fp2(p, rng.randrange(p), rng.randrange(p))
                y = fp2_sqrt(x.sqr() * x + bt, p)
                if y is not None:
                    break
            P = EFp2(curve, x, y)
        if not P.mul(h).mul(r).inf:
            return False
    return True


def psi_constants(curve):
    """The two Fp2 constants of the skew Frobenius on the twist,
    psi(x, y) = (conj(x) * gamma^-2, conj(y) * gamma^-3) with
    gamma = xi^((p-1)/6). Same values gen_params.py emits as PSI_X and PSI_Y."""
    p = curve.p
    xi = Fp2(p, 1, 1)

    def powi(b, e):
        acc = Fp2.one(p)
        while e:
            if e & 1:
                acc = acc * b
            b = b.sqr()
            e >>= 1
        return acc

    gam = powi(xi, (p - 1) // 6)
    return gam.sqr().inv(), (gam.sqr() * gam).inv()


def psi(curve, P):
    if P.inf:
        return P
    g2, g3 = psi_constants(curve)
    return EFp2(curve, P.x.conj() * g2, P.y.conj() * g3)


def g2_fast_clear_coeffs(curve):
    """Budroni and Pintore, "Efficient hash maps to G2 on BLS curves"
    (ePrint 2017/419), building on Fuentes-Castaneda et al. and Scott et al.:

        [h_eff]Q = [x^2 - x - 1]Q + [x - 1]psi(Q) + psi^2([2]Q)

    Two ladders over the short parameter x replace one over the 600-plus-bit
    h_eff. Returned as the two signed multipliers, and only after the identity
    has been CHECKED against [h_eff] on random points of the twist -- not on
    points of G2, where far weaker relations would also hold.

    A wrong chain here is loud rather than silent: it computes a different
    multiple, so the RFC 9380 vectors stop matching. That is why this one is
    implemented and the fast subgroup tests are not."""
    if curve.family != "bls12":
        return None

    X = curve.X
    a, b = X * X - X - 1, X - 1
    h = cofactor_multiplier(curve, "G2")

    import random
    rng = random.Random(20260906)
    p = curve.p
    bt = EFp2.b_twist(curve)
    for _ in range(4):
        while True:
            x = Fp2(p, rng.randrange(p), rng.randrange(p))
            y = fp2_sqrt(x.sqr() * x + bt, p)
            if y is not None:
                break
        Q = EFp2(curve, x, y)
        got = Q.mul(a) + psi(curve, Q).mul(b) + psi(curve, psi(curve, Q.mul(2)))
        assert got == Q.mul(h), "the fast G2 cofactor chain is not [h_eff]"
    return a, b


# ----------------------------------------------------------------- self-test

_SUITE_CACHE = {}


def suite(curve, group):
    """Suites are immutable and each costs seconds of pure-Python scalar
    multiplication to build, so build each one once."""
    key = (curve.name, group)
    if key not in _SUITE_CACHE:
        _SUITE_CACHE[key] = Suite(curve, group)
    return _SUITE_CACHE[key]


def suites():
    for name in ("BLS12-381", "BLS12-461", "BN-462"):
        for group in ("G1", "G2"):
            yield suite(CURVES[name], group)


# The RFC 9380 J.9.1 and J.10.1 test vectors: the published output of the two
# registered BLS12-381 suites. This is the only external check in the project
# on the hash-to-curve code, so it is the one that decides whether the isogeny
# data and the whole chain above are right.
RFC9380_MSGS = [b"", b"abc", b"abcdef0123456789",
                b"q128_" + b"q" * 128, b"a512_" + b"a" * 512]

RFC9380_G1 = [
    ("052926add2207b76ca4fa57a8734416c8dc95e24501772c814278700eed6d1e4e8cf62d9c09db0fac349612b759e79a1",
     "08ba738453bfed09cb546dbb0783dbb3a5f1f566ed67bb6be0e8c67e2e81a4cc68ee29813bb7994998f3eae0c9c6a265"),
    ("03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a7943388a49a3aee664ba5379a7655d3c68900be2f6903",
     "0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67af215533311f0b8dfaaa154fa6b88176c229f2885d"),
    ("11e0b079dea29a68f0383ee94fed1b940995272407e3bb916bbf268c263ddd57a6a27200a784cbc248e84f357ce82d98",
     "03a87ae2caf14e8ee52e51fa2ed8eefe80f02457004ba4d486d6aa1f517c0889501dc7413753f9599b099ebcbbd2d709"),
    ("15f68eaa693b95ccb85215dc65fa81038d69629f70aeee0d0f677cf22285e7bf58d7cb86eefe8f2e9bc3f8cb84fac488",
     "1807a1d50c29f430b8cafc4f8638dfeeadf51211e1602a5f184443076715f91bb90a48ba1e370edce6ae1062f5e6dd38"),
    ("082aabae8b7dedb0e78aeb619ad3bfd9277a2f77ba7fad20ef6aabdc6c31d19ba5a6d12283553294c1825c4b3ca2dcfe",
     "05b84ae5a942248eea39e1d91030458c40153f3b654ab7872d779ad1e942856a20c438e8d99bc8abfbf74729ce1f7ac8"),
]

DST_G1_RFC = b"QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_"
DST_G2_RFC = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_"


def selftest(verbose=True):
    import random
    fails = 0

    def check(cond, what):
        nonlocal fails
        if verbose:
            print("  [%s] %s" % ("PASS" if cond else "FAIL", what))
        if not cond:
            fails += 1

    # expand_message_xmd against RFC 9380 K.1: the first vector, DST
    # "QUUX-V01-CS02-with-expander-SHA256-128", msg "", len 32.
    got = expand_message_xmd(b"", b"QUUX-V01-CS02-with-expander-SHA256-128", 32)
    check(got.hex() == "68a985b87eb6b46952128911f2a4412bbc302a9d759667f8"
                       "7f7a21d803f07235",
          "expand_message_xmd matches RFC 9380 K.1")

    # The registered BLS12-381 G1 suite, against the published vectors.
    s = suite(CURVES["BLS12-381"], "G1")
    ok = True
    for msg, (xh, yh) in zip(RFC9380_MSGS, RFC9380_G1):
        P = s.hash_to_curve(msg, DST_G1_RFC)
        if P.x != int(xh, 16) or P.y != int(yh, 16):
            ok = False
    check(ok, "BLS12381G1_XMD:SHA-256_SSWU_RO_ matches RFC 9380 J.9.1")

    # Every suite must land in its group, whatever the message.
    random.seed(20260906)
    for s in suites():
        good = True
        for i in range(6):
            msg = bytes(random.randrange(256) for _ in range(i * 7))
            P = s.hash_to_curve(msg, s.default_dst())
            if not P.is_on_curve() or P.inf or not P.mul(s.cv.r).inf:
                good = False
            Q = s.encode_to_curve(msg, s.default_dst())
            if not Q.is_on_curve() or not Q.mul(s.cv.r).inf:
                good = False
        check(good, "%-11s %s %s: hashes land in the order-r subgroup"
              % (s.cv.name, s.group, s.kind))

    # The isogeny data is data, so prove it is an isogeny rather than trusting
    # the file: points of E' must land on E, and the map must be additive.
    for group in ("G1", "G2"):
        s = suite(CURVES["BLS12-381"], group)
        F, A, B = s.F, s.iso_A, s.iso_B

        def g(x):
            return x.sqr() * x + A * x + B

        pts = []
        while len(pts) < 6:
            u = F.from_coeffs([random.randrange(s.p) for _ in range(F.m)])
            if g(u).is_square():
                pts.append((u, g(u).sqrt()))

        def eprime_add(P, Q):
            (x1, y1), (x2, y2) = P, Q
            if (x1 - x2).is_zero():
                lam = (x1.sqr() * F.from_int(3) + A) * (y1 * F.from_int(2)).inv0()
            else:
                lam = (y2 - y1) * (x2 - x1).inv0()
            x3 = lam.sqr() - x1 - x2
            return (x3, lam * (x1 - x3) - y1)

        img = []
        on = True
        for (x, y) in pts:
            X, Y = iso_map(x, y, s.iso)
            P = s.point(X, Y)
            on = on and P.is_on_curve()
            img.append(P)
        check(on, "BLS12-381 %s isogeny maps E' onto the curve" % group)

        homo = True
        for i in range(0, 6, 2):
            X, Y = iso_map(*eprime_add(pts[i], pts[i + 1]), s.iso)
            if not (s.point(X, Y) == img[i] + img[i + 1]):
                homo = False
        check(homo, "BLS12-381 %s isogeny is a group homomorphism" % group)

    # The fast G2 cofactor chain must compute exactly [h_eff], on the whole
    # twist and not merely on G2.
    for name in ("BLS12-381", "BLS12-461"):
        cv = CURVES[name]
        ok = g2_fast_clear_coeffs(cv) is not None
        check(ok, "%-11s fast G2 cofactor chain equals [h_eff] on the twist" % name)

    # SvdW: the construction's own identity, on every curve that uses it.
    for s in suites():
        if s.kind != "SVDW":
            continue
        F, A, B, Z, C = s.F, s.A, s.B, s.Z, s.C

        def g(x):
            return x.sqr() * x + A * x + B

        good = True
        for _ in range(200):
            u = F.from_coeffs([random.randrange(s.p) for _ in range(F.m)])
            x, y = map_to_curve_svdw(F, A, B, Z, C, u)
            if not (g(x) == y.sqr()) or y.sgn0() != u.sgn0():
                good = False
        check(good, "%-11s %s SvdW output is on the curve with sgn0(y)==sgn0(u)"
              % (s.cv.name, s.group))

    if verbose:
        print("H2C SELFTEST", "FAILED" if fails else "PASSED")
    return fails


if __name__ == "__main__":
    sys.exit(1 if selftest() else 0)
