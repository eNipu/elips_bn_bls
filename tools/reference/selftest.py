"""Field-axiom self-test for the reference implementation (issue #9).

An oracle that is itself wrong is worse than no oracle, so this runs before any
vector is generated.
"""
import random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from elips_ref import CURVES, Fp2, Fp6, Fp12, EFp, EFp2, fp12_to_list, fp12_from_list

FAIL = 0

def fp2_sqrt(a, p):
    """Square root in Fp2 for p = 3 mod 4, by the complex method."""
    if a.is_zero():
        return Fp2.zero(p)
    nrm = (a.a*a.a + a.b*a.b) % p
    if pow(nrm, (p-1)//2, p) != 1:
        return None
    s = pow(nrm, (p+1)//4, p)               # sqrt of the norm in Fp
    for cand in (s, (-s) % p):
        t = (a.a + cand) * pow(2, p-2, p) % p
        if pow(t, (p-1)//2, p) == 1:
            x0 = pow(t, (p+1)//4, p)
            x1 = a.b * pow(2*x0, p-2, p) % p if x0 else 0
            r = Fp2(p, x0, x1)
            if r.sqr() == a:
                return r
    return None


def check(cond, msg):
    global FAIL
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)

def rnd_fp2(p, rng):  return Fp2(p, rng.randrange(p), rng.randrange(p))
def rnd_fp6(p, rng):  return Fp6(rnd_fp2(p, rng), rnd_fp2(p, rng), rnd_fp2(p, rng))
def rnd_fp12(p, rng): return Fp12(rnd_fp6(p, rng), rnd_fp6(p, rng))

def axioms(name, rnd, one, p, rng, n=40):
    for _ in range(n):
        a, b, c = rnd(p, rng), rnd(p, rng), rnd(p, rng)
        check(a*b == b*a,                     f"{name} mul commutative")
        check((a*b)*c == a*(b*c),             f"{name} mul associative")
        check(a*(b+c) == a*b + a*c,           f"{name} distributive")
        check(a*one == a,                     f"{name} mul identity")
        check(a.sqr() == a*a,                 f"{name} sqr == mul self")
        check(a + (-a) == a - a,              f"{name} negation")
        if not (hasattr(a, 'is_zero') and a.is_zero()):
            check(a * a.inv() == one,         f"{name} inverse")

for cname, cv in CURVES.items():
    print(f"\n=== {cname} (p = {cv.p.bit_length()} bits) ===")
    rng = random.Random(0xE11B5 ^ cv.p)
    p = cv.p

    axioms("Fp2",  rnd_fp2,  Fp2.one(p),  p, rng)
    axioms("Fp6",  rnd_fp6,  Fp6.one(p),  p, rng)
    axioms("Fp12", rnd_fp12, Fp12.one(p), p, rng, n=12)

    # mul_xi must equal multiplication by the element 1+u
    xi = Fp2(p, 1, 1)
    for _ in range(20):
        a = rnd_fp2(p, rng)
        check(a.mul_xi() == a * xi, "Fp2.mul_xi == *(1+u)")

    # mul_v must equal multiplication by v = (0,1,0)
    v = Fp6(Fp2.zero(p), Fp2.one(p), Fp2.zero(p))
    for _ in range(20):
        a = rnd_fp6(p, rng)
        check(a.mul_v() == a * v, "Fp6.mul_v == *v")

    # v^3 must equal xi, and w^2 must equal v
    check(v*v*v == Fp6(xi, Fp2.zero(p), Fp2.zero(p)), "v^3 == xi")
    w = Fp12(Fp6.zero(p), Fp6.one(p))
    check(w.sqr() == Fp12(v, Fp6.zero(p)), "w^2 == v")

    # Frobenius: x^p is a field homomorphism, and x^(p^12) == x
    for _ in range(3):
        a, b = rnd_fp12(p, rng), rnd_fp12(p, rng)
        check(a.frob(1) * b.frob(1) == (a*b).frob(1), "Frobenius multiplicative")
        check(a.frob(12) == a,                        "x^(p^12) == x")
    # conj is the p^6 Frobenius: verify on the cyclotomic subgroup, where it holds
    for _ in range(2):
        a = rnd_fp12(p, rng)
        cyc = a.frob(6) * a.inv()          # norm-1 element
        check(cyc.conj() == cyc.inv(),     "conj == inverse on cyclotomic subgroup")

    # flat coordinate round-trip must be exact
    for _ in range(10):
        a = rnd_fp12(p, rng)
        check(fp12_from_list(p, fp12_to_list(a)) == a, "fp12 flat round-trip")

    # curve over Fp
    G = None
    for xi_ in range(1, 400):
        rhs = (xi_**3 + cv.b_signed) % p
        if pow(rhs, (p-1)//2, p) == 1:
            G = EFp(cv, xi_, pow(rhs, (p+1)//4, p)); break
    check(G is not None and G.is_on_curve(), "found a point on E(Fp)")
    if G:
        check((G + (-G)).inf,            "P + (-P) == O")
        check(G.dbl() == G + G,          "dbl == add self")
        check(G.mul(5) == G+G+G+G+G,     "SCM 5 == repeated add")
        check(G.mul(0).inf,              "[0]P == O")
        # cofactor-clear into the r-torsion, then check order
        h = (p + 1 - (cv.X + 1 if cv.family == "bls12" else 6*cv.X**2 + 1)) // cv.r
        P = G.mul(h)
        if not P.inf:
            check(P.mul(cv.r).inf,       "[r]P == O for P in G1")

    # curve over Fp2 (the sextic twist)
    bt = EFp2.b_twist(cv)
    Q = None
    for k in range(1, 400):
        x = Fp2(p, k, 1)
        rhs = x.sqr()*x + bt
        # try a square root of the Fp2 element via norm test
        nrm = (rhs.a*rhs.a + rhs.b*rhs.b) % p
        if pow(nrm, (p-1)//2, p) != 1:
            continue
        y = fp2_sqrt(rhs, p)
        if y is not None:
            Q = EFp2(cv, x, y); break
    check(Q is not None and Q.is_on_curve(), "found a point on the twist E'(Fp2)")
    if Q:
        check((Q + (-Q)).inf,        "Q + (-Q) == O")
        check(Q.dbl() == Q + Q,      "twist dbl == add self")
        check(Q.mul(7) == Q+Q+Q+Q+Q+Q+Q, "twist SCM 7 == repeated add")

print()
if FAIL:
    print(f"SELFTEST FAILED: {FAIL} check(s)")
    sys.exit(1)
print("SELFTEST PASSED")
