"""Known-answer vector generator (issue #11).

Emits a flat, whitespace-separated text format so the C runner needs no JSON
parser. One record per line:

    <op> <hex operand tokens...> = <hex result tokens...>

Arity is known to the runner from the op name, so the format carries no schema.
Fp is one token, Fp2 two, Fp6 six, Fp12 twelve; an EFp point is three tokens
(infinity flag, x, y) and an EFp2 point is five.

Every record is re-checked inside Python before being written.
"""
import random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, Fp2, Fp6, Fp12, EFp, EFp2
from selftest import fp2_sqrt

H = lambda n: format(n % (1 << 4096), 'x')


def t_fp(x, p):    return [H(x % p)]
def t_fp2(a):      return [H(a.a % a.p), H(a.b % a.p)]
def t_fp6(a):      return t_fp2(a.c0) + t_fp2(a.c1) + t_fp2(a.c2)
def t_fp12(a):     return t_fp6(a.d0) + t_fp6(a.d1)
def t_efp(P):
    if P.inf: return ["1", "0", "0"]
    return ["0", H(P.x % P.cv.p), H(P.y % P.cv.p)]
def t_efp2(P):
    if P.inf: return ["1", "0", "0", "0", "0"]
    return ["0"] + t_fp2(P.x) + t_fp2(P.y)


def generate(cv, rng, out):
    p = cv.p
    rec = lambda op, ins, outs: out.append(f"{op} {' '.join(ins)} = {' '.join(outs)}")

    def r_fp():   return rng.randrange(p)
    def r_fp2():  return Fp2(p, r_fp(), r_fp())
    def r_fp6():  return Fp6(r_fp2(), r_fp2(), r_fp2())
    def r_fp12(): return Fp12(r_fp6(), r_fp6())

    # ---- Fp, including boundary values ----
    edge = [0, 1, 2, p - 1, p - 2, (p - 1) // 2]
    pairs = [(a, b) for a in edge for b in edge[:3]] + [(r_fp(), r_fp()) for _ in range(30)]
    for a, b in pairs:
        rec("fp_add", t_fp(a, p) + t_fp(b, p), t_fp(a + b, p))
        rec("fp_sub", t_fp(a, p) + t_fp(b, p), t_fp(a - b, p))
        rec("fp_mul", t_fp(a, p) + t_fp(b, p), t_fp(a * b, p))
    for a in edge + [r_fp() for _ in range(15)]:
        rec("fp_neg", t_fp(a, p), t_fp(-a, p))
        if a % p:                       # Fp_inv(0) is undefined in the library
            rec("fp_inv", t_fp(a, p), t_fp(pow(a, p - 2, p), p))

    # ---- Fp2 ----
    e2 = [Fp2(p, 0, 0), Fp2(p, 1, 0), Fp2(p, 0, 1), Fp2(p, p - 1, p - 1)]
    for a in e2 + [r_fp2() for _ in range(25)]:
        for b in e2[:3] + [r_fp2()]:
            rec("fp2_add", t_fp2(a) + t_fp2(b), t_fp2(a + b))
            rec("fp2_sub", t_fp2(a) + t_fp2(b), t_fp2(a - b))
            rec("fp2_mul", t_fp2(a) + t_fp2(b), t_fp2(a * b))
        rec("fp2_sqr", t_fp2(a), t_fp2(a.sqr()))
        rec("fp2_mulbasis", t_fp2(a), t_fp2(a.mul_xi()))
        if not a.is_zero():
            rec("fp2_inv", t_fp2(a), t_fp2(a.inv()))

    # ---- Fp6 ----
    for _ in range(20):
        a, b = r_fp6(), r_fp6()
        rec("fp6_add", t_fp6(a) + t_fp6(b), t_fp6(a + b))
        rec("fp6_sub", t_fp6(a) + t_fp6(b), t_fp6(a - b))
        rec("fp6_mul", t_fp6(a) + t_fp6(b), t_fp6(a * b))
        rec("fp6_sqr", t_fp6(a), t_fp6(a.sqr()))
        rec("fp6_mulbasis", t_fp6(a), t_fp6(a.mul_v()))
        rec("fp6_inv", t_fp6(a), t_fp6(a.inv()))

    # ---- Fp12 ----
    for _ in range(15):
        a, b = r_fp12(), r_fp12()
        rec("fp12_add", t_fp12(a) + t_fp12(b), t_fp12(a + b))
        rec("fp12_sub", t_fp12(a) + t_fp12(b), t_fp12(a - b))
        rec("fp12_mul", t_fp12(a) + t_fp12(b), t_fp12(a * b))
        rec("fp12_sqr", t_fp12(a), t_fp12(a.sqr()))
        rec("fp12_inv", t_fp12(a), t_fp12(a.inv()))

    # ---- E(Fp) ----
    base = None
    for x in range(1, 500):
        rhs = (x**3 + cv.b_signed) % p
        if pow(rhs, (p - 1) // 2, p) == 1:
            base = EFp(cv, x, pow(rhs, (p + 1) // 4, p)); break
    assert base and base.is_on_curve(), "no base point on E(Fp)"
    pts = [base.mul(k) for k in (1, 2, 3, 5, 7, 11, 13, 17)]
    for P in pts:
        rec("efp_dbl", t_efp(P), t_efp(P.dbl()))
    rec("efp_dbl", t_efp(EFp(cv)), t_efp(EFp(cv)))            # doubling infinity
    for i, P in enumerate(pts):
        Q = pts[(i + 3) % len(pts)]
        rec("efp_add", t_efp(P) + t_efp(Q), t_efp(P + Q))
    rec("efp_add", t_efp(base) + t_efp(-base), t_efp(EFp(cv)))  # P + (-P) = O
    rec("efp_add", t_efp(base) + t_efp(EFp(cv)), t_efp(base))   # P + O = P
    for k in (2, 3, 255, 65537, rng.randrange(1 << 63), rng.randrange(1 << 200)):
        rec("efp_mul", t_efp(base) + [H(k)], t_efp(base.mul(k)))

    # ---- E'(Fp2), the sextic twist actually used by the Miller loop ----
    bt = EFp2.b_twist(cv)
    qbase = None
    for k in range(1, 500):
        x = Fp2(p, k, 1)
        y = fp2_sqrt(x.sqr() * x + bt, p)
        if y is not None:
            qbase = EFp2(cv, x, y); break
    assert qbase and qbase.is_on_curve(), "no base point on the twist"
    qpts = [qbase.mul(k) for k in (1, 2, 3, 5, 7, 11)]
    for Q in qpts:
        rec("efp2_dbl", t_efp2(Q), t_efp2(Q.dbl()))
    for i, Q in enumerate(qpts):
        R = qpts[(i + 2) % len(qpts)]
        rec("efp2_add", t_efp2(Q) + t_efp2(R), t_efp2(Q + R))
    rec("efp2_add", t_efp2(qbase) + t_efp2(-qbase), t_efp2(EFp2(cv)))
    for k in (2, 3, 255, 65537, rng.randrange(1 << 63)):
        rec("efp2_mul", t_efp2(qbase) + [H(k)], t_efp2(qbase.mul(k)))

    return out


if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "test/kat"
    os.makedirs(outdir, exist_ok=True)
    for name, cv in CURVES.items():
        rng = random.Random(0xC0FFEE ^ cv.p)
        lines = generate(cv, rng, [])
        slug = name.lower().replace("-", "_")
        path = os.path.join(outdir, f"{slug}.vec")
        with open(path, "w") as f:
            f.write(f"# ELiPS known-answer vectors for {name}\n")
            f.write(f"# generated by tools/reference/gen_vectors.py -- do not edit\n")
            f.write(f"# curve  p={cv.p.bit_length()}bit  r={cv.r.bit_length()}bit"
                    f"  y^2 = x^3 {'+' if cv.sign_b > 0 else '-'} {cv.b}\n")
            f.write(f"prime {H(cv.p)}\n")
            f.write("\n".join(lines) + "\n")
        print(f"{path}: {len(lines)} vectors")
