"""BLS12-381 prototype: prepared Miller replay with an 8M2 line kernel.

Run with Python 3, no third-party packages or C build required:
    python3 tools/reference/normalized_miller_ref.py
    python3 tools/reference/normalized_miller_ref.py --samples 2000

API (P and Q are affine elips_ref.EFp / EFp2 points):
    table = prepare_g2(Q)             # once for a fixed Q
    raw = miller_product_scaled8([(P, table)])
    e = final_exponentiate(raw)        # exact e, not raw Miller output
    e3 = final_exponentiate(raw, cubed=True)  # ELiPS fast convention

Reuses elips_ref.py's textbook Fp2/Fp6/Fp12 and affine point arithmetic.
The new kernel, preparation, and shared Miller loop are implemented here.
Preparation deliberately uses slow affine arithmetic for readability; its
timing is NOT a model of the C implementation's projective preparation.

For L = a*yP + c*w^3 + b*xP*w^5, prepare beta=c/a and gamma=b/a.
Online, B=beta/yP and C=gamma*xP/yP give Lhat=1+B*w^3+C*w^5.
The kernel computes 2*f*Lhat in EIGHT Fp2 multiplications, without halving.
All discarded factors lie in Fp2*, which final exponentiation annihilates.
The raw Miller result changes; only the reduced pairing is preserved.

Research only: variable-time Python, deterministic test randomness, no
production security claim, and no claim of publication-level novelty.
"""
import argparse
from dataclasses import dataclass
from pathlib import Path
import random
import sys
from unittest.mock import patch

from elips_ref import CURVES, EFp, EFp2, Fp2, Fp6, Fp12, fp12_from_list


CURVE = CURVES["BLS12-381"]
P_MOD, ORDER = CURVE.p, CURVE.r
FINAL_EXP = (P_MOD**12 - 1) // ORDER
TOP = max(CURVE.loop_digits)
DIGITS = tuple(CURVE.loop_digits.get(i, 0) for i in range(TOP - 1, -1, -1))
LINE_COUNT = sum(1 + (d != 0) for d in DIGITS)
DEFAULT_VECTORS = Path(__file__).resolve().parents[2] / "test/kat/pairing_bls12_381.vec"


@dataclass(frozen=True)
class RawLine:
    a: Fp2
    b: Fp2
    c: Fp2


@dataclass(frozen=True)
class NormalizedLine:
    beta: Fp2                     # c/a, scales 1/yP
    gamma: Fp2                    # b/a, scales xP/yP


def _twice_mul_linear(t, B, C, plus, minus):
    """2*t*(B+C*v) mod v^3-xi, using four Fp2 multiplications.

    Evaluate the degree-three product at 0, +1, -1 and infinity.
    Interpolate TWICE the coefficients, so no division by 2 is needed.
    plus=B+C and minus=B-C are shared between the two calls per line.
    """
    m0 = t.c0 * B
    m3 = t.c2 * C
    s = t.c0 + t.c2
    mp = (s + t.c1) * plus
    mm = (s - t.c1) * minus
    return Fp6(
        (m0 + m3.mul_xi()).mul_int(2),
        mp - mm - m3.mul_int(2),
        mp + mm - m0.mul_int(2),
    )


def mul_line_scaled8(f, B, C):
    """Return 2*f*(1+B*w^3+C*w^5), NOT the exact unscaled product.

    All operands must belong to the same field tower as elips_ref.py.
    For f=A+D*w, w^2=v:
      2*f*Lhat = (2*A + v^2*H(D)) + (2*D + v*H(A))*w
      H(t) = 2*t*(B+C*v).
    Multiplication by xi/v is just additions and coordinate permutations.
    """
    plus, minus = B + C, B - C
    hA = _twice_mul_linear(f.d0, B, C, plus, minus)
    hD = _twice_mul_linear(f.d1, B, C, plus, minus)
    return Fp12(f.d0 + f.d0 + hD.mul_v().mul_v(),
                f.d1 + f.d1 + hA.mul_v())


def batch_invert(values):
    """Invert nonzero Fp2 values with 1I2 + 3*(n-1)M2; empty is allowed.

    Reject zero explicitly: the reference field's inv(0) returns zero,
    which would silently invalidate line normalization.
    """
    if not values:
        return []
    p = values[0].p
    if any(v.p != p or v.is_zero() for v in values):
        raise ValueError("batch inversion needs nonzero values in one field")
    prefix = [values[0]]
    for v in values[1:]:
        prefix.append(prefix[-1] * v)
    inv = prefix[-1].inv()
    out = [None] * len(values)
    for i in range(len(values) - 1, 0, -1):
        out[i] = inv * prefix[i - 1]
        inv = inv * values[i]
    out[0] = inv
    return out


def _validate_point(point, cls):
    if not isinstance(point, cls) or point.cv is not CURVE:
        raise ValueError("expected a BLS12-381 point in the correct group")
    if point.inf or not point.is_on_curve() or not point.mul(ORDER).inf:
        raise ValueError("point must be nonidentity and in the order-r subgroup")


def _raw_lines(Q):
    """Generate Q-only line coefficients using affine textbook arithmetic.

    At a tangent: (a,b,c) = (xi*2y, -3x^2, 3x^3-2y^2).
    At a chord:   (a,b,c) = (xi*H, -R, R*x-y*H).
    These clear denominators without the C code's projective scale factors.
    Q must already be validated.
    """
    T = Q if CURVE.loop_digits[TOP] > 0 else -Q
    lines = []
    for d in DIGITS:
        xx = T.x.sqr()
        lines.append(RawLine(T.y.mul_int(2).mul_xi(),
                             -xx.mul_int(3),
                             (xx * T.x).mul_int(3) - T.y.sqr().mul_int(2)))
        T = T.dbl()
        if d:
            S = Q if d > 0 else -Q
            H, R = S.x - T.x, S.y - T.y
            lines.append(RawLine(H.mul_xi(), -R, R * T.x - T.y * H))
            T = T + S
    return tuple(lines)


def _normalize_lines(lines):
    inverses = batch_invert([line.a for line in lines])
    return tuple(NormalizedLine(line.c * ai, line.b * ai)
                 for line, ai in zip(lines, inverses))


def prepare_g2(Q):
    """Validate a fixed Q and return its reusable normalized line table."""
    _validate_point(Q, EFp2)
    return _normalize_lines(_raw_lines(Q))


def miller_product_scaled8(terms):
    """Shared Miller loop for iterable (P, prepared_Q) pairs.

    Tables must come from prepare_g2; they are trusted local objects, not an
    encoding that can safely be accepted from an attacker. P is validated
    here. The empty product is one. No eight-term chunk limit is imposed.
    """
    prepared = []
    for P, table in terms:
        _validate_point(P, EFp)
        if len(table) != LINE_COUNT:
            raise ValueError("wrong prepared-table length")
        iy = pow(P.y, -1, P_MOD)
        prepared.append((table, iy, P.x * iy % P_MOD))
    f, k = Fp12.one(P_MOD), 0
    for d in DIGITS:
        f = f.sqr()                     # shared by every term
        for table, iy, xy in prepared:
            for index in range(k, k + 1 + (d != 0)):
                line = table[index]
                f = mul_line_scaled8(f, line.beta.mul_int(iy),
                                     line.gamma.mul_int(xy))
        k += 1 + (d != 0)
    return f


def final_exponentiate(f, cubed=False):
    """Definition-level final exponent, not an optimized addition chain.

    Default: exact pairing e. cubed=True: e^3, matching ELiPS's fast
    BLS12 convention. The extra nonzero Fp2 factors vanish in both cases.
    """
    if f.p != P_MOD or (f.d0.is_zero() and f.d1.is_zero()):
        raise ValueError("final exponentiation needs nonzero BLS12-381 Fp12")
    return f.pow(FINAL_EXP * (3 if cubed else 1))


def _dense_line(a, b, c):
    z = Fp2.zero(a.p)
    return Fp12(Fp6(a, z, z), Fp6(z, b, c))


def _dense_miller(P, lines):
    """Control: unnormalized lines and ordinary dense Fp12 multiplication."""
    f, k = Fp12.one(P_MOD), 0
    for d in DIGITS:
        f = f.sqr()
        for _ in range(1 + (d != 0)):
            line = lines[k]
            f = f * _dense_line(line.a.mul_int(P.y), line.c,
                                line.b.mul_int(P.x))
            k += 1
    return f


def _check(condition, message):
    # Unlike assert, these self-tests remain effective under python -O.
    if not condition:
        raise AssertionError(message)


def _expect_value_error(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("invalid input was accepted")


def _load_vectors(path):
    records, prime = [], None
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "prime":
            prime = int(fields[1], 16)
        elif fields[0] == "pair":
            if len(fields) != 20 or fields[7] != "=":
                raise ValueError("malformed pairing vector")
            vals = [int(s, 16) for s in fields[1:7] + fields[8:]]
            if any(v < 0 or v >= P_MOD for v in vals):
                raise ValueError("noncanonical vector coefficient")
            px, py, qx0, qx1, qy0, qy1 = vals[:6]
            records.append((EFp(CURVE, px, py),
                            EFp2(CURVE, Fp2(P_MOD, qx0, qx1),
                                 Fp2(P_MOD, qy0, qy1)),
                            fp12_from_list(P_MOD, vals[6:])))
        else:
            raise ValueError("unexpected vector record")
    if prime != P_MOD or len(records) != 4:
        raise ValueError("expected the four BLS12-381 pairing vectors")
    return records


def selftest(samples, vector_path):
    rng = random.Random(20260909)
    one, zero = Fp2.one(P_MOD), Fp2.zero(P_MOD)

    def random2():
        return Fp2(P_MOD, rng.randrange(P_MOD), rng.randrange(P_MOD))

    def random12():
        return Fp12(Fp6(random2(), random2(), random2()),
                    Fp6(random2(), random2(), random2()))

    # Exact integer identity, independent of all finite-field code.
    for i in range(3):
        for j in range(2):
            t = [int(k == i) for k in range(3)]
            B, C = int(j == 0), int(j == 1)
            m0, m3 = t[0] * B, t[2] * C
            mp = sum(t) * (B + C)
            mm = (t[0] - t[1] + t[2]) * (B - C)
            got = [2*m0, mp-mm-2*m3, mp+mm-2*m0, 2*m3]
            expected = [0] * 4
            expected[i+j] = 2
            _check(got == expected, "integer interpolation identity")
    _check(FINAL_EXP % (P_MOD**2 - 1) == 0, "Fp2 factor is not annihilated")
    _check(LINE_COUNT == 69 and TOP == 64, "unexpected loop parameters")
    T = CURVE.loop_digits[TOP]
    for d in DIGITS:
        _check(T % ORDER != 0 and (2*T) % ORDER != 0, "degenerate tangent")
        T *= 2
        if d:
            _check((T-1) % ORDER != 0 and (T+1) % ORDER != 0, "degenerate chord")
            T += d
    _check(T == CURVE.X, "signed loop does not compute x")
    print("PASS integer identities and all normalization preconditions")

    for i in range(samples):
        f, B, C = random12(), random2(), random2()
        if i % 5 == 0:
            B = zero
        if i % 7 == 0:
            C = zero
        if i % 11 == 0:
            C = B
        if i % 13 == 0:
            C = -B
        if i % 17 == 0:
            f = Fp12.zero(P_MOD)
        if i % 19 == 0:
            f = Fp12.one(P_MOD)
        expected = f * _dense_line(one, B, C)
        _check(mul_line_scaled8(f, B, C) == expected + expected,
               f"kernel sample {i}")

    # Count actual operator calls rather than merely printing a claimed cost.
    count = 0
    original_mul = Fp2.__mul__

    def counted_mul(a, b):
        nonlocal count
        count += 1
        return original_mul(a, b)

    f, B, C = random12(), random2(), random2()
    with patch.object(Fp2, "__mul__", counted_mul):
        mul_line_scaled8(f, B, C)
    _check(count == 8, f"kernel used {count}M2 instead of 8M2")
    _check(batch_invert([]) == [], "empty batch inversion")
    for n in (1, 2, 69):
        values = [Fp2(P_MOD, rng.randrange(1, P_MOD), rng.randrange(P_MOD))
                  for _ in range(n)]
        count = 0
        with patch.object(Fp2, "__mul__", counted_mul):
            inverses = batch_invert(values)
        _check(count == 3*(n-1), "batch-inversion multiplication count")
        _check(all(a*b == one for a, b in zip(values, inverses)),
               "batch inversion")
    inv_count = 0
    original_inv = Fp2.inv

    def counted_inv(a):
        nonlocal inv_count
        inv_count += 1
        return original_inv(a)

    lines = tuple(RawLine(a, random2(), random2()) for a in values)
    count = 0
    with patch.object(Fp2, "__mul__", counted_mul), patch.object(Fp2, "inv", counted_inv):
        normalized = _normalize_lines(lines)
    _check(count == 342 and inv_count == 1, "table normalization cost")
    _check(all(n.beta * line.a == line.c and n.gamma * line.a == line.b
               for line, n in zip(lines, normalized)), "normalized coefficients")
    _expect_value_error(lambda: batch_invert([one, zero]))
    _expect_value_error(lambda: final_exponentiate(Fp12.zero(P_MOD)))
    print(f"PASS {samples} kernel checks; counted 8M2 and 1I2+342M2 normalization")

    records = _load_vectors(vector_path)
    terms, expected_values = [], []
    for i, (P, Q, expected) in enumerate(records, 1):
        table = prepare_g2(Q)
        raw = miller_product_scaled8([(P, table)])
        dense = _dense_miller(P, _raw_lines(Q))
        ratio = raw * dense.inv()
        _check(not ratio.d0.c0.is_zero() and ratio.d0.c1.is_zero()
               and ratio.d0.c2.is_zero() and ratio.d1.is_zero(),
               f"KAT {i}: raw discrepancy is not in Fp2*")
        _check(final_exponentiate(raw) == expected, f"KAT {i}: exact pairing")
        if i == 1:
            _check(raw != dense, "test did not witness changed raw Miller value")
            _check(not expected.is_one() and expected.pow(ORDER).is_one(),
                   "degenerate or wrong-order pairing")
            _check(final_exponentiate(raw, cubed=True) == expected.pow(3),
                   "cubed pairing convention")
        terms.append((P, table))
        expected_values.append(expected)
        print(f"PASS pairing KAT {i}/4 and nonzero Fp2 discrepancy")

    for n in (0, 2, 4, 9):
        selected = [terms[i % len(terms)] for i in range(n)]
        expected = Fp12.one(P_MOD)
        for i in range(n):
            expected = expected * expected_values[i % len(terms)]
        _check(final_exponentiate(miller_product_scaled8(selected)) == expected,
               f"shared loop with {n} terms")
    P, Q, e = records[0]
    table = terms[0][1]
    _check(final_exponentiate(miller_product_scaled8([(P, table), (-P, table)])).is_one(),
           "cancelling pairing product")
    s, t = rng.randrange(1, ORDER), rng.randrange(1, ORDER)
    raw = miller_product_scaled8([(P.mul(s), prepare_g2(Q.mul(t)))])
    _check(final_exponentiate(raw) == e.pow(s*t % ORDER),
           "full-width random scalar bilinearity")
    _expect_value_error(lambda: prepare_g2(EFp2(CURVE)))
    _expect_value_error(lambda: prepare_g2(EFp2(CURVE, zero, zero)))
    _expect_value_error(lambda: prepare_g2(EFp(CURVE, P.x, P.y)))
    _expect_value_error(lambda: miller_product_scaled8([(EFp(CURVE), table)]))
    _expect_value_error(lambda: miller_product_scaled8([(EFp(CURVE, 0, 0), table)]))
    # (0,2) is on y^2=x^3+4 but has order 3, not order r.
    _expect_value_error(lambda: miller_product_scaled8([(EFp(CURVE, 0, 2), table)]))
    _expect_value_error(lambda: miller_product_scaled8([(P, table[:-1])]))
    print("PASS shared products (0,2,4,9 terms), cancellation, bilinearity, invalid inputs")
    print("All checks passed. Raw Miller values differ; reduced pairings agree.")
    print("C cost model (1M2=3Mp): 15M2 -> 8M2 per line; 1449Mp saved over 69 lines.")
    print("Charge an extra 1Ip + 1Mp per online P, and 1I2 + 342M2 per table normalization.")
    print("Affine preparation and Python timings are not C performance estimates.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", type=int, default=1000,
                        help="full-width randomized kernel tests (default: 1000)")
    parser.add_argument("--vectors", type=Path, default=DEFAULT_VECTORS,
                        help="four-record BLS12-381 pairing KAT file")
    args = parser.parse_args()
    if args.samples < 32:
        parser.error("--samples must be at least 32 to exercise edge cases")
    try:
        selftest(args.samples, args.vectors)
    except (AssertionError, ValueError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
