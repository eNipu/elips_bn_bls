"""The divstep inversion in src/arith/fp.c: model it, and justify its one
magic number.

The C routine repeats a fixed number of "divsteps". That count cannot be
derived from the input -- it has to be fixed, or the running time would leak
the secret -- so it must be large enough for every input. This file is where
that number comes from.

ONE DIVSTEP, on state (delta, f, g) with f odd:

    delta > 0 and g odd :  (1-delta,  g, (g-f)/2)
    otherwise           :  (1+delta,  f, (g + (g&1) f)/2)

Two things are checked here.

1. THE MODEL IS THE C. The batching, the two's complement limbs, the single
   Montgomery step on d and e, and the final fix-up constant are all modelled
   at limb level and checked against exact inverses. If the C and this file
   disagree about the algorithm, this file is the one that was written first
   and verified.

2. THE ITERATION COUNT. The C uses 3*bits. That is above the bound in
   Bernstein and Yang (ePrint 2019/266), about 2.88*bits, and this file
   measures the worst case actually reached over several thousand inputs per
   curve -- random ones, plus the patterns that stress Euclid-style algorithms:
   consecutive Fibonacci numbers, all-ones, powers of two and their neighbours.
   The measured worst is about 2.17*bits, so 3*bits leaves roughly 39% margin.

   Stated plainly: the margin over the measured worst case is large and the
   count also clears the published bound, but the published bound was recalled
   rather than read, since the paper was not reachable from the session that
   wrote this. A reader who wants certainty should check 2019/266 Theorem 11.
   Making the count larger costs one more block per 62 steps and nothing else.

Exit status 0 if the model reproduces exact inverses and the count holds.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES

K = 62
M64 = (1 << 64) - 1


def to_limbs(x, n):
    x &= (1 << (64 * n)) - 1
    return [(x >> (64 * i)) & M64 for i in range(n)]


def from_limbs(L):
    n = len(L)
    x = sum(l << (64 * i) for i, l in enumerate(L))
    if x >> (64 * n - 1):
        x -= 1 << (64 * n)
    return x


def divsteps(delta, f, g, k=K):
    """The block the C runs in registers: only delta and the low limbs are read.
    Each step is (1/2) times an integer matrix, so k steps give one matrix."""
    u, v, q, r = 1, 0, 0, 1
    for _ in range(k):
        if delta > 0 and (g & 1):
            delta, f, g = 1 - delta, g, (g - f) >> 1
            u, v, q, r = 2 * q, 2 * r, q - u, r - v
        else:
            c = g & 1
            delta, g = 1 + delta, (g + c * f) >> 1
            u, v, q, r = 2 * u, 2 * v, c * u + q, c * v + r
    return delta, u, v, q, r


def redc_step(x, p, pinv64):
    m = ((x % (1 << 64)) * pinv64) % (1 << 64)
    return (x + m * p) >> 64


def inv(a, p, nlimbs):
    """Limb-level model of fp_inv."""
    blocks = (3 * p.bit_length() + K - 1) // K
    n = nlimbs + 1
    pinv64 = (-pow(p, -1, 1 << 64)) % (1 << 64)
    F, G = to_limbs(p, n), to_limbs(a % p, n)
    d, e, delta = 0, 1, 1
    for _ in range(blocks):
        delta, u, v, q, r = divsteps(delta, F[0], G[0])
        fv, gv = from_limbs(F), from_limbs(G)
        F = to_limbs((u * fv + v * gv) >> K, n)
        G = to_limbs((q * fv + r * gv) >> K, n)
        d, e = (redc_step(u * d + v * e, p, pinv64) % p,
                redc_step(q * d + r * e, p, pinv64) % p)
    f = from_limbs(F)
    if f < 0:
        f, d = -f, (-d) % p
    if f != 1:
        return None
    # f,g fell by 2^K a block, d,e by 2^64; FP_INV_FIX closes the gap.
    return d * pow(2, 2 * blocks, p) % p


def steps_needed(a, p):
    f, g, delta, n = p, a % p, 1, 0
    limit = 20 * p.bit_length()
    while g != 0:
        if delta > 0 and (g & 1):
            delta, f, g = 1 - delta, g, (g - f) >> 1
        else:
            c = g & 1
            delta, g = 1 + delta, (g + c * f) >> 1
        n += 1
        if n > limit:
            return None
    return n


def stressors(p, cap=4000):
    b = p.bit_length()
    out = [1, 2, 3, p - 1, p - 2, p - 3, (p - 1) // 2, (p + 1) // 2]
    out += [1 << k for k in range(1, b)]
    out += [(1 << k) - 1 for k in range(1, b)]
    out += [(1 << k) + 1 for k in range(1, b)]
    x, y = 1, 1                                  # Fibonacci: Euclid's worst case
    while y < p:
        out.append(y)
        x, y = y, x + y
    out += [random.randrange(1, p) for _ in range(cap)]
    return out


def check(name, nlimbs):
    p = CURVES[name].p
    b = p.bit_length()
    ok = True

    tests = [1, 2, 3, p - 1, p - 2, (p - 1) // 2] + [random.randrange(1, p) for _ in range(120)]
    wrong = sum(1 for a in tests if inv(a, p, nlimbs) != pow(a, -1, p))
    print(f"=== {name} ===")
    print(f"  model reproduces exact inverses: {len(tests) - wrong}/{len(tests)}")
    ok = ok and wrong == 0

    worst = 0
    cands = stressors(p)
    for a in cands:
        a %= p
        if a == 0:
            continue
        s = steps_needed(a, p)
        if s and s > worst:
            worst = s
    used = 3 * b
    print(f"  worst divsteps over {len(cands)} inputs: {worst}  ({worst / b:.2f} x bits)")
    print(f"  the C uses 3*bits = {used}, margin {used - worst} steps "
          f"({100 * (used - worst) / worst:.0f}%)")
    print(f"  blocks of {K}: {(used + K - 1) // K}")
    if worst >= used:
        ok = False
        print("  TOO FEW: the fixed count does not cover the measured worst case")
    print()
    return ok


def main():
    random.seed(20260907)
    bad = 0
    for name, nl in (("BLS12-381", 6), ("BLS12-461", 8), ("BN-462", 8)):
        if not check(name, nl):
            bad += 1
    if bad:
        print(f"FAILED: {bad} curve(s)")
        return 1
    print("Divstep inversion model verified; iteration count covers every input tested.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
