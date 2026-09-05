"""Is the library's fe_optimal a fixed power of the correct value?

A final exponentiation that computes e^k for a fixed k is still a valid
pairing and is a legitimate optimization -- several published chains do exactly
that. So before calling fe_optimal wrong, find out whether it is simply a
different, consistent normalization.
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, fp12_from_list, fp12_to_list

d = json.load(open(sys.argv[1]))
cv = CURVES[d["curve"]]
p, r, X = cv.p, cv.r, cv.X
L = lambda k: [int(v, 16) % p for v in d[k]]
plain = fp12_from_list(p, L("fe_plain"))
opt   = L("fe_optimal")

cands = {}
for k in list(range(-140, 141)):
    if k: cands[str(k)] = k % r
for i in range(1, 12):
    cands[f"p^{i}"]  = pow(p, i, r)
    cands[f"-p^{i}"] = (-pow(p, i, r)) % r
for nm, v in (("X", X), ("X^2", X*X), ("X-1", X-1), ("(X-1)^2", (X-1)**2),
              ("3(X-1)^2", 3*(X-1)**2), ("X^3", X**3), ("3", 3), ("-3", -3)):
    cands[nm] = v % r
# a few products that show up in published chains
for a in (3, -3, 1, -1):
    for i in range(0, 4):
        cands[f"{a}*p^{i}"] = (a * pow(p, i, r)) % r

hit = None
for nm, k in cands.items():
    if k and fp12_to_list(plain.pow(k)) == opt:
        hit = (nm, k); break

print(f"curve {d['curve']}")
if hit:
    print(f"  fe_optimal == fe_plain ^ ({hit[0]})")
    print(f"  exponent mod r = {hit[1]}")
    print(f"  gcd(k, r) = {__import__('math').gcd(hit[1], r)}  "
          f"(1 means still a non-degenerate pairing)")
else:
    print("  fe_optimal is NOT any tested fixed power of fe_plain")
    print("  tested: +-1..140, +-p^i, X, X^2, (X-1)^2, 3(X-1)^2, small*p^i")
