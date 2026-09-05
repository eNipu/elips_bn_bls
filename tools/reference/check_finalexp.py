"""Which of the library's two final exponentiations is correct? (issues #13, #3)

Takes the library's own raw Miller output and raises it to (p^12-1)/r directly.
That is the definition, so whichever library routine matches it is the correct
one. Also recomputes the whole pairing from the dumped points, to say whether
the Miller loop itself agrees.
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES, Fp2, fp12_from_list, fp12_to_list
from pairing_ref import optimal_ate

d = json.load(open(sys.argv[1]))
cv = CURVES[d["curve"]]
p, r = cv.p, cv.r
assert int(d["p"], 16) == p, "prime mismatch between dump and reference"
assert int(d["r"], 16) == r, "order mismatch between dump and reference"

L = lambda k: [int(v, 16) % p for v in d[k]]
raw, plain, opt = L("miller_raw"), L("fe_plain"), L("fe_optimal")

e = (p**12 - 1) // r
direct = fp12_to_list(fp12_from_list(p, raw).pow(e))

print(f"curve {d['curve']}   exponent (p^12-1)/r is {e.bit_length()} bits")
print(f"  library fe_plain    == raw^((p^12-1)/r) : {plain == direct}")
print(f"  library fe_optimal  == raw^((p^12-1)/r) : {opt   == direct}")
print(f"  fe_plain == fe_optimal                  : {plain == opt}")

D = fp12_from_list(p, direct)
print(f"  raw^((p^12-1)/r) in mu_r                : {D.pow(r).is_one()}")
print(f"  raw^((p^12-1)/r) non-degenerate         : {not D.is_one()}")

for nm, val in (("fe_plain", plain), ("fe_optimal", opt)):
    V = fp12_from_list(p, val)
    print(f"  {nm:10s} in mu_r: {V.pow(r).is_one()}   non-degenerate: {not V.is_one()}")

# Does the reference Miller loop agree, once both use the same final exponent?
ref = fp12_to_list(optimal_ate(cv, int(d["px"],16), int(d["py"],16),
                               Fp2(p, int(d["qx0"],16), int(d["qx1"],16)),
                               Fp2(p, int(d["qy0"],16), int(d["qy1"],16))))
print(f"  reference pairing == raw^((p^12-1)/r)   : {ref == direct}")
print(f"  reference pairing == fe_plain           : {ref == plain}")
print(f"  reference pairing == fe_optimal         : {ref == opt}")
