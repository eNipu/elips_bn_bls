# Progress log

Running record of what is done, what it proved, and where to pick up.
Newest phase last. Plan lives in `MODERNIZATION_PLAN.md`; tasks are GitHub issues.

---

## Phase 0 — Build the test oracle and baseline (issue #1)

**Status: substantially complete. Exit gate met, with one finding that changes Phase 2.**

### What exists now

| Artifact | Purpose |
|---|---|
| `tools/reference/elips_ref.py` | Independent Fp2/Fp6/Fp12 tower and curve arithmetic, written from the defining equations |
| `tools/reference/selftest.py` | Field axioms, tower relations, `[r]P = O`; gates vector generation |
| `tools/reference/gen_vectors.py` | Emits 857 known-answer vectors per curve |
| `tools/reference/pairing_ref.py` | Textbook optimal ate, computed entirely in Fp12 on untwisted points |
| `tools/reference/check_finalexp.py` | Decides which library final exponentiation is correct |
| `tools/reference/trace_finalexp.py` | Traces each chain over the exponent ring to get its exact exponent |
| `test/kat_runner.c` | Drives the library from a vector file; real exit code |
| `test/dump_pairing.c`, `test/diag_pairing.c` | Export the library's own inputs and outputs for comparison |
| `test/kat/*.vec` | Committed vectors for BN-462, BLS12-461, BLS12-381 |
| `bench/baseline.json` | Baseline timings |

### Results

**Field and curve layers are correct.** 857 vectors x 3 curves, all passing:

```
test/kat/bn_462.vec:     857 passed, 0 failed, 1 non-canonical
test/kat/bls12_461.vec:  857 passed, 0 failed, 1 non-canonical
test/kat/bls12_381.vec:  857 passed, 0 failed, 1 non-canonical
```

The runner was verified to actually fail — corrupt one vector and it exits 1 naming the record; missing file or unknown op exits 2. An oracle that cannot fail is not an oracle.

The single non-canonical value per curve is `Fp_neg(0)`, which returns `p` rather than `0`. Congruent but unreduced, and `Fp_cmp` compares integers rather than residues. Logged for Phase 2.

**Miller loops are correct.** The independent reference agrees with the library's raw Miller output on both BN-462 and BLS12-461, once both use the same final exponent.

**Final exponentiation is not.** See issue #16. `finalexp_plain` is correct; `finalexp_optimal` — which every public pairing entry point actually calls — computes `e^3` on BLS12 and `e^(~12X^3)` on BN. Both remain bilinear and non-degenerate, `gcd(q,r)=1`, which is exactly why the existing suite never caught it. Confirmed two independent ways: empirically against the reference, and symbolically by tracing the chain over the exponent ring, with the trace then re-validated against measured values.

### Enabling result worth remembering

The library's tower collapses to a single polynomial. From `v = w^2` and `1+u = v^3 = w^6` we get `u = w^6 - 1`, hence

```
Fp12 = Fp[w] / (w^12 - 2*w^6 + 2)
```

verified irreducible over all three primes. Sage can therefore build **the same field**, not an abstract `GF(p^12)`, making its pairing directly comparable coefficient by coefficient. Sage computes a BLS12-381 Tate pairing in ~0.1 s, so it is cheap enough to use routinely.

### Hypotheses formed and disproved

Four suspected bugs did not survive checking. Each would have caused a regression if "fixed":

1. Miller loop bounds on both curves — the leading signed digit is absorbed by the accumulator initialization.
2. `Fp2_set_ui` broadcasting to both coefficients — relied on by `set_basis()` to build `1+u`.
3. `char str[5]` in the split-scalar routines — fits its four characters and terminator exactly.
4. Global mutation in `bls12_finalexp_optimal` (defect A4) — measured to be idempotent and to restore `bls12_X_length` and `bls12_X_binary[]` correctly. Still non-reentrant and thread-hostile, but it does not corrupt results on the normal path.

### Remaining in Phase 0

- [x] #8 — superseded by the CMake build in Phase 1; tag `phase0-reference` created at `987ed77`
- [x] #14 — property tests are now a committed CTest target (`test/property_tests.c`)
- [ ] #15 — replace the single-sample baseline with repetitions and spread (deferred to Phase 4, where the speedup gates need it)

### Reproduce

```bash
python3 tools/reference/selftest.py
python3 tools/reference/gen_vectors.py test/kat
clang -O2 -Iinclude -I/opt/homebrew/include src/*.c test/kat_runner.c \
      -L/opt/homebrew/lib -lgmp -o kat_runner
./kat_runner test/kat/bn_462.vec
python3 tools/reference/trace_finalexp.py
```

### Hand-off note

Phase 2 gains a defect beyond the 20 in the plan: issue #16. It needs a maintainer
decision before it can be fixed, since the BLS12 factor of 3 may be a deliberate
Hayashida–Hayasaka–Teruya style optimization. The BN factor of roughly `12·X^3` is
harder to read as intentional.

Phase 1 can start regardless — it touches the build system, not the mathematics.


---

## Phase 1 — CMake build system and repository hygiene (issue #2)

**Status: complete. Exit gate met.**

### What exists now

| Artifact | Purpose |
|---|---|
| `CMakeLists.txt` | Library, tests, install, package config |
| `cmake/FindGMP.cmake` | Locates GMP, exports `GMP::GMP` |
| `cmake/CorruptCheck.cmake` | Negative control driving the vector runner |
| `test/property_tests.c` | Subgroup, non-degeneracy, bilinearity, additivity |
| `test/finalexp_agreement.c` | Pins issue #16 as an expected failure |
| `.github/workflows/ci.yml` | Matrix build and test |
| `.github/workflows/docs.yml` | Doxygen to gh-pages |

### Results

Eight CTest targets, all green on a clean configure:

```
kat.bn_462            kat.bls12_461         kat.bls12_381
kat.detects_corruption
property.bn           property.bls12
finalexp.agreement.bn finalexp.agreement.bls12   (WILL_FAIL, issue #16)
```

Property tests pass on both curves, including `[r]P = O` and `[r]Q = O`, which
nothing in the repository checked before.

CI covers `ubuntu-latest` (x86-64) and `macos-latest` (AArch64) across Release,
Asan and Ubsan, so both Phase 6 assembly targets are exercised from the start.
A separate job re-runs the reference self-test and regenerates the vectors, so
the committed ones cannot drift unnoticed.

### Decisions taken

**MemorySanitizer omitted.** It needs every dependency instrumented and GMP is
not, so it would produce false positives rather than findings. Revisit in
Phase 3 when the arithmetic moves off `mpz_t`.

**`-DELIPS_CURVE` deferred to Phase 3.** Compile-time curve selection only means
something once `FP_LIMBS` is a constant. Today the curve is chosen at runtime by
`init_bn()` / `bls12_inits()`, so the switch would be scaffolding with nothing
behind it.

**Autotools inputs deleted, not kept alongside.** `configure.ac`, both
`Makefile.am` and `INSTALL` are gone. Two build systems guarantee one rots.

### Corrections to earlier assumptions

**ASan does not catch defect M7.** The Phase 1 gate assumed sanitizers would
surface the memory defects. They do not surface this one, because ASan does not
instrument variable-length arrays. M7 was instead confirmed directly:
`mpz_get_str` writes its NUL terminator at index `length`, exactly one past a
`char binary[length]` VLA, for every value tested. The defect is real; the
detection method in the plan was wrong. Phase 2 must fix it by review, not by
waiting for a sanitizer report.

### Open question for the maintainer

`docs/` is still tracked: 1280 of the 1406 remaining files. GitHub Pages serves
from `master:/docs`, so untracking it takes
https://enipu.github.io/elips_bn_bls/ offline. `.github/workflows/docs.yml`
regenerates and publishes to `gh-pages`; once Pages is repointed there,

```bash
git rm -r --cached docs && git commit -m "Untrack the generated Doxygen site"
```

drops the repository to about 126 tracked files. Not done unilaterally because
it is outward-facing.

### Next

Phase 2 (issue #3) — fix the confirmed defects on the existing architecture.
Now unblocked, with an oracle and CI in place to catch regressions. Note that
Phase 2 inherits issue #16, which needs a maintainer decision first.

---

## Phase 2 — Fix the confirmed defects (issue #3)

**Status: complete. Exit gate met.**

### Fixed

All 20 defects from the audit, minus A3 and A4 which are structural and belong
to Phase 3, plus three found during the work.

| ID | Defect |
|---|---|
| M7 | `mpz_get_str` writes its NUL at index `length`; the `char binary[length]` VLA overflowed by one byte at 14 sites |
| **M8** | **new** — split-scalar `memmove` copied `sizeof(buffer)` instead of the digit count, writing past the row for any short sub-scalar. Six sites |
| M1 | `bls12_generate_prime` called `mpz_init` where it meant `mpz_clear`, on all three exit paths |
| M2 | `EFpd_total` written by `weil()` and read by both G2 generators, never initialised |
| M3 | `curve_a` read and released without ever being initialised |
| M4 | `f_ltq` leaked an `Fp2` per call, inside the Miller loop |
| M5 | `bls12_finalexp_optimal` initialised two unused `Fp12` and leaked them |
| M6 | `bls12_4split_G2_scm` called `EFp2_init` in its cleanup block |
| **M9** | **new** — `clear_parameters` missed three fields; `init_precoms` had no teardown at all |
| **M10** | **new** — 18 shadowing `gmp_randstate_t`, none released; each call leaked a Mersenne Twister state |
| A1 | `Fp2_mul_Fp` read its output instead of its input |
| A2 | Three `Fp12` routines did the same |
| A5 | `Fp_inv` and `Fp_div` discarded `mpz_invert`'s status; `Fp_div` had no declaration |
| A6 | `bls12_generate_order` inverted the success convention; callers ignored it |
| A7 | `gmp_printf` argument with no conversion specifier |
| A8 | Duplicate declaration across two headers |
| X1–X5 | Retired wholesale by deleting the KSS16 tree |

### Measured

| Check | Before | After |
|---|---|---|
| Build warnings | 8 | **0** |
| Leaks, 10 / 40 / 160 pairings | 48 / 168 / 648 | **0 / 0 / 0** |
| Leaks, 2 / 8 / 32 point generations | 5 / 17 / 65 | **0 / 0 / 0** |
| Pairing time | 8.614 ms | 8.678 ms (within noise, as required) |
| Vectors | 857 x 3 pass | 857 x 3 pass |
| CTest | 8/8 | 8/8 |
| ASan + UBSan | clean | clean |

Leaks previously grew linearly with work. They are now zero at every scale.

The KSS16 deletion removed 1576 lines, none of it reachable from any pairing.

### New tests

`property_tests.c` gained a split-scalar consistency check: 2-split and 4-split
G1/G2 scalar multiplication and G3 exponentiation must agree with the plain
versions over random scalars. This guards the M8 fix, which changed the
semantics of the memmove that builds the window indices.

### Corrections to earlier assumptions

**ASan does not catch M7 or M8.** Both are VLA overflows and ASan does not
instrument VLAs. M7 was confirmed by showing `mpz_get_str` writes at index
`length` for every value tested; M8 by computing the write range against the row
size and measuring how often real scalars trigger it (about two thirds). The
plan's assumption that sanitizers would surface the memory defects was wrong for
exactly the two that matter most.

**M10 was invisible to the pairing-scaling leak test.** Leaks looked constant at
43 because point generation happens once per run. Scaling *point generation*
rather than pairings exposed it immediately. Worth remembering: scale the thing
that allocates, not the thing that looks expensive.

### Deliberately not done

**Issue #16** (final exponentiation raises to the wrong exponent) is untouched.
It needs a maintainer decision first, because the BLS12 factor of 3 may be a
deliberate Hayashida–Hayasaka–Teruya optimization. The expected-failure test in
CTest keeps it visible.

**A3 and A4** stay for Phase 3, which deletes the globals that cause them.

### Next

Phase 3 (issue #4) — fixed-width `mpn` representation and Montgomery arithmetic.
This is the phase that has to deliver a measured 3x or the plan's model of where
time goes is wrong.

---

## Phase 3 — Fixed-width mpn representation and Montgomery arithmetic (issue #4)

**Status: foundation complete and measured. The rest is blocked on a plan
correction, described below. Exit gate NOT yet evaluable.**

### What exists now

| Artifact | Purpose |
|---|---|
| `include/elips/fp_params.h` | Generated Montgomery parameters for all three curves |
| `include/elips/fp.h` | The new field API |
| `src/arith/fp.c` | CIOS Montgomery multiply, constant-time add/sub/select, inversion |
| `tools/reference/gen_params.py` | Emits the parameters; hand-typing an `R^2` limb is not survivable |
| `test/fp_difftest.c` | Differential test against GMP's mpz |

Built once per curve in CMake (`elips_arith_BLS12_381` and friends) with a
`fp.diff.<curve>` test each, so all three are exercised on every CI run.

### Correctness

1,204,813 checks per curve, zero failures, on BLS12-381, BLS12-461 and BN-462.
Includes the boundary values Montgomery code actually gets wrong: 0, 1, `p-1`,
`p-2`, `p/2`, every pairing of those, and the Montgomery round-trip identity.
ASan and UBSan clean.

### Measured, BLS12-461, Apple Silicon, clang -O2

| Operation | old mpz | new Montgomery | change |
|---|---|---|---|
| multiply | 0.1527 us | 0.0716 us | **2.13x faster** |
| add | 0.0425 us | 0.0062 us | **6.90x faster** |
| inverse | 1.4462 us | 26.7978 us | **19x slower** |
| `Fp_init`+`Fp_clear` | 0.0028 us | none | eliminated |
| `Fp12_init`+`Fp12_clear` | 0.0361 us | none | eliminated |

### Two corrections to the plan, both from measurement

**1. Constant-time inversion is not affordable in an affine Miller loop.**

Three inversion strategies, measured rather than assumed:

| Method | Time | Constant time |
|---|---|---|
| `mpz_invert` | 1.4 us | no |
| `mpn_sec_invert` | 26.4 us | yes |
| `a^(p-2)` Fermat chain | 76.3 us | yes |

The plan (Appendix B) said to start with the Fermat chain and only reach for
something better if a measurement demanded it. The measurement demands it: the
chain is the worst of the three, and GMP's `mpn_sec_invert` is both faster and a
tenth of the code, so `fp_inv` uses that.

The ordering consequence is the real point. The BLS12 Miller loop performs about
79 inversions, one per iteration plus the correction steps. At 26.8 us each that
is **2.1 ms of inversion alone**, against a current whole-Miller-loop cost of
2.49 ms. Turning on constant-time inversion before inversions leave the loop
would roughly double the pairing time, not improve it.

So **Phase 4's projective coordinates are a prerequisite for Phase 3's gate, not
a follow-on.** The plan has them in the wrong order. Recommended fix: pull
projective coordinates into Phase 3, or explicitly allow `fp_inv_vartime` in the
Miller loop until Phase 4 lands. `fp_inv_vartime` exists and is documented as
never-for-secrets.

**2. The Phase 0 model over-weighted inversion.**

Section 3 of the plan estimated modular inversion at 20-25% of pairing runtime.
Measured: 79 inversions x 1.45 us is about 115 us against a 2490 us Miller loop,
so roughly **5%**. The estimate was out by a factor of four or five.

Allocation and raw multiply count dominate instead, which is consistent with the
2.13x and 6.90x above and with allocation disappearing entirely. The direction
of the plan is unaffected; the attribution was wrong.

### Projection for the gate

An Fp12 multiplication is roughly 54 Fp multiplies and 100 Fp additions.

| | old | new |
|---|---|---|
| 54 multiplies | 8.25 us | 3.87 us |
| 100 additions | 4.25 us | 0.62 us |
| total | 12.50 us | 4.49 us |

That is **2.8x from arithmetic alone**, before counting eliminated allocation or
lazy reduction, neither of which is implemented yet. The >= 3x gate looks
reachable but is not comfortable, and it depends on lazy reduction actually
delivering. This is exactly the situation the gate was written for: do not
assume, measure, and stop if it comes in under.

### Ordering decision

Resolved: projective coordinates were pulled into Phase 3. See below.


---

## Phase 3, part 2 — tower and Jacobian curve arithmetic

**Status: field tower and curve layer complete and verified. Miller loop and
final exponentiation remain, so the >= 3x pairing gate is still not evaluated.**

### What exists now

| Artifact | Purpose |
|---|---|
| `src/arith/fpx.c` | fp2/fp6/fp12 on fixed-width Montgomery `fp_t` |
| `src/arith/ec.c`, `src/arith/ec_tmpl.h` | Jacobian group law, one template, two instantiations |
| `test/kat_runner_new.c` | Drives the NEW layer through the Phase 0 vectors |
| `test/ec_test.c` | Group-law properties and the subgroup check |

The curve template is deliberate. The original library kept two hand-copied
versions of every curve routine and they had already drifted apart; generating
`ep` and `ep2` from one text removes that failure mode.

### Verification

- 857 vectors pass on all three curves **against the new layer**, using the same
  files that validated the legacy one. Both implementations are therefore
  checked against the independent oracle, not against each other.
- 12 group-law properties pass, including `[r]P == O` for a real G2 generator.
- ASan and UBSan clean. 15 CTest targets green.

### Measured, BLS12-461

| Operation | legacy | new | change |
|---|---|---|---|
| Fp12 multiply | 20.03 us | 5.62 us | **3.56x** |
| Fp12 square | 14.24 us | 5.63 us | 2.55x |
| Fp12 inverse | 52.81 us | 37.61 us | 1.41x |
| EC double | 6.14 us | 1.55 us | **4.00x** |
| EC add, complete | 5.39 us | 5.53 us | 0.98x |
| EC add, generic | 5.39 us | 3.91 us | 1.38x |

Fp12 multiply beat the 2.8x projection because eliminating allocation helps more
than the arithmetic alone predicted.

### Two things measurement changed

**Complete addition is not worth its cost in the Miller loop.** To resolve the
coincident-point case without branching it computes a full doubling on every
call and selects, which gives back everything Jacobian coordinates gained: 0.98x
against the affine code it replaces. `ep_add_generic` drops the fixup. It is
still branch-free and therefore still constant time, but it is wrong for
coincident or opposite points, so it is documented for the Miller loop only,
where the operands are distinct by construction. The EC test initially failed
against it, correctly, because the test itself violated that precondition.

**Squaring needs its own formula.** `fp12_sqr` currently calls `fp12_mul`, so it
measures 2.55x against a legacy layer that has a dedicated squaring. Worth
writing before the final measurement, since the Miller loop squares once per
iteration.

### Another latent defect found in the legacy layer

`EFp2_rational_point` generates points on `y^2 = x^3 - b` over Fp2, not on the
sextic twist `y^2 = x^3 + b*xi` that actually carries G2. Measured residual is
`p-4`, where a real G2 point gives `4+4u`. Nothing broke, because the group law
never references `b` and the function is only used by tests. Same class of
latent error as the rest of the audit; recorded rather than fixed, since Phase 3
replaces the routine.

### Remaining in Phase 3

- [ ] Dedicated `fp12_sqr`, and cyclotomic squaring for the final exponentiation
- [ ] Line functions in Jacobian coordinates, and sparse fp12 multiplication
- [ ] Frobenius constants for the new tower
- [ ] One parameterised Miller loop for both families
- [ ] Final exponentiation (and the issue #16 decision feeds in here)
- [ ] Curve context struct, deleting the globals and with them A3 and A4
- [ ] Re-measure the full pairing against the >= 3x gate

All the primitives the Miller loop needs now exist and are measured, so the
remaining work is assembly of known parts rather than exploration.

---

## Phase 3, part 3 — Miller loop and the gate

**Status: BLS12 pairing complete on the new stack. The >= 3x gate is met on a
like-for-like comparison. Absolute wall-clock is not yet better; see below.**

### What exists now

| Artifact | Purpose |
|---|---|
| Frobenius constants in `fp_params.h` | Generated in Montgomery form, so no init step and no globals |
| `fp12_frobenius`, `fp12_exp` | The p, p^2, p^3 maps and public-exponent exponentiation |
| dedicated `fp12_sqr` | Two fp6 multiplies instead of three |
| `src/pairing/miller.c` | Line functions, sparse multiply, Miller loop, final exponentiation |
| `test/pairing_test.c` | The new pairing against the Phase 0 reference value |

### Correctness

The new pairing equals the value Phase 0 established as correct, coefficient by
coefficient, on the first run. Also non-degenerate and in mu_r. 16 CTest targets
green, ASan and UBSan clean, zero warnings on a clean build.

Frobenius was checked against explicit exponentiation by `p^k` on all three
curves, plus the homomorphism property and twelve-applications-is-the-identity.

The line-function derivation independently reproduced the `(0,3,5)` sparsity the
legacy code uses. Arriving at the same pattern from the defining equations is
good evidence the twist conventions were read correctly.

### Measured, BLS12-461, like for like

| | legacy | new | speedup |
|---|---|---|---|
| Miller loop | 2610.8 us | 856.0 us | **3.05x** |
| final exponentiation | 37977.1 us | 10095.9 us | **3.76x** |
| whole pairing | 40587.9 us | 10951.9 us | **3.71x** |

Both sides of the final exponentiation row use the same definitionally-correct
exponent, so this compares implementations rather than algorithms.

### What this does not say

**The new stack is not yet faster in absolute terms.** The legacy fast path runs
a pairing in about 8.5 ms, but it computes `e^3` rather than `e` (issue #16),
while the new code computes the correct value by direct exponentiation. Fair
ratios, unfair clock.

Closing that gap needs a fast final exponentiation chain and cyclotomic
squaring. Both are well understood; neither is written.

### Bearing on issue #16

The hard exponent's base-p digits are full width (381, 254, 381, 126 bits for
BLS12-381), so the cheap decomposition that would give `e` directly does not
exist. The factor of 3 the standard BLS12 chain introduces is a property of that
algorithm, not a mistake in this codebase. That answers the question the issue
put to the maintainer for the BLS12 half. The BN factor of roughly `12*X^3` is
still unexplained and still looks unintentional.

### Remaining in Phase 3

- [x] Fast final exponentiation chain — done, see part 4 below
- [ ] Cyclotomic squaring, for a further cut in the final exponentiation
- [ ] BN support in the new Miller loop
- [ ] Curve context struct, retiring the globals and with them A3 and A4
- [ ] Retire the legacy layer once the new one covers every entry point


---

## Phase 3, part 4 — fast final exponentiation. Gate met on the clock.

**Status: the >= 3x gate is met on an absolute, like-for-like wall-clock
comparison.**

### Measured, BLS12-461, both sides computing e^3

| | legacy | new | speedup |
|---|---|---|---|
| Miller loop | 2576.8 us | 858.4 us | **3.00x** |
| final exponentiation | 6065.2 us | 1588.5 us | **3.82x** |
| **whole pairing** | **8642.0 us** | **2446.9 us** | **3.53x** |

For reference, the exact-`e` path costs 9861.8 us, which is why the fast chain
matters rather than being a nicety.

### The chain

Rests on an identity verified numerically on both BLS12 curves:

```
3*lambda = (x-1)^2 * (x+p) * (x^2+p^2-1) + 3
```

Five parameter exponentiations replace a 1534-bit ladder. Every inverse in the
hard part is a conjugation, because the easy part has already put the element in
the cyclotomic subgroup.

### This settles the BLS12 half of issue #16

The factor of three is **a property of the standard algorithm, not a defect**.
`lambda` alone has no short evaluation — its base-p digits are full width — and
recovering `e` from `e^3` needs a cube root in mu_r, which is a full-width
exponentiation. The legacy library and RELIC both land in the same place.

Since `gcd(3, r) = 1`, `e^3` is still bilinear and non-degenerate, so any
protocol that only compares pairings is unaffected. Raw values will not match an
implementation that emits `e`.

Both paths are kept and documented: `pairing_final_exp_fast` for speed,
`pairing_final_exp_plain` when the exact value is needed. The suite pins the
relationship between them, so neither can drift.

**The BN factor of roughly `12*X^3` remains unexplained and still looks
unintentional.** Nothing here changes that half of the issue.

### Phase 3 scorecard against its gate

| Gate condition | Result |
|---|---|
| Vectors pass on all three curves | yes, 857 each, new layer |
| Differential test, 10^6 inputs per level | yes, 1,204,813 checks per curve |
| Sanitizers clean | yes |
| **Speedup >= 3x** | **3.53x on the whole pairing** |

BN landed shortly after; see part 5.


---

## Phase 3, part 5 — BN support

**Status: all three curves now have a working, verified pairing on the new
stack.**

BN's optimal ate needs two correction lines after the loop, and those need the
skew Frobenius on the twist:

```
psi(x, y) = (conj(x) * gamma^-2, conj(y) * gamma^-3)
```

with `psi^2` multiplying by the Fp2 norms of those. The generator asserts their
imaginary parts are zero, which is what the derivation predicts — a cheap
independent check that it was done right.

The new BN pairing matched the reference coefficient by coefficient on the first
run, as BLS12 had.

| | legacy | new | speedup |
|---|---|---|---|
| BN-462 Miller loop | 4077.9 us | 1326.8 us | **3.07x** |

### No fast final exponentiation for BN, on purpose

The identity the BLS12 chain rests on is specific to that family. Applying it to
BN gives a wrong exponent, which the suite caught the moment BN was wired in.
BN has its own standard chain and it is simply not written yet, so BN callers
get the exact slow path or nothing.

Declining to provide one is a deliberate choice, not an oversight. The legacy
library ships a BN "optimal" final exponentiation that raises to roughly
`12*X^3` times the correct exponent, and nothing noticed for years because the
result is still bilinear. **A missing function fails better than a plausible
wrong one.**

### Current state of the new stack

| Curve | field | curve | Miller | final exp (exact) | final exp (fast) |
|---|---|---|---|---|---|
| BLS12-381 | yes | yes | yes | yes | yes |
| BLS12-461 | yes | yes | yes | yes | yes |
| BN-462 | yes | yes | yes | yes | **not written** |

17 CTest targets green, ASan and UBSan clean, zero warnings on a clean build.

### Remaining in Phase 3

- [ ] BN fast final exponentiation chain
- [ ] Cyclotomic squaring, a further cut for both families
- [ ] Retire the legacy layer, which is what actually deletes the globals and
      defects A3 and A4. The new layer never had them: its constants are
      compile-time, so the planned "curve context struct" turned out to be
      unnecessary rather than merely deferred.
