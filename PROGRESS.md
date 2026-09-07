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

### BN fast chain: initially withheld, then written properly

The BLS12 identity is family-specific; applying it to BN gave a wrong exponent,
which the suite caught the moment BN was wired in. Rather than ship a plausible
wrong chain — the exact defect the legacy library has — BN was left with the
exact slow path until the right chain existed. It now does; see part 6.

### Current state of the new stack

| Curve | field | curve | Miller | final exp (exact) | final exp (fast) |
|---|---|---|---|---|---|
| BLS12-381 | yes | yes | yes | yes | yes |
| BLS12-461 | yes | yes | yes | yes | yes |
| BN-462 | yes | yes | yes | yes | yes (exact e) |

17 CTest targets green, ASan and UBSan clean, zero warnings on a clean build.

### Remaining in Phase 3

- [x] BN fast final exponentiation chain — done, part 6
- [ ] Cyclotomic squaring, a further cut for both families
- [ ] Retire the legacy layer, which is what actually deletes the globals and
      defects A3 and A4. The new layer never had them: its constants are
      compile-time, so the planned "curve context struct" turned out to be
      unnecessary rather than merely deferred.


---

## Phase 3, part 6 — BN fast final exponentiation

Derived symbolically rather than from memory:

```
lambda = d0 + d1 p + d2 p^2 + d3 p^3
d0 = -36x^3 - 30x^2 - 18x - 2
d1 = -36x^3 - 18x^2 - 12x + 1
d2 =            6x^2      + 1
d3 =                        1
```

Regrouped by power of x, that is three parameter exponentiations plus a few
shared small powers. **Unlike BLS12 the decomposition is exact**, so BN returns
`e`, not `e^3`.

| BN-462 | legacy | new | speedup |
|---|---|---|---|
| Miller loop | 4013.5 us | 1286.2 us | 3.12x |
| final exponentiation | 5457.5 us | 1499.9 us | 3.64x |
| **whole pairing** | **9471 us** | **2786 us** | **3.40x** |

The final exponentiation row flatters the legacy side, which is computing the
wrong exponent while the new one computes `e` exactly. The new path is both
faster and correct.

### A bug the suite caught, worth recording

The first version of this chain used `ELIPS_LOOP` to compute `f^x`. Those
constants coincide on BLS12, but **the BN Miller loop runs over `6x+2` while its
final exponentiation needs `x`**, so BN silently got a wrong exponent. The
generator now emits `ELIPS_PARAM` for the mother parameter separately from
`ELIPS_LOOP` for the Miller loop, and the comment at the use site says why they
must not be confused.

This is the same shape as the defect the whole project started from: a
wrong-but-bilinear pairing that passes casual inspection. The difference is that
this time a test caught it within minutes.

### Where Phase 3 stands

| Curve | whole pairing, legacy | new | speedup | new value |
|---|---|---|---|---|
| BLS12-461 | 8642 us | 2447 us | 3.53x | `e^3` |
| BN-462 | 9471 us | 2786 us | 3.40x | `e` (exact) |


---

## Phase 3, part 7 — public API, generators, standalone

**Status: the new layer is self-contained. Phase 3's arithmetic goals are met.**

### What changed

Verified group generators for all three curves, **generated rather than
transcribed**. The script asserts each is on its curve and of order exactly `r`
before emitting it, so a generator that quietly landed in the wrong subgroup
cannot ship.

On top of that: `ep_generator`, `ep2_generator`, `ep_in_subgroup`,
`ep2_in_subgroup`, and `elips_pairing`, which validates both inputs and returns
0 rather than a subtly wrong value when a point is the identity or outside the
order-`r` subgroup.

**The subgroup check is new to this codebase.** Nothing in the legacy layer ever
performed one, which is precisely what small-subgroup attacks rely on.

### Standalone

`test/standalone_test.c` links no legacy code. Twelve checks per curve, all
passing on all three:

```
generators on curve and of order r; pairing non-degenerate and in mu_r;
bilinearity in both arguments; the identity rejected rather than paired
```

**BLS12-381 has a working pairing here for the first time** — the legacy layer
never supported it. This is also what unblocks retiring that layer, since tests
no longer need to borrow points from it.

20 CTest targets green on both build types, ASan and UBSan clean, zero warnings.

### Phase 3 final scorecard

| Gate condition | Result |
|---|---|
| Vectors pass, all three curves | yes, 857 each against the new layer |
| Differential test, 10^6 inputs per level | yes, 1,204,813 checks per curve |
| Sanitizers clean | yes |
| Speedup >= 3x | **3.53x BLS12-461, 3.40x BN-462** |

| Curve | pairing, legacy | new | speedup | value |
|---|---|---|---|---|
| BLS12-381 | not supported | works | — | `e^3` |
| BLS12-461 | 8642 us | 2447 us | 3.53x | `e^3` |
| BN-462 | 9471 us | 2786 us | 3.40x | `e` exact |

### What is left, and where it belongs

- **Cyclotomic squaring** — this is a Phase 4 item in the plan, not Phase 3.
  Worth roughly a further 1.4x on the pairing.
- **Retiring the legacy layer** — now unblocked. It is what finally deletes the
  globals and defects A3 and A4. The planned "curve context struct" turned out
  to be unnecessary: the new layer's constants are compile-time, so it never had
  globals to remove.
- **Windowing and GLV** for scalar multiplication — Phase 4.

---

## Phase 4 — optimization and the RELIC comparison (issue #5)

Three of this phase's six items had already landed in Phase 3 (projective
coordinates, sparse multiplication). This round added the rest bar GLV.

### Cyclotomic squaring and a dedicated fp6 squaring

`fp6_sqr` no longer delegates to `fp6_mul`. Taking each doubled cross term from
a squaring via `2ab = (a+b)^2 - a^2 - b^2` costs six fp2 squarings where the
multiplication path costs six fp2 multiplications.

`fp12_sqr_cyc` uses the defining property of the cyclotomic subgroup. There
`conj(a) = a^-1`, so `d0^2 - v*d1^2 = 1`, and the `d0^2 + v*d1^2` that squaring
needs rewrites as `2*d0^2 - 1`:

```
a^2 = (2*d0^2 - 1) + (2*d0*d1) w
```

Derived here rather than transcribed, so the tower's own conventions apply and
there is no index mapping to get wrong. Valid only inside the subgroup, and the
test asserts that the easy part really lands there rather than assuming it.

| whole pairing | before | after |
|---|---|---|
| BLS12-461 | 2447 us | **2014 us** (3.91x vs legacy) |
| BN-462 | 2786 us | **2450 us** (3.51x vs legacy) |

### Constant-time fixed-window scalar multiplication

Replaces double-and-add-always. A 4-bit window does 256 doublings and 64
additions on a 256-bit scalar where the old routine did 256 of each. The table
lookup scans all 16 entries under a mask, so no address depends on the scalar.

### Head to head with RELIC, BLS12-381

Same machine, same compiler, RELIC built from its `gmp-pbc-bls381` preset.

| | RELIC | ELiPS | ratio |
|---|---|---|---|
| pairing | 955.5 us | **1145.2 us** | 1.20x slower |
| G1 scalar mult | 128.1 us | 181.1 us | 1.41x slower |
| G2 scalar mult | 204.5 us | 467.1 us | 2.28x slower |

Our breakdown: Miller loop 441.6 us, final exponentiation 719.1 us.

**Three caveats, all of which matter for reading those numbers.**

1. RELIC is using `ARITH=gmp`, its portable backend. Its fastest configurations
   are x86-64 assembly (`x64-asm-6l` for this curve) and cannot run on this
   AArch64 machine. So this is portable-C against portable-C, and on x86-64
   with assembly RELIC would pull further ahead. That gap is what Phase 6 exists
   to close.
2. RELIC's `g1_mul`/`g2_mul` defaults are not necessarily constant time and
   likely use GLV endomorphisms and wNAF. Ours is constant-time fixed-window.
   Some of the scalar multiplication gap is that choice, not implementation
   quality — but not all of it, and G2 at 2.28x is the clear weak point.
3. RELIC's BLS12 pairing uses the same standard chain, so it also returns `e^3`.
   Like for like on that point.

**Reading it plainly:** within 20% of RELIC on the pairing, from a library that
was roughly nine times slower than that and did not support BLS12-381 at all.
"Comparable to RELIC" is a fair description on this hardware and configuration.
It is not yet true against RELIC's assembly builds.

### Remaining in Phase 4

- [ ] GLV / split-scalar decomposition, constant-time. The obvious next target:
      G2 scalar multiplication is the worst gap at 2.28x, and GLV with the
      psi endomorphism is exactly what closes it.
- [ ] Compressed squaring for the `f^x` chains
- [ ] Revisit the >= 8x gate. Against the legacy layer we are at 3.9x; the gate
      as written is not met, and the RELIC comparison suggests the remaining
      headroom on portable C is modest. Worth deciding whether 8x was the right
      target or whether "within X% of RELIC" is the more meaningful gate.

### GLV on G2, and the revised gate

The Phase 4 gate is now **RELIC-relative** rather than a multiple of the legacy
layer: BLS12-381 pairing within 25% of RELIC on comparable arithmetic, scalar
multiplication within 1.5x. The old 8x figure was set before anyone knew where
the time went, and it measured the old code's weaknesses rather than this one's
quality.

`psi` acts on G2 as multiplication by `p mod r` — verified on all three curves,
not assumed. For BLS12 that reduces to the mother parameter `x`, only 64 to 77
bits against a 255 to 308 bit order, so a base-`|x|` decomposition gives four
short digits and quarters the ladder.

| G2 scalar mult | fixed window | GLV | speedup |
|---|---|---|---|
| BLS12-381 | 461.5 us | **312.6 us** | 1.48x |
| BLS12-461 | 903.6 us | **603.0 us** | 1.50x |

### Standing against the gate, BLS12-381

| | RELIC | ELiPS | ratio | gate | met |
|---|---|---|---|---|---|
| pairing | 955.5 us | 1145.2 us | **1.20x** | <= 1.25x | **yes** |
| G1 scalar mult | 128.1 us | 181.1 us | **1.41x** | <= 1.5x | **yes** |
| G2 scalar mult | 204.5 us | 312.6 us | **1.53x** | <= 1.5x | marginally no |

### Two honest limits on the GLV

**The decomposition is variable time.** It uses GMP division; the ladder around
it is constant time. So `ep2_mul_glv` is opt-in and documented as unsafe for
secret scalars, and `ep2_mul` remains the default. A constant-time Barrett
decomposition is the remaining piece.

**BN gets no GLV.** There `psi` acts as `6x^2`, 231 bits against a 462-bit
order, so only a two-dimensional split exists and the win would be about half.
Worth doing; not done.

### The remaining lever

G2 is 0.03x outside its gate, and the cause is identifiable: the complete
addition computes a full doubling on every call so the coincident-point case can
be selected without branching. That is roughly a 1.4x tax on every addition, and
the GLV ladder is addition-dominated.

The principled fix is Renes–Costello–Batina complete formulas in homogeneous
projective coordinates, which are genuinely complete in one formula at about 12
multiplications, against the current 16 plus a wasted doubling. That would speed
up the ladder and possibly the Miller loop. It is a coordinate-system change for
the EC layer, so it is a real piece of work rather than a tweak.

### RCB complete formulas, and the constant-time GLV

The curve layer moved from Jacobian to homogeneous projective with the
Renes–Costello–Batina complete formulas (a = 0), and the Miller loop's line
functions were re-derived to match. **All 20 tests passed on the first run after
a coordinate-system change** — the Phase 0 vectors earning their keep.

It is a tradeoff, not a clean win.

| | Jacobian | RCB |
|---|---|---|
| ep2_add (381) | ~4.00 us | **2.24 us** |
| ep2_dbl (381) | ~0.80 us | 1.44 us |
| pairing (381) | 1145.2 us | 1175.1 us |
| G2 via GLV (381) | 312.6 us | **283.2 us** |
| G2 window ladder (381) | 461.5 us | 552.7 us |

Addition is 1.8x cheaper and doubling 1.8x dearer, because Jacobian
`dbl-2009-l` at 2M+5S is simply cheaper than any homogeneous doubling. So
addition-heavy work (the GLV ladder) wins and doubling-heavy work (the pairing,
the plain window ladder) loses a little.

Kept on three grounds: it is what carries G2 over the gate, the pairing cost is
2–5% and still inside its gate, and **the incomplete `add_generic` and its
precondition are gone** — one formula is now correct for every input pair, so
there is no longer a routine in the codebase that is wrong for inputs a caller
might plausibly supply.

One measurement worth keeping: splitting the line evaluation out of the doubling
cost 11% of the Miller loop, because the shared squares are most of the work.
Re-fusing them recovered it.

The GLV decomposition is now constant time — restoring division by mask rather
than GMP's `mpz_tdiv_qr` — for 1.6% overhead, verified against the plain ladder
on 3000 scalars per curve. It is no longer opt-in.

**A precondition that was implicit and should not have been:** `psi` acts as
multiplication by `x` only on G2. On an arbitrary twist point it does not, so
`ep2_mul_glv` returns a wrong answer rather than failing. That is why `ep2_mul`
does not dispatch to it, and why the known-answer vectors keep using the
general routine.

### Phase 4 standing against the RELIC-relative gate, BLS12-381

| | RELIC | ELiPS | ratio | gate | met |
|---|---|---|---|---|---|
| pairing | 955.5 us | 1175.1 us | **1.23x** | <= 1.25x | **yes** |
| G1 scalar mult | 128.1 us | 176.2 us | **1.38x** | <= 1.5x | **yes** |
| G2 scalar mult (GLV) | 204.5 us | 287.9 us | **1.41x** | <= 1.5x | **yes** |

All three met, with G2 now constant time on the fast path.

Caveat that stays true: RELIC is on its portable GMP backend because its x86-64
assembly cannot run on this AArch64 machine. This is portable C against portable
C. Phase 6 is what addresses the rest.

---

## Phase 5 — Constant-time completion and API hardening (issue #6)

**Status: complete. Exit gate met.**

Phases 3 and 4 had already delivered constant-time field arithmetic, scalar
multiplication, inversion and subgroup checks, so this phase was the three
things the hand-off named — CSPRNG, dudect, serialization — plus the trust
boundary that ties them together.

### What exists now

| Artifact | Purpose |
|---|---|
| `include/elips/sysrand.h`, `src/util/sysrand.c` | The system CSPRNG, with no curve dependency so both layers can use it |
| `include/elips/random.h`, `src/util/random.c` | Uniform field elements and scalars, and the scalar range predicate |
| `include/elips/serialize.h`, `src/arith/serialize.c` | Compressed and uncompressed point encodings, with validating readers |
| `fp_sqrt`, `fp2_sqrt`, `fp_exp`, `fp2_exp` | What decompression rests on |
| `test/serialize_test.c` | Round trips, the format's fixed points, and every rejection |
| `test/dudect_test.c` | Timing-leakage tests, with a negative control |

### The CSPRNG

`elips_random_bytes` reads the operating system: `getrandom(2)` on Linux,
`arc4random_buf` on macOS and the BSDs, `/dev/urandom` as the fallback. No
userspace generator, no seeding, no state. Every failure mode of the old design
— guessable seed, two processes agreeing, a stale state after fork — comes from
having state to get wrong.

`fp_rand` and `elips_random_scalar` reduce 128 extra bits into range with
`mpn_sec_div_r`, so the bias is below 2^-128 and the running time does not
depend on the sample. Rejection sampling would be exactly uniform and would
leak.

The legacy layer's `gmp_randstate_t` now takes a 256-bit seed from the same
source instead of `time(NULL)`. **Stated plainly: this fixes the guessable
seed, not the generator.** GMP's Mersenne Twister is still predictable from its
own output. Nothing on that path produces key material — it generates test
points — and issue #17 retires it.

### Serialization, and what it took to make it interoperate

Format is the BLS12-381 convention: big-endian, three flag bits in the top of
the first byte, Fp2 written imaginary part first. The field width is derived
from the requirement rather than from `ceil(FP_BITS/8)`, which is what lets one
implementation serve all three curves:

| Curve | FP_BITS | naive width | spare bits | width used |
|---|---|---|---|---|
| BLS12-381 | 381 | 48 | 3 | 48 (the standard's own size) |
| BLS12-461 | 461 | 58 | 3 | 58 |
| BN-462 | 462 | 58 | **2** | **59** |

BN-462 has only two spare bits at 58 bytes, so the three-flag scheme does not
fit and it takes 59. Choosing `(FP_BITS + 3 + 7) / 8` makes that fall out
rather than needing a special case.

**Byte-exact against the specification, which took a parameter change.** The
G1 generator already matched the published encoding. The G2 generator did not:
`gen_params.py` searched for the first suitable point, and found a perfectly
valid generator of the same subgroup that no other implementation uses. Format
compatibility would not have saved it — signatures verify against the standard
generator or not at all. `gen_params.py` now takes BLS12-381's generators from
the specification, in the published encoded form, and re-derives and re-asserts
them (on curve, order exactly r) before emitting. The two curves with no
specification keep the search.

The three published encodings are pinned as known answers in
`serialize_test.c`. That is the only place in this project where the new layer
is checked against something outside it.

**Readers validate.** Non-canonical coordinates, inconsistent flags, points off
the curve and points outside the order-r subgroup each get their own return
code, and each has its own negative test — asserting the specific code, so a
decoder cannot pass by rejecting everything. The subgroup check is the one that
stops small-subgroup attacks and the one an implementation is most likely to
skip, because everything appears to work without it.

### dudect

Two input classes, fixed and random, timed and compared with Welch's t-test
over a ladder of percentile crops; the reported figure is the largest |t|.
Following dudect's own thresholds: under 5 is clean, over 10 is leaking, and
between the two the sample is too small to decide. A first reading over 10 is
re-measured on a fresh sample before it is reported, because a preempted run on
a shared machine throws a large t with nothing wrong.

```
  [PASS ] fp_mul           max|t| =    2.41
  [PASS ] fp_add           max|t| =    2.22
  [PASS ] fp_cselect       max|t| =    1.08
  [PASS ] fp_inv           max|t| =    1.08
  [PASS ] ep_mul           max|t| =    0.93
  [PASS ] ep2_mul          max|t| =    2.27
  [PASS ] ep2_mul_glv      max|t| =    1.31
  [PASS ] miller           max|t| =    2.68
  [PASS ] control_vartime  max|t| =  779.34   <- negative control, must leak
```

**The negative control is the point.** `control_vartime` calls
`fp_inv_vartime` — documented as variable time, never used on secrets — on
secret inputs, and the test fails if no leak is found. A leakage detector that
cannot detect a leak says nothing about the targets it passes. Same argument as
`kat.detects_corruption` in Phase 0.

Registered in CTest for Release builds only. Under Asan or Ubsan the
measurement is of the instrumentation, not of the routine, so a clean result
there would be evidence of nothing. The binary is still built in every
configuration so the harness itself is sanitized.

### A measurement that changed the harness, worth recording

The first version reported **|t| = 98 for `fp_add`** — an addition with no
branch in it. The leak was in the test, not the code: preparing a random field
element calls into the kernel, and doing that inside the timed loop for one
class and not the other leaves the cache and branch predictor in visibly
different states. Preparing every input before any measurement made the same
test read 1.4.

That is the failure mode of this whole technique. A timing test that measures
its own setup will happily accuse correct code, and the accusation looks exactly
like a real finding.

### One real leak found and fixed

`fp_inv` early-returned on a zero input. Correct, and a branch on the operand —
the one thing that routine is not allowed to do. dudect did not catch it
because random inputs are never zero. It now runs `mpn_sec_invert`
unconditionally and selects, and `fp_difftest` pins `fp_inv(0) == 0`.

### Cleared out while here

- **`ep_add_generic` / `ep2_add_generic` deleted.** Since the move to RCB they
  were aliases for the complete addition, but the header still documented them
  as "wrong for p==q, p==-q or infinity". A public routine advertised as unsafe
  for inputs a caller can plausibly supply is a defect waiting for its first
  careless call site.
- **Two leaks in the legacy 4-split routines** (`bls12_4split_G2_scm`,
  `bls12_4split_G3_exp`): `A`/`B` and `C`/`D` were initialised beside `x_1` and
  `x_2` and never released. Phase 2's "zero leaks" was measured with macOS
  `leaks`; LeakSanitizer on Linux sees them.
- **Test-harness leaks** in `kat_runner_new.c`, `pairing_test.c` and
  `ec_test.c`, for the same reason.

### Corrections to earlier claims

**"Zero warnings, ASan and UBSan clean" was true on macOS and false on Linux.**
Both were measured with clang on Apple Silicon. On this branch, before any
Phase 5 change:

- `test/kat_runner.c` and `test/kat_runner_new.c` used `getline`, `strtok_r`
  and `ssize_t` without a feature-test macro. Under `-std=c11` glibc hides all
  three. The new-layer runner failed to compile; the legacy one compiled with
  `getline` and `strtok_r` implicitly declared as returning `int`, which
  truncates their pointer results on any 64-bit target. **The Linux CI job
  cannot have been passing.**
- GCC reported six `-Wstringop-overflow` warnings in `fp12_mul_sparse035`. A
  false positive: given an `fp12_t` parameter, GCC narrows what it believes
  `f[1]` to be as soon as the body indexes into it, then reports every later
  whole-`fp6` access as overflowing. The object really is 288 bytes and the
  code was correct, but a false warning that cannot be told apart from a true
  one is worth a signature change, so the routine now takes its two `fp6`
  halves separately. No runtime cost.
- LeakSanitizer runs by default alongside AddressSanitizer on Linux and is not
  supported on Apple Silicon, which is why the leaks above were invisible.

All three are fixed. Both sanitizer builds are green on Linux now.

### Measured

x86-64, gcc -O2, this machine. Not comparable with the Apple Silicon figures
earlier in this file — the same Miller loop measures 1474 us here against 442
us there — so these are for internal proportions only.

| | BLS12-381 | BLS12-461 | BN-462 |
|---|---|---|---|
| Miller loop | 1474 us | 2702 us | 4075 us |
| final exponentiation (fast) | 2384 us | 3833 us | 3713 us |
| `elips_pairing`, with subgroup checks | 6033 us | 10782 us | 14157 us |
| G1 decompress and validate | 550 us | 1070 us | 1571 us |
| G2 decompress and validate | 2243 us | 3777 us | 5321 us |

**The subgroup checks cost more than a third of `elips_pairing`** — 2175 us of
6033 on BLS12-381 — because each one is a full scalar multiplication by r.
That is the price of the guarantee, and it is being paid on every call.

Deliberately not optimised. The known fast tests (Scott's `psi(Q) == [x]Q` for
G2, the GLV-endomorphism test for G1) are each valid only under conditions on
the curve that have to be checked, not recalled, and a subgroup test that is
wrong is a security hole rather than a slow path. Correct and slow now;
Phase 6 or the backlog can make it fast, with the derivation written down.

### Phase 5 scorecard against its gate

| Gate condition | Result |
|---|---|
| dudect reports no leakage on secret-dependent paths | yes, 8 targets, max abs t 2.68 |
| ...on both architectures | CTest targets, so both CI legs run them |
| the detector can detect | yes, negative control at abs t 779 |
| variable-time paths only behind `_vartime` | yes; `fp_inv_vartime` is the only one, and it is the control |
| malformed and off-curve inputs rejected | yes, six distinct failure modes, each with its own test |

31 CTest targets green on Release, 23 on Asan and on Ubsan (dudect is not
registered under sanitizers), zero warnings on gcc and clang.

### What Phase 5 did not do

- **The new headers are not installed.** `include/elips/*.h` and the per-curve
  `elips_arith_*` libraries are outside the install and export set, so the new
  API cannot be consumed via `find_package`. That needs a decision on how a
  downstream project selects its curve, which is a packaging question rather
  than a hardening one, so it is left for whoever makes it.
- **No hash-to-curve.** Serialization is the half of interoperability this
  phase was asked for; a full BLS signature implementation also needs
  `hash_to_curve`, which is a specification of its own.
- **Scalar range is offered, not enforced.** `elips_scalar_is_reduced` exists
  and is tested, but `ep_mul` cannot reject k >= r: the subgroup checks
  legitimately call it with exactly r. Callers with secret scalars reduce
  first; documented at the declaration.

---

## Phase 5b — hash to curve, curve selection, and the decisions on what is left

**Status: complete. RFC 9380 is implemented on all three curves and BLS12-381
reproduces the specification's published test vectors.**

### What exists now

| Artifact | Purpose |
|---|---|
| `include/elips/hash_to_curve.h`, `src/hash/hash_to_curve.c` | RFC 9380 `hash_to_curve` and `encode_to_curve` for G1 and G2 |
| `src/hash/h2c_tmpl.h` | The two maps, written once and instantiated per group |
| `include/elips/sha256.h`, `src/hash/sha256.c` | SHA-256, so the build stays one library plus GMP |
| `tools/reference/h2c_ref.py` | An independent RFC 9380 implementation: the oracle |
| `tools/reference/h2c_iso_bls12_381.json` | The isogeny data, with provenance |
| `tools/reference/gen_h2c_params.py`, `gen_h2c_vectors.py` | Constants and vectors, both gated on the reference self-test |
| `test/kat/h2c_*.vec`, `test/h2c_test.c` | 46 records per curve, plus a corruption control |

### Two maps, and why

| Curve | Map | Status |
|---|---|---|
| BLS12-381 | simplified SWU over an 11-isogenous (G1) / 3-isogenous (G2) curve | The registered RFC 9380 suite. **Byte-exact with the specification's vectors.** |
| BLS12-461, BN-462 | Shallue–van de Woestijne | No registered suite exists; SvdW is what RFC 9380 uses for BN254 and needs no isogeny. |

Both curves have `A = 0`, which simplified SWU cannot handle, hence the isogeny
on BLS12-381 and SvdW elsewhere. The plan's §10.1 records the alternatives that
were rejected and why.

### The one place data had to be trusted, and what was done about it

The 11- and 3-isogeny coefficients are the only numbers in this project that
cannot be derived from the curve parameters: finding an isogenous curve needs
Vélu's formulas or modular polynomials, and the result is only useful if it is
*the same model the specification chose*.

They were extracted mechanically from py_ecc 8.0.0 — never retyped — committed
as data with a provenance note, and then **verified three independent ways**
before the generator will emit them:

1. points of `E'` map onto `E`;
2. the map is a group homomorphism, so it really is an isogeny;
3. the whole suite reproduces RFC 9380's published vectors.

A mistranscribed coefficient fails all three.

Everything else is derived: the SvdW `Z` values come from the Appendix H.1
search, `c1..c4` from their defining equations (each asserted), and the cofactor
multipliers from the curve parameter — `1-x` on BLS12 G1 and `3(x²-1)h2` on G2,
both landing on RFC 9380's own values without being told them.

### Verification

```
BLS12-381   46 passed, 0 failed     <- includes RFC 9380 J.9.1 and J.10.1
BLS12-461   46 passed, 0 failed
BN-462      46 passed, 0 failed
```

Each layer is pinned separately — SHA-256, `expand_message_xmd`,
`hash_to_field`, then the map — because they all fail the same way, with a wrong
point, and one end-to-end check would not say which broke. There is a corruption
control, as for the Phase 0 vectors.

The reference's own self-test also passes `expand_message_xmd` against RFC 9380
Appendix K.1, and checks that every suite's output lands in the order-`r`
subgroup for random messages.

### A correction the search itself produced

`cofactor_multiplier` first asserted that the multiplier is a multiple of the
cofactor. It is not: on BLS12-381 G1, `h_eff = 1-x` is 64 bits against a 128-bit
cofactor, and multiplying by it still lands in the subgroup because `E(Fp)` is
not cyclic there and `(1-x)` kills both of its cofactor components. That is
precisely why RFC 9380 can specify such a short multiplier. The assertion would
have rejected the correct answer; it now checks the property that matters —
`[h]P` has order dividing `r` for points drawn without reference to the subgroup.

### Fast G2 cofactor clearing

`h_eff` is 636 bits on BLS12-381, so cofactor clearing, not the map, dominated
`hash_to_g2`. Budroni–Pintore (ePrint 2017/419):

```
[h_eff]Q = [x^2 - x - 1]Q + [x - 1]psi(Q) + psi^2([2]Q)
```

| hash_to_g2 | before | after |
|---|---|---|
| BLS12-381 | 4070 us | **2133 us** (1.91x) |
| BLS12-461 | 7041 us | **3794 us** (1.86x) |

The identity was **verified numerically against the derived `h_eff` on random
points of the twist** — not of G2, where much weaker relations hold and would
hide a wrong chain — and then again by the RFC vectors still matching. BN keeps
the plain multiplication; no chain for it was verified.

### Measured

| | BLS12-381 | BLS12-461 | BN-462 |
|---|---|---|---|
| `hash_to_g1` | 469 us | 696 us | 563 us |
| `encode_to_g1` | 271 us | 417 us | 290 us |
| `hash_to_g2` | 2133 us | 3794 us | 4987 us |

x86-64, gcc -O2, the same slow VM as the Phase 5 numbers. G2 is still dominated
by the two short ladders rebuilding their window tables; a shared table would
help and is not done.

### Curve selection, and an API that can actually be installed

```bash
cmake -B build -DELIPS_CURVE=BLS12_381    # default; also BLS12_461, BN_462
```

The choice decides which library is installed and what `find_package(ELiPS)`
hands back as `ELiPS::arith`, with `ELiPS_CURVE` in the package config so a
consumer can read it back. An unknown value is a configure error listing the
valid ones. The test suite still builds all three.

This also closes the gap Phase 5 recorded as remaining: `include/elips/*.h` was
not installed at all, so the API carrying the pairing, the serialization and now
the hash-to-curve could not be consumed by anyone. CI now builds a downstream
project against the installed package with `-Wall -Wextra -Werror`, hashes to
G1, pairs, serializes and round-trips — and configures and installs all three
curve selections.

### Constant time

`hash_to_g1` and `hash_to_g2` are dudect targets, because the message is not
always public: an OPRF or a PAKE hashes a secret. Both branches of every
value-dependent choice are computed and selected under a mask, square roots
included.

```
  [PASS ] hash_to_g1       max|t| =     1.93
  [PASS ] hash_to_g2       max|t| =     1.67
  [PASS ] control_vartime  max|t| =  1030.35   <- negative control, must leak
```

### State of the art, surveyed and decided

`MODERNIZATION_PLAN.md` §10 now records a decision for every item Phases 0–5
left open, with citations: fast subgroup membership tests (Scott 2021/1130 as
corrected by El Housni–Guillevic–Piellard 2022/352, generalised by Dai et al.
2022/348), GLV for BN and for G1, Karabina compressed squaring, multi-pairing
with fixed-argument precomputation, safegcd, and why these three curves.

**Two of those changed a Phase 6 assumption.**

- **The x86-64 assembly target is wrong in §6.** The state of the art there is
  AVX-512 IFMA applied to the extension fields, not a scalar limb multiply:
  ePrint 2025/1283 (TCHES 2025) reports ~1.2M cycles for a full BLS12-381
  pairing, about 400 µs at 3 GHz against this library's 1175 µs of portable C.
  So the headroom is roughly **3x, not the ~1.3x the plan assumed**. AArch64
  keeps the scalar plan; there is no IFMA equivalent.
- **Fast subgroup tests are the largest single win available** — the checks
  inside `elips_pairing` cost 2175 µs of 6033 µs — and are still not
  implemented, deliberately. See below.

### Deliberately not done, and the reason is the same each time

**Fast subgroup membership tests.** Worth about a third of the pairing. Not
built because their correctness conditions are exactly the kind of thing that is
easy to get wrong, the literature has already had to correct one of those proofs
once, and the full texts were not reachable from this session (the egress policy
blocks `eprint.iacr.org`). The distinction that decided it:

> A wrong cofactor chain computes a different multiple, so the RFC vectors stop
> matching and it cannot ship. A wrong subgroup test accepts points it should
> reject, and no vector notices.

The first was built this session; the second was not. The Phase 6 task is to
have `gen_params.py` derive and assert the conditions per curve — the way it
already asserts that generators have order `r` — and only then switch it on.

### Two Linux-only findings, again

Not new defects, but worth recording that both surfaced only off macOS:

- `tools/reference/gen_h2c_vectors.py` takes about a minute because the
  reference does 462-bit scalar multiplications in pure Python. Memoising the
  suites cut it from "times out" to 66 s.
- The `-Wstringop-overflow` lesson from Phase 5 applied again while writing the
  Horner evaluation over the isogeny tables: a table of `limb_t[n][2][FP_LIMBS]`
  cannot be indexed generically by casting to the field type, because those are
  array types. The template takes an `H2C_ELEM` macro from its includer instead.

---

## Phase 5c — edge cases, examples, and a real warning gate

**Status: complete. One latent bug found and fixed; CI now gates what it only
claimed before.**

### The CI pipeline was already green

Worth stating plainly, since it was the starting question: the last two runs on
this branch passed all 8 jobs, both architectures, all three build types. The
one failing run is `206b08a`, the commit before Phase 5 — which is the run that
was broken by the `getline`/`ssize_t` portability defect Phase 5 fixed.

What was wrong was not that CI failed, but what it did **not** check:

| Gap | Fix |
|---|---|
| `ELIPS_WERROR` existed but no job set it, so "zero warnings" was a claim | CI configures with `-DELIPS_WERROR=ON`, and the option now reaches the per-curve libraries and every test target, not just the legacy library |
| The downstream consumer returned 0/1 with no output, so a regression that still exited 0 passed silently | It now asserts thirteen things and prints each: hash-to-curve, both subgroup checks, bilinearity through the installed package, a serialization round trip, and a corrupted encoding being refused |
| The matrix only ever exercised the default curve | The curve-selection job builds, installs and runs the examples for all three |
| `actions/checkout@v4` was being force-upgraded to Node 24 with a deprecation warning | bumped to `@v5` in both workflows |

### One real bug, found by looking for it

`fp2_sqrt` did not tolerate its output aliasing its input:

```c
fp2_cselect(r, cand_i, cand_b, mask);   /* writes r */
fp2_sqr(chk, r);
int ok = fp2_eq(chk, a);                /* reads a -- the same storage */
```

With `r == a` the verification compared the root against itself squared, decided
the input was not a residue, and returned zero. `fp_sqrt` has the same shape but
reads `a` before writing `r`, so it was fine; `fp2_sqrt` was not.

Nothing in the library aliases these, which is why every existing suite passed.
The header promises aliasing is tolerated, so it was a contract violation and a
trap for the next caller. It now accumulates into a local and writes `r` once,
at the end.

It first appeared on BLS12-461 only, because whether it shows depends on the
random draw — which is a good argument for the test using real randomness rather
than a fixed seed.

### test/edge_test.c

84 checks per curve, on all three, in the places implementations actually break:

- **Aliasing.** Every binary and unary routine in fp, fp2, fp12, and the curve
  and scalar-multiplication layers, with the output aliasing each input in turn.
  This section exists because two of the original audit's confirmed defects were
  exactly aliasing bugs, and no vector file would catch one: vectors never alias.
- **SHA-256 padding**, at sixteen lengths around the block boundaries, each also
  fed one byte at a time to exercise the buffered update.
- **expand_message_xmd limits**: zero bytes, one byte, exactly one block, the
  255-block maximum, one past it, and a tag of 255 against 256 bytes — the
  latter two must differ, since the oversize tag is hashed down rather than
  truncated.
- **Degenerate group law**: `P + (-P)`, `P + O`, `O + O`, `2O`, and `P + P`
  against `2P`.
- **Scalars at the ends**: 0, 1, `r-1`, `r`, and GLV against the plain ladder at
  each.
- **The trust boundary**: the identity refused by the pairing, every one of the
  eight flag combinations on a compressed encoding (exactly two may decode), an
  all-zero buffer, a point off the curve.

Everything except the `fp2_sqrt` aliasing case passed on the first run.

### The domain separation tag is no longer optional

`elips_hash_to_g1`, `hash_to_g2`, `encode_to_g1` and `encode_to_g2` now return
`int` and refuse `dst_len == 0`, leaving the identity behind.

RFC 9380 requires a non-empty tag, and hashing without one silently removes the
domain separation the argument exists to provide. It is also an easy mistake —
an uninitialised length, `sizeof` on a pointer — and the library already returns
a status everywhere else that a caller can get it wrong. Consistency, and one
line at each call site.

### examples/

Three programs, built for the selected curve, run by CTest and by CI on all
three curves so they cannot drift away from the API:

| Example | Shows |
|---|---|
| `01_pairing.c` | bilinearity, non-degeneracy, mu_r membership, input validation, and where the time goes |
| `02_hash_and_serialize.c` | hash-to-curve, both encodings, domain separation demonstrated rather than asserted, and each thing the deserializer refuses with the message it gives |
| `03_bls_signature.c` | keygen, sign, verify, tamper detection, wrong key, the identity offered as a signature, and aggregation |

The BLS one is the flagship and is explicit about what it is not: it uses its
own tag rather than the IETF ciphersuite, and it has **no proof of possession**,
so its aggregation section names the rogue-key attack rather than leaving a
reader to copy an insecure pattern. It also makes the §10.6 multi-pairing case
concrete — verification computes two pairings and so runs the final
exponentiation twice.

The pairing example prints the cost of the subgroup checks directly: 1256 µs of
Miller loop plus 1639 µs of final exponentiation against 4410 µs for
`elips_pairing`, so the validation is more than a third of the call. That is the
§10.3 argument, visible without reading the plan.

### What the gate caught immediately

Turning `-Werror` on across the matrix failed the macOS jobs on the first run,
which is the gate doing its job on its first day:

```
src/bls12_scm.c:107: error: 'sprintf' is deprecated ... [-Werror,-Wdeprecated-declarations]
```

Apple's SDK marks `sprintf` deprecated; glibc does not, so gcc and Linux clang
had been silent about five call sites in the legacy split-scalar and G3
exponentiation routines. Replaced with `snprintf`, and the two files now include
`<stdio.h>` themselves rather than relying on a transitive include.

Not an overflow -- Phase 0 checked that `char str[5]` holds its four characters
and terminator, and that finding still stands. It was a deprecation, invisible
on one platform and fatal on the other, which is exactly the class of problem a
two-architecture matrix exists to find and could not while nothing set
`-Werror`.

### Where the suite stands

43 CTest targets on Release (40 tests plus the 3 examples), 27 under each
sanitizer, zero warnings with `-Werror` on gcc and clang, all three curves.

---

## Phase 5d — the legacy layer is gone (issue #17)

**Status: complete. Exit criteria met. Defects A3 and A4 are gone with it, and
issue #16 no longer has anything to describe.**

### What went

| Deleted | Size |
|---|---|
| `src/*.c` — the runtime-curve mpz layer | 44 files, 7241 lines |
| `include/ELiPS_bn_bls/*.h` | 47 headers |
| `test/kat_runner.c` (legacy vector runner) | superseded |
| `test/property_tests.c` | its properties are covered by `standalone_*` and `edge_*` |
| `test/finalexp_agreement.c` | pinned issue #16 as an expected failure; the code it described is gone |
| `test/main.c`, `dump_pairing.c`, `diag_pairing.c` | orphans, in no build |

**12,271 lines deleted against 659 added.** Tracked files, excluding the
generated Doxygen site, drop from 169 to 76. The whole library is now 5408
lines, of which 2495 are generated parameter tables.

### The one piece of real work: replacing the reference

`test/pairing_test.c` used to say "the new pairing agrees with the old one".
That is the weaker of the two available statements — two implementations can
share a misreading of the twist conventions and agree on a wrong answer — and it
was the last thing keeping the old layer alive.

`tools/reference/gen_pairing_vectors.py` now emits `test/kat/pairing_*.vec` from
`pairing_ref.py`, the independent optimal ate written from the defining
equations, with the final exponentiation computed as one exponentiation by
`(p^12-1)/r` rather than by any addition chain. The points are read back out of
the generated `fp_params.h`, so the vectors pin the pairing on exactly the
generators the library ships, and the generator script asserts they are on their
curves and of order `r` before using them.

Four pairings per curve — `e([a]G1, [b]G2)` for `(a,b)` in `(1,1) (7,1) (1,11)
(5,9)` — so the file witnesses bilinearity with concrete values rather than only
asserting it. Eight checks per record: both inputs in their subgroups, the exact
value against the reference, non-degeneracy, `mu_r` membership, the fast chain's
relationship to the exact one, and `elips_pairing` agreeing with the pieces it is
built from.

**Only the value after the full final exponentiation is recorded.** A raw Miller
value depends on how the line functions are normalised, and those constants
cancel in the final exponentiation, so pinning one would pin a convention rather
than the pairing — and that convention already changed once, when the curve
layer moved to RCB.

**BLS12-381 has a pairing reference for the first time.** The legacy layer never
supported it, so while it was the reference side that curve's pairing value was
checked against nothing.

`test/ec_test.c` also stopped borrowing a G2 point from the old layer; it uses
`ep2_generator` and its multiples, which the parameter generator asserts is on
the twist and of order exactly `r` — stronger provenance than the borrowed point
had.

### A3 and A4, confirmed rather than assumed

Both were structural defects of the legacy final exponentiation: A3 wrote
through its input pointer, A4 mutated global curve parameters mid-computation.
Phase 2 left them alone deliberately.

They are gone, and the check is mechanical rather than a reading: every object
file in the library was inspected for non-const file-scope symbols, and there
are none. The library has no mutable global state left for that class of defect
to live in. The "curve context struct" the plan called for was never needed —
the constants are compile-time.

### The suite after

38 CTest targets, down from 43 but covering more:

- `kat.*` now runs the Phase 0 field and curve vectors on all three curves
  against the only remaining implementation, with its corruption control.
- `pairing.*` is new on all three curves, with its own corruption control.
- `property.*` and `finalexp.agreement.*` are gone with the code they tested.

Every suite is now driven by the Python oracle. Nothing in `test/` compares one
C implementation against another.

28 under each sanitizer, zero warnings under `-Werror` on gcc and clang.

### Packaging

One library, `ELiPS::arith`, built for the curve `-DELIPS_CURVE` selects.
`ELiPS::elips` is kept as an alias of it, so a build file written against the
old name still configures — verified by building a consumer against each name.
`include/ELiPS_bn_bls/` is no longer installed because it no longer exists.

### The Python tooling the deletion stranded

Three reference scripts were written to study the legacy final exponentiation.
Removing that code changed what each is worth, and they did not all land in the
same place.

`check_finalexp.py` and `find_exponent.py` both read a JSON dump produced by
`test/dump_pairing.c`. That producer is gone, so neither script can be run at
all: there is no longer any way to give them an input. They are deleted. What
they established is not lost, because it is what the current pairing test
asserts on every run.

`trace_finalexp.py` was the third, and it was the one worth keeping. It traced
`src/bn_final_exp.c` and `src/bls12_finalexp.c`, so after the deletion it still
passed while describing nothing. Rather than drop it or leave a CI step that
guards nothing, it now traces the two `pairing_final_exp_fast` chains in
`src/pairing/miller.c`, and it asserts rather than reports:

| curve | chain claims | trace says |
|---|---|---|
| BN-462 | `lambda` | exact |
| BLS12-461 | `3*lambda` | exact |
| BLS12-381 | `3*lambda` | exact |

The check is cheap and independent in a way the vectors are not. A vector says
the pairing produced the right number for four specific inputs. This says the
chain computes the exponent its comment claims, for every input, by tracking the
exponent through the chain instead of the field element: `fp12_mul` adds,
`fp12_sqr_cyc` doubles, `fp12_exp_param` multiplies by the signed x,
`fp12_frobenius` multiplies by a power of p, `fp12_conj` negates. The last of
those holds only because the easy part already put the element in the
cyclotomic subgroup, which is the same fact both C chains rely on.

Confirmed it fails when it should: changing one `frobenius(b,1)` to
`frobenius(b,2)` in the traced BLS12 chain makes both BLS12 curves report a
mismatch and the script exit 1. The CI step is renamed accordingly — it used to
say "documents issue #16", and issue #16 is moot now that the code it described
is deleted.

### Note for whoever repoints GitHub Pages

`docs/` is still tracked and is now 1280 of the 76+1280 tracked files —
94% of the repository is a generated Doxygen site describing code that in large
part no longer exists. Untracking it is a one-line change once Pages serves from
`gh-pages`; the workflow for that already exists.

## A dudect finding: the fixed class must not be degenerate

`dudect.ep2_mul_glv` failed on macOS Release with max|t| = 57.44, confirmed at
70.82 on an independent second sample. It passed on every Linux job and on
macOS Asan and Ubsan, and the commit that exposed it changed no C code at all.

**It was the harness, and the reason is worth writing down.** Class 0 was
`k = 0`. With a zero scalar every digit of the GLV decomposition is zero, the
ladder selects the identity at every step, and the accumulator sits at infinity
from the first doubling to the last. Class 0 therefore ran every field
multiplication on all-zero operands while class 1 ran on random ones.

`ep2_mul_glv` was read line by line before the harness was touched, and it is
structurally constant time: every loop bound is a public bit length, the
restoring division's conditional subtraction is applied by mask, and the
16-entry table is scanned in full with `fp2_cselect` on each entry. The
instruction counts of the two classes are identical. Checked the compiler too,
since a mask-select turning into a branch at -O2 is the classic way this goes
wrong: clang 18 at -O2 emits `cmove` for the division step, no branch.

So the two classes differed only in operand values, and what the t-test found
was operand-value dependence in the hardware. That reading is consistent with
where it appeared: Apple silicon at -O2 only, quiet at the lower optimisation
levels the sanitizer builds use, and 0.8 to 2.3 across ten runs on x86-64.

That comparison is not the threat model. The question a leakage test has to ask
is whether timing separates one realistic secret from another, so both classes
must be realistic. Class 0 is now a pinned pseudorandom scalar, `splitmix64`
from a fixed seed, masked below `2^(ORDER_BITS-1)` so it is in range without a
reduction. Verified nonzero, less than r, and roughly half its bits set on all
three curves.

This is the second time this file measured something other than the routine.
The first was `fp_add` reading |t| = 98 because random inputs were generated
inside the timed loop for one class only. Both failures share a shape: the
routine was fine and the experiment was not. The negative control is what makes
the difference visible either way, and it still fires at 909, so the harness has
not been blunted.

Ten consecutive local runs after the change: 0.79 to 3.03. Whether macOS Release
agrees is what the next CI run answers; if it still reports a leak the finding is
real and the hardware reading above is wrong.

---

## Phase 6 — fast subgroup tests (issue for plan §10.3)

The subgroup checks were `[r]P`, a full-length ladder, and they were the
largest single cost inside `elips_pairing`. They are now endomorphism tests.

**What made this safe to do.** A subgroup test that is too strict fails the
first time a valid point is checked. One that is too permissive fails silently:
it accepts attacker-chosen points of small order while every other test in the
suite still passes. That asymmetry is why this was deferred, and it decided how
it got built.

The full texts of 2021/1130, 2022/352 and 2022/348 were not reachable from this
session, so the deciding argument is machine-checked rather than recalled:

Each test asks whether `E(P) == [m]P` for an endomorphism satisfying a known
quadratic. A point of prime order `l` that satisfies the test forces `m` to be a
root of that quadratic mod `l` — `m²+m+1` for `phi` on G1, `m²−t·m+p` for `psi`
on G2. So the test is exact when no prime dividing the cofactor divides that
integer, and finding a shared prime is a gcd, not a factorisation.

Three results fall out, and two of them are better than expected:

- **BLS12 G1 is exact for every seed.** With `m = −x²`, `m²+m+1 = x⁴−x²+1 = r`
  exactly. `r` is prime and misses the cofactor, so nothing outside G1 can pass.
  The "no exceptional BLS12 seeds" result is a one-line identity here.
- **BN G1 needs no test at all.** `#E(Fp) = r`, cofactor 1, so on-curve is the
  whole condition. 939.6 µs to 0.9 µs.
- **G2 has no such identity** and is the case the literature had to correct
  (2022/352 fixing 2021/1130). Checked per curve by gcd; all three pass.

**Where the checks live, in order of how much they would catch.**
`tools/reference/subgroup_ref.py` derives each test, proves exactness, and then
confirms it numerically on points that really are outside the subgroup.
`gen_params.py` re-asserts the exactness conditions on every regeneration, so a
curve added later cannot ship a permissive test. `test/kat/subgroup_*.vec`
carries points to accept and points to reject; `test/subgroup_test.c` checks
both directions and refuses to run on a file with nothing to reject. The
negative control relabels a must-reject point as must-accept and requires the
runner to fail.

**Measured, x86-64 gcc -O2:**

| curve | G1 before | G1 after | G2 before | G2 after | pairing |
|---|---|---|---|---|---|
| BLS12-381 | 331.5 µs | 187.7 µs | 1131.4 µs | 333.1 µs | 4069 to 3127 µs, −23% |
| BLS12-461 | 613.7 µs | 350.4 µs | 2078.2 µs | 580.9 µs | 6900 to 5139 µs, −26% |
| BN-462 | 939.6 µs | 0.9 µs | 2936.3 µs | 1502.0 µs | 8701 to 6328 µs, −27% |

Short of the plan's "roughly a third", which assumed the checks disappeared
rather than becoming shorter ladders. BN's G2 multiplier is `6x²` at half the
width of `r`, so BN gains least on G2 and most on G1.

Also new: `ep_eq`/`ep2_eq`, projective equality by cross-multiplication with no
inversion, correct for infinity in both slots; and `ep_phi`, the GLV
endomorphism on G1.

42 CTest on Release, 32 under each sanitizer, no warnings under `-Werror`.

---

## Phase 6 — multi-pairing (plan §10.6, first half)

`elips_pairing_multi(out, P, Q, n)` computes the product of n pairings with one
shared Miller loop and one final exponentiation.

Two savings, and they are not equal. The per-iteration `fp12_sqr` now happens
once instead of n times, which helps more as n grows. The final exponentiation
runs once instead of n times, which is the larger saving at small n and is
exactly where a verification equation spends its wasted work: the final
exponentiation is 49% of a pairing on this profile.

**The correctness statement is exact.** Raising to a fixed exponent is a
homomorphism, so

    final_exp(f_1 * ... * f_n) == final_exp(f_1) * ... * final_exp(f_n)

bit for bit, including the stray cube the BLS12 fast chain contributes, because
both sides carry it. That makes the multi-pairing checkable against the
single-pairing path with exact equality rather than a tolerance.

`test/pairing_test.c` checks two independent things, since either can pass while
the other fails:

- against the **oracle**, the product of the values recorded in
  `test/kat/pairing_*.vec`, a number this library never computed;
- against **itself**, the product of individual `elips_pairing` results.

It also runs a ten-term product to cross the internal chunk boundary, checks
that the empty product is one and is not an error, and checks that an identity
input is refused with the result left at one rather than a partial product.

**Confirmed the tests can fail.** Squaring once per term instead of once per
iteration, and dropping the last term of each chunk, were each injected and each
made three checks fail. A test for a faster path that only ever agrees with the
faster path proves nothing.

**Chunking.** Per-term state is one `ep2_t`, so it sits in a fixed-size buffer
(`ELIPS_MULTI_CHUNK`, 8) and longer inputs run in chunks whose Miller values are
multiplied together. Any n works and the single final exponentiation survives
for any n; only the sharing of squarings stops at a boundary. At 8 a chunk
already captures 7/8 of the squaring saving, and a two-term verification never
reaches a boundary.

**Measured, x86-64 gcc -O2, against the same pairings computed separately:**

| terms | BLS12-381 | BN-462 |
|---|---|---|
| 2 | 1.45x | 1.32x |
| 4 | 1.78x | 1.54x |
| 8 | 2.08x | 1.74x |
| 16 | 2.24x | 1.79x |

The n=16 row crosses a chunk boundary, which is why it gains little over n=8.
That is the design behaving as described, not a surprise.

`examples/03_bls_signature.c` now verifies as `e(-sigma, G2) * e(H(m), pk) == 1`
in one call. Every forgery case it demonstrates still rejects, which is the
check that matters when rearranging a verification equation.

Also corrected: the comment on `ep_in_subgroup` in the public header still said
"is [r]P the identity?", which stopped being true when the endomorphism tests
landed.

Still open from §10.6: fixed-argument precomputation (ePrint 2010/342), worth a
further 25 to 37 percent of the Miller loop when one argument is fixed. It is
separable from this and untouched.

42 CTest on Release, 32 under each sanitizer, no warnings under `-Werror`.

---

---

# HAND-OFF — next session starts at Phase 6

**Phases 0 to 5b are complete.** Phase 6 (assembly) is next, and §10.10 of the
plan changed its target: read that before starting.

## Where things stand

Branch `claude/pairing-crypto-modernize-2faff2`.
38 CTest targets green on Release (including the three examples), 28 on Asan and
Ubsan (dudect is not registered under sanitizers), zero warnings under
`-Werror` on gcc and clang. CI sets `-DELIPS_WERROR=ON`, so that is a gate
rather than a claim.

**One layer.** The legacy mpz layer was retired by issue #17; every test is now
driven by the Python oracle rather than by a second C implementation.

| Curve | pairing value | serialized G1/G2 | hash-to-curve | GLV |
|---|---|---|---|---|
| BLS12-381 | `e^3` | 48/96 | RFC 9380 `SSWU_RO_`, **byte-exact** | G2 |
| BLS12-461 | `e^3` | 58/116 | RFC 9380 `SVDW_RO_`, no registered suite | G2 |
| BN-462 | `e` exact | 59/118 | RFC 9380 `SVDW_RO_`, no registered suite | none |

For BLS12-381 the generators, the point encodings and the hash-to-curve outputs
are all pinned against published values. The other two curves have no standard,
so their vectors pin self-consistency.

`cmake -B build -DELIPS_CURVE=BLS12_381` selects the installed curve;
`find_package(ELiPS)` gives `ELiPS::arith` and sets `ELiPS_CURVE`.

## Start here

1. **Read `MODERNIZATION_PLAN.md` §10.** Every open item now has a decision, a
   citation and a status. §10.10 revises the Phase 6 assembly target and §10.3
   is the largest available win.
2. **Take a profile on an x86-64 machine.** Every measurement through Phase 4
   was on Apple Silicon and Phase 5's on a slow VM. Phase 6 is conditional on a
   profile and that profile has not been taken on the target architecture.

## Ranked by value, from the survey

| Work | Worth | Risk |
|---|---|---|
| Fast subgroup tests (§10.3) | ~1/3 of the pairing | conditions must be derived and asserted, not recalled |
| AVX-512 IFMA on x86-64 (§10.10) | ~3x on that architecture | large, and needs the C fallback beside it |
| Multi-pairing + fixed-argument precomputation (§10.6) | ~40% of a signature verification | low; it is bookkeeping over the existing loop |
| Karabina compressed squaring (§10.5) | 10–15% of the final exponentiation | low |
| GLV on G1, and on BN G2 (§10.4) | ~1.4x on those ladders | low; the constant-time decomposition already exists |

## Open issues

- **#17 retire the legacy mpz layer** — done, Phase 5d. Ready to close.
- **#16 final exponentiation exponent** — the code it described is deleted.
  Ready to close.
- **#15 baseline measurements.** Superseded by the RELIC comparison.

## Things a fresh session should not re-derive

Everything in the previous hand-offs still holds. Added by Phase 5b:

- All three primes are `3 mod 4`; square roots are one exponentiation by
  `(p+1)/4`, and Fp2 roots follow Adj and Rodríguez-Henríquez.
- The serialized field width is `(FP_BITS + 3 + 7) / 8`, not `ceil(FP_BITS/8)`.
- BLS12-381's G2 generator is the published one and its encoding is a pinned
  known answer. **Do not regenerate it from a search.**
- `h_eff` for BLS12 G1 is `1-x` and is **not** a multiple of the cofactor. That
  is correct — `E(Fp)` is not cyclic there. Do not "fix" it.
- The Budroni–Pintore G2 cofactor chain is verified against `[h_eff]` on random
  points of the **twist**. Checking it on G2 points proves nothing.
- dudect must prepare every input before it times anything, and its fixed
  class must be a representative secret. `k = 0` is not one: it holds the
  accumulator at infinity and compares all-zero operands against random
  ones, which is a hardware question, not a control-flow one.
- A whole-`fp6` or whole-`fp2` access to an element of a nested array parameter
  trips GCC's `-Wstringop-overflow`. Pass the halves, or take an element macro.
- Every routine in `fp.h`, `fpx.h` and `ec.h` must tolerate its output aliasing
  its inputs, and `test/edge_test.c` checks it. Write the result into a local
  and store it once at the end; `fp2_sqrt` was written the other way and was
  wrong for `r == a` until the aliasing cases went looking.
- The four hash-to-curve entry points return `int` and refuse an empty domain
  separation tag. That is deliberate, not an oversight to simplify away.

## Known limitations, stated plainly

- Subgroup checks are full scalar multiplications: 2175 µs of `elips_pairing`'s
  6033 µs on BLS12-381. Correct, and slow. §10.3 is the fix.
- No multi-pairing, so a two-pairing product costs two final exponentiations.
- BN-462 has no GLV and no fast G2 cofactor chain.
- `hash_to_g2` rebuilds a window table for each of its two short ladders; a
  shared table would help.
- No signature layer. §10.8 explains why that line is where it is.
- `docs/` is still tracked: 1280 files against 76 of everything else, so 94% of
  the repository is a generated Doxygen site, much of it describing code that no
  longer exists. GitHub Pages serves from `master:/docs`; a workflow to publish
  from `gh-pages` exists, and once Pages is repointed `git rm -r --cached docs`
  drops the repository to 76 tracked files.
- The IACR survey behind §10 was done through search abstracts: the egress
  policy blocks `eprint.iacr.org`, so no full text was read. Citations are
  pointers, not sources of copied formulas.
