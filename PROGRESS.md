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
