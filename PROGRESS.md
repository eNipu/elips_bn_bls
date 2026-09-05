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

- [ ] #8 — fold the ad-hoc build commands into `tools/build_reference.sh`; push tag `phase0-reference`
- [ ] #14 — property tests as a standing suite (`[r]P = O`, `e(P,Q)^r = 1`, additivity, four-variant cross-check). The checks exist and have been run; they are not yet a committed test target.
- [ ] #15 — replace the single-sample baseline with repetitions and spread

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
