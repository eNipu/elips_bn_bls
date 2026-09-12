# Contributing

**Status: under active modernization.** See
[`MODERNIZATION_PLAN.md`](MODERNIZATION_PLAN.md) for the roadmap,
[`PROGRESS.md`](PROGRESS.md) for what has landed and why, and
[`PRD-developer-api.md`](PRD-developer-api.md) for the bindings and demo.

## Build and test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The test suite needs [GMP](https://gmplib.org/) (`libgmp-dev` on
Debian/Ubuntu, `brew install gmp` on macOS). It is kept deliberately, as an
independent oracle: the tests check the Montgomery layer, the wide reduction
and both inversions against it, and agreement with a separate implementation
is evidence in a way that agreement with ourselves is not. The library itself
does not link it -- configure with `-DELIPS_BUILD_TESTS=OFF` to build without
GMP present.

The test suite always builds all three curves, because the vectors have to
cover all three.

## Sanitizer builds

```bash
cmake -B build-asan -DCMAKE_BUILD_TYPE=Asan     # address + undefined
cmake -B build-ubsan -DCMAKE_BUILD_TYPE=Ubsan   # undefined only
```

MemorySanitizer is not wired up: it needs every dependency instrumented, and
the GMP the tests link is not, so it would report noise rather than findings.

## Constant time

The parts that handle a secret key are constant time, and it is tested rather
than asserted.

```bash
ctest --test-dir build -R dudect        # statistical, Welch t-test on timings
python3 tools/verify/ct_branch_scan.py build   # structural, reads the object code
```

`dudect` carries two negative controls that must leak, so a harness that has
stopped detecting anything fails instead of passing quietly. It also measures
its own noise ceiling per target per run, by reshuffling the class labels on
the samples it just took, so the verdict does not depend on a constant
measured on somebody else's machine.

`ct_branch_scan.py` reads the disassembly for masked selects the compiler
turned into branches. This has happened three times here (`wide.c`, then
issue #30, then issue #32) and a statistical test found each one only on the
toolchain that did it. Run both. They fail differently on purpose.

Every mask built as `0 - condition` should go through `ct_mask` in
`src/arith/ct.h`. The cost is one register move and the failure mode is
silent.

## Test oracle

`tools/reference/` holds an independent Python implementation of the field
tower, curve arithmetic and the optimal ate pairing, written from the defining
equations rather than from `src/`. It generates the vectors in `test/kat/` and
is what found the final exponentiation defect.

```bash
python3 tools/reference/selftest.py        # field axioms and tower relations
python3 tools/reference/h2c_ref.py         # RFC 9380 maps, isogenies, SvdW
python3 tools/reference/gen_vectors.py test/kat          # field and curve
python3 tools/reference/gen_pairing_vectors.py test/kat  # the pairing itself
python3 tools/reference/gen_h2c_vectors.py test/kat      # hash to curve
python3 tools/reference/subgroup_ref.py    # subgroup tests, derived and proved
python3 tools/reference/divstep_ref.py     # inversion model and its iteration count
python3 tools/reference/karabina_ref.py    # compressed squaring: derive, then decide
python3 tools/reference/gen_subgroup_vectors.py test/kat  # accept/reject points
python3 tools/reference/trace_finalexp.py  # exponent of each final-exp chain
python3 tools/reference/normalized_miller_ref.py  # experimental 8M2 prepared-line kernel
```

`normalized_miller_ref.py` is a Python-only research prototype, not a change to
the C pairing. It tests a normalized, subfield-scaled line multiplication that
uses eight Fp2 multiplications instead of fifteen. Raw Miller values change;
final-exponentiated values must match the pairing vectors. It also checks
multi-pairing, normalization costs, and invalid inputs. It makes no claim of
novelty or constant-time execution.

`trace_finalexp.py` also produces no vectors. It walks the
two `pairing_final_exp_fast` chains in `src/pairing/miller.c` symbolically,
tracking the exponent rather than the field element, and asserts that each
computes what its comment claims: `lambda` for BN, `3*lambda` for BLS12. A
numeric vector says the answer is right today; this says the chain is the
right chain, and fails if an edit changes its exponent.

Every suite in `test/` is driven by these vectors, so the library is checked
against a Python implementation written from the defining equations rather than
against a second C implementation of the same ideas. That distinction is why the
pairing reference survived retiring the old layer: the oracle was never the old
code.

A useful structural fact: the tower collapses to a single polynomial. Since
`v = w²` and `1+u = v³ = w⁶`, we have `u = w⁶−1` and therefore

```
Fp12 = Fp[w] / (w¹² − 2w⁶ + 2)
```

which is irreducible over all three primes. Sage can build exactly this field
rather than an abstract `GF(p¹²)`, so its pairing is comparable coefficient by
coefficient.
