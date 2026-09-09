# Pairing research archive

Research only. These variable-time test programs use deterministic, public
inputs. They are not production APIs or constant-time implementation claims.
Keep exploratory source, including unsuccessful approaches, in this directory.
Do not remove experiments when a newer candidate replaces them.

**Start here:** [the beginner research note](RESEARCH_NOTE.md) develops the
field basics, both kernel derivations, worked examples, code map, and
reproduction steps. Run its companion with
`python3 tools/research/pairing/tutorial_examples.py` from the repository root.

## Contents

- `normalized-line/`: original C 15M2/10M2/8M2 comparison, tangent identity
  checks, original report, and two original timing logs. The Python prototype
  remains in `tools/reference/normalized_miller_ref.py`.
- `norm-and-exponent/`: prescribed-norm experiment, original primitive-ratio
  benchmark, reconstructed LLL/Babai search, and historical out-of-box notes.
  **The old notes contain retracted claims; read `CORRECTIONS.md` before
  citing them.**
- `scaled-square/`: Gaussian and Fourier interpolation squarings, a published
  asymmetric-squaring control, independent Python checks, C benchmarks, and
  a timing-log analyzer. See `scaled-square/RESULTS.md` for this round's result.
- `inline_checks.py`: reconstructed symbolic identities, bounded polynomial
  searches, seed-weight checks, and line-support checks from inline commands.
- `IMPORT_MANIFEST.json`: SHA-256 hashes at import from the temporary research
  directories. Subsequent portability edits are intentional; these hashes
  record provenance, not the current files' contents.

All surviving research source files and notes were copied, not moved. The
original temporary files remain untouched. Generated binaries/build trees and
downloaded third-party papers are not included. Earlier inline shell-only
experiments are not covered by the import manifest.

## Reproduce

From the repository root, configure this opt-in standalone project. Choose a
new build directory; nothing needs to be installed into the system:

```sh
cmake -S tools/research/pairing -B /tmp/elips-research-build \
  -DCMAKE_BUILD_TYPE=Release -DELIPS_BUILD_TESTS=OFF -DELIPS_WERROR=ON
cmake --build /tmp/elips-research-build -j 4
ctest --test-dir /tmp/elips-research-build --output-on-failure

PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/scaled-square/gaussian_ref.py
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/norm-and-exponent/norm_state.py
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/norm-and-exponent/lattice_search.py
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/scaled-square/analyze_bench.py

/tmp/elips-research-build/scaled_square_experiment --bench
/tmp/elips-research-build/normalized_line_experiment --bench
/tmp/elips-research-build/primitive_ratios
```

These Python tests require only the standard library. The legacy
`norm_state.py` uses `assert`; do not run it with `python -O`. C requires the
same compiler, CMake, and GMP development files as the main library.

Optional symbolic checks need SymPy:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/inline_checks.py
```

The lattice script defaults to quick unit tests. To repeat the slow heuristic
search, pass `--search --max-x 5 --max-p 3` (positive control), or
`--search --max-x 4 --max-p 3`. It uses exact Fraction-based LLL and can take
minutes. Do not interpret its returned candidate as an optimality proof.

For ASan+UBSan, configure a separate directory with `-DCMAKE_BUILD_TYPE=Asan`,
then build and run CTest there. These programs are deliberately not added to
the normal library build or CI.

Scaled kernels preserve the result **after final exponentiation**, not the
raw Miller value. In particular, they must not silently replace the existing
bit-identical `pairing_miller_prec` interface.
