# ELiPS — Efficient Library for Pairing-based Systems

Pairing-based cryptography over BN and BLS12 curves, built on GMP.

> **Status: under active modernization.** See [`MODERNIZATION_PLAN.md`](MODERNIZATION_PLAN.md)
> for the roadmap and [`PROGRESS.md`](PROGRESS.md) for what has landed.
>
> **Known defect:** the optimal final exponentiation raises to the wrong
> exponent — `e^3` on BLS12 and `e^(~12·X³)` on BN. The result is still a valid
> bilinear pairing, but it does not match `finalexp_plain` in this same library
> and will not interoperate with other implementations.
> See [issue #16](https://github.com/eNipu/elips_bn_bls/issues/16).

## Requirements

- CMake 3.16 or newer
- A C11 compiler
- [GMP](https://gmplib.org/) (`libgmp-dev` on Debian/Ubuntu, `brew install gmp` on macOS)

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Test

```bash
ctest --test-dir build --output-on-failure
```

The suite runs known-answer vectors for BN-462, BLS12-461 and BLS12-381,
property tests including subgroup order checks, and a negative control that
confirms the vector runner actually detects corrupted input.

## Install

```bash
cmake --install build --prefix /usr/local
```

Then from a downstream project:

```cmake
find_package(ELiPS REQUIRED)
target_link_libraries(your_target PRIVATE ELiPS::elips)
```

## Sanitizer builds

```bash
cmake -B build-asan -DCMAKE_BUILD_TYPE=Asan     # address + undefined
cmake -B build-ubsan -DCMAKE_BUILD_TYPE=Ubsan   # undefined only
```

MemorySanitizer is not wired up: it needs every dependency instrumented, and
GMP is not, so it would report noise rather than findings.

## Test oracle

`tools/reference/` holds an independent Python implementation of the field
tower, curve arithmetic and the optimal ate pairing, written from the defining
equations rather than from `src/`. It generates the vectors in `test/kat/` and
is what found the final exponentiation defect.

```bash
python3 tools/reference/selftest.py        # field axioms and tower relations
python3 tools/reference/gen_vectors.py test/kat
python3 tools/reference/trace_finalexp.py  # exact exponent of each chain
```

A useful structural fact: the tower collapses to a single polynomial. Since
`v = w²` and `1+u = v³ = w⁶`, we have `u = w⁶−1` and therefore

```
Fp12 = Fp[w] / (w¹² − 2w⁶ + 2)
```

which is irreducible over all three primes. Sage can build exactly this field
rather than an abstract `GF(p¹²)`, so its pairing is comparable coefficient by
coefficient.

## Curves

| Curve | p | r | Status |
|---|---|---|---|
| BN-462 | 462 bits | 462 bits | supported |
| BLS12-461 | 461 bits | 308 bits | supported |
| BLS12-381 | 381 bits | 255 bits | field layer only; curve support in Phase 3 |

## Licence

LGPL. See [`COPYING.LESSER`](COPYING.LESSER).

## Contact

khandaker@s.okayama-u.ac.jp
