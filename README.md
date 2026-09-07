# ELiPS — Efficient Library for Pairing-based Systems

Pairing-based cryptography over BN and BLS12 curves, built on GMP.

> **Status: under active modernization.** See [`MODERNIZATION_PLAN.md`](MODERNIZATION_PLAN.md)
> for the roadmap and [`PROGRESS.md`](PROGRESS.md) for what has landed.
>
> The API is `include/elips/*.h`: fixed-width Montgomery arithmetic,
> constant-time scalar multiplication, standard serialization and RFC 9380
> hash-to-curve. The original runtime-curve `mpz_t` layer under
> `include/ELiPS_bn_bls/` was retired by
> [issue #17](https://github.com/eNipu/elips_bn_bls/issues/17), taking with it
> the wrong-exponent final exponentiation of
> [issue #16](https://github.com/eNipu/elips_bn_bls/issues/16) and the two
> global-state defects that depended on it.

## Requirements

- CMake 3.16 or newer
- A C11 compiler
- [GMP](https://gmplib.org/) (`libgmp-dev` on Debian/Ubuntu, `brew install gmp` on macOS)

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Choosing a curve

The modern arithmetic is compiled for one curve: the field width has to be a
compile-time constant for the loops to unroll. Pick it with `ELIPS_CURVE`:

```bash
cmake -B build -DELIPS_CURVE=BLS12_381    # default, and the only standardised one
cmake -B build -DELIPS_CURVE=BLS12_461
cmake -B build -DELIPS_CURVE=BN_462
```

That choice decides which library is installed and what
`find_package(ELiPS)` hands back as `ELiPS::arith`; `ELiPS_CURVE` is set in the
package config so a consumer can read it back. The test suite always builds all
three regardless, because the vectors have to cover all three.

**Use BLS12-381 unless you have a reason not to.** It is the curve with a
specification, so it is the one whose generators, encodings and hash-to-curve
outputs match other implementations byte for byte.

## Test

```bash
ctest --test-dir build --output-on-failure
```

The suite runs known-answer vectors for BN-462, BLS12-461 and BLS12-381
(field, curve, pairing and RFC 9380 hash-to-curve), property tests including
subgroup order checks, edge cases and aliasing contracts, timing-leakage tests
with a negative control, the three examples, and negative controls that confirm
the vector runners actually detect corrupted input.

`-DELIPS_WERROR=ON` turns warnings into errors. CI sets it, so a warning cannot
reach the branch; it is off by default so a contributor is not blocked by
whatever their compiler happens to emit.

## Install

```bash
cmake --install build --prefix /usr/local
```

Then from a downstream project:

```cmake
find_package(ELiPS REQUIRED)
target_link_libraries(your_target PRIVATE ELiPS::arith)   # the modern layer
```

```c
#include "elips/hash_to_curve.h"
#include "elips/serialize.h"
#include "elips/pairing.h"

static const char DST[] = "MY-PROTOCOL-V01-CS01-" ELIPS_H2C_SUITE_G1;

ep_t P;  ep2_t Q;  fp12_t z;
elips_hash_to_g1(&P, msg, msg_len, (const uint8_t *)DST, sizeof DST - 1);
ep2_generator(&Q);
if (!elips_pairing(z, &P, &Q)) { /* a point failed validation */ }

uint8_t enc[EP_SER_COMPRESSED_BYTES];
ep_write_compressed(enc, &P);
```

A verification equation should be one call, not several pairings compared
afterwards. `elips_pairing_multi` shares one Miller loop and one final
exponentiation across every term, which is 1.45x for the two-term case on
BLS12-381 and 2.24x at sixteen terms:

```c
/* e(-sigma, G2) * e(H(m), pk) == 1 */
ep_t  Ps[2];  ep2_t Qs[2];
ep_neg(&Ps[0], &sig);   ep2_generator(&Qs[0]);
ep_copy(&Ps[1], &h);    ep2_copy(&Qs[1], &pk);

fp12_t prod, one;
fp12_set_one(one);
int ok = elips_pairing_multi(prod, Ps, Qs, 2) && fp12_eq(prod, one);
```

When the G2 arguments are reused, which in a verification they are, precompute
the line functions once and replay them. That removes the G2 point arithmetic
from the loop entirely:

```c
ep2_prec_t pc[2];                 /* built once, per fixed Q */
if (!ep2_precompute(&pc[0], &G2)) { /* Q is not in G2 */ }
if (!ep2_precompute(&pc[1], &pk)) { /* ... */ }

/* then, for each verification */
int ok = elips_pairing_multi_prec(prod, Ps, pc, 2) && fp12_eq(prod, one);
```

`ep2_precompute` checks the subgroup itself, because a table has no point left
to check afterwards. A table costs about two Miller loops to build and 20 KB
(BLS12-381) to hold, so it pays from the second use onward.

| | BLS12-381 | BN-462 |
|---|---|---|
| Miller loop, precomputed vs direct | 1.58x | 1.48x |
| one pairing | 1.33x | 1.55x |
| 2-term product vs 2 separate pairings | 2.02x | 2.52x |
| 8-term product vs 8 separate pairings | 3.85x | 4.58x |

`ELiPS::elips` is kept as an alias of the same library, so a build file written
against the old name still configures.

## Examples

Three runnable programs under `examples/`, built for the selected curve and run
by CTest so they cannot rot:

```bash
./build/examples/elips_example_pairing   # bilinearity, validation, timing
./build/examples/elips_example_hash      # hash-to-curve, encodings, what the
                                         # deserializer refuses and why
./build/examples/elips_example_bls       # a BLS signature end to end, with
                                         # aggregation and its rogue-key caveat
```

The BLS one is a demonstration of this API, not a signature implementation to
deploy: it uses its own domain separation tag rather than the IETF ciphersuite,
and it has no proof of possession, which aggregation needs. The file says so at
the top and explains what a real one would add.

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
python3 tools/reference/h2c_ref.py         # RFC 9380 maps, isogenies, SvdW
python3 tools/reference/gen_vectors.py test/kat          # field and curve
python3 tools/reference/gen_pairing_vectors.py test/kat  # the pairing itself
python3 tools/reference/gen_h2c_vectors.py test/kat      # hash to curve
python3 tools/reference/subgroup_ref.py    # subgroup tests, derived and proved
python3 tools/reference/divstep_ref.py     # inversion model and its iteration count
python3 tools/reference/karabina_ref.py    # compressed squaring: derive, then decide
python3 tools/reference/gen_subgroup_vectors.py test/kat  # accept/reject points
python3 tools/reference/trace_finalexp.py  # exponent of each final-exp chain
```

`trace_finalexp.py` is the odd one out: it produces no vectors. It walks the
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

## Curves

| Curve | p | r | Generators | Serialized G1 / G2 | hash-to-curve |
|---|---|---|---|---|---|
| BLS12-381 | 381 bits | 255 bits | from the specification | 48 / 96 bytes | RFC 9380 `SSWU_RO_`, byte-exact |
| BLS12-461 | 461 bits | 308 bits | generated | 58 / 116 bytes | RFC 9380 `SVDW_RO_`, no registered suite |
| BN-462 | 462 bits | 462 bits | generated | 59 / 118 bytes | RFC 9380 `SVDW_RO_`, no registered suite |

Sizes are the compressed encodings; uncompressed is twice each. BN-462 needs one
byte more than its field width because the three flag bits do not fit otherwise.

Only BLS12-381 has a specification to conform to. For it, the generators, the
point encodings and the hash-to-curve outputs are pinned against the published
values and match any conforming implementation. The other two curves have no
standard, so their vectors pin self-consistency rather than interoperability.

**On BLS12 the pairing returns `e³`, not `e`.** That is a property of the
standard final-exponentiation chain, which RELIC and the original library also
use, and not a defect: since `gcd(3, r) = 1` it is still bilinear and
non-degenerate, so any protocol that only compares pairings is unaffected. Raw
values will not match an implementation that outputs `e`; use
`pairing_final_exp_plain` when the exact value is needed. On BN the value is
`e` exactly.

## Licence

LGPL. See [`COPYING.LESSER`](COPYING.LESSER).

## Contact

khandaker@s.okayama-u.ac.jp
