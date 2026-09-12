# ELiPS

Pairing-based cryptography over BLS12 and BN curves, in C, with no runtime
dependencies.

If you are here to **verify a BLS signature**, start below and you will have
one working in a couple of minutes. If you are here for the **pairing itself**,
skip to [Using the C library](#using-the-c-library).

**[Try it in your browser](https://enipu.github.io/elips_bn_bls/demo/)** — many
signatures collapsing into one, computed in the page, no backend.

---

## Signatures in two minutes

BLS signatures on BLS12-381: 32-byte keys, 48-byte public keys, 96-byte
signatures, and any number of signatures aggregate into one that is still 96
bytes.

> The `elips` packages are **not on PyPI or npm yet**. Until they are, install
> from a clone. Everything below is a command you can paste.

### Python

```bash
git clone https://github.com/eNipu/elips_bn_bls
pip install ./elips_bn_bls/bindings/python
```

A C compiler is the only requirement. There is no system library to install.

```python
import elips

sk  = elips.keygen()
pk  = elips.sk_to_pk(sk)
sig = elips.sign(sk, b"hello")

elips.verify(pk, b"hello", sig)          # returns None, or raises
```

**Verification raises; it does not return a boolean.** A verify that returns a
value invites `if verify(...): accept()`, which is wrong in the direction that
matters. Returning `None` makes that branch never fire, so the mistake shows up
on the first *valid* signature rather than on the first forged one.

```python
try:
    elips.verify(pk, msg, sig)
except elips.ElipsError:
    reject()
```

Aggregation, which is the reason to choose BLS at all:

```python
agg = elips.aggregate([elips.sign(sk, m) for sk, m in zip(sks, msgs)])
len(agg)                                 # 96, whatever the number of signers
elips.aggregate_verify(pks, msgs, agg)
```

Full reference: [`bindings/python/README.md`](bindings/python/README.md).

### JavaScript and the browser

```bash
cd elips_bn_bls/bindings/js
. /path/to/emsdk/emsdk_env.sh && ./build.sh
```

```js
import { load } from "./elips.mjs";
const elips = await load();

const sk  = elips.keygen();
const pk  = elips.skToPk(sk);
const msg = new TextEncoder().encode("hello");
const sig = elips.sign(sk, msg);

elips.verify(pk, msg, sig);              // returns undefined, or throws
```

One 68 KB file with the WebAssembly embedded, so there is no `.wasm` path to
configure and it works the same from a `<script type="module">` tag, a bundler
and node.

Full reference: [`bindings/js/README.md`](bindings/js/README.md).

## What this is, and what it is not

**It is** an implementation of `draft-irtf-cfrg-bls-signature` on BLS12-381,
proof-of-possession scheme, minimal-pubkey-size — the variant Ethereum uses.
Every operation is checked against
[py_ecc](https://github.com/ethereum/py_ecc), an independent implementation, on
every build, in both directions, because an implementation can be
self-consistently wrong and pass a one-way check. The field arithmetic,
hash-to-curve and pairing underneath are pinned against RFC 9380's published
vectors and against a Python oracle written from the defining equations.

The parts that handle a secret key are constant time, and that is tested rather
than asserted: `dudect` runs on every build with negative controls that must
leak, so a harness that has stopped detecting anything fails instead of passing
quietly.

**It is not audited.** Implementing the draft and agreeing with another
implementation is a different claim from having been reviewed by a
cryptographer. There are also two known open issues worth reading before
trusting it with anything: [#30](https://github.com/eNipu/elips_bn_bls/issues/30),
an unexplained timing signal in the GLV scalar multiplication under clang, and
the fact that no formal review has happened at all.

## How fast

Measured on an Intel Xeon, same machine. The pairing figure is the one
recorded in [`bench/baseline.json`](bench/baseline.json); the sign and verify
figures were taken in one run against each other. `bench/compare.py` alternates
two builds rather than comparing separate runs, because on this machine the
absolute numbers wander by about 20% between runs while the ratios between
operations do not.

| | native | WebAssembly | |
|---|---|---|---|
| sign | 1.77 ms | 7.41 ms | 4.2x |
| verify | 3.43 ms | 18.56 ms | 5.4x |
| pairing | 1.39 ms | | |

The browser is slower for a structural reason rather than a missing
optimisation: `wasm32` has no 64×64 → 128 bit multiply, so every field
multiplication goes through a software helper and a pairing is tens of
thousands of them. Native x86-64 and AArch64 do it in one instruction. The
number is published rather than omitted because omitting an unflattering
measurement is the one thing this repository has consistently refused to do.

---

# Using the C library

The layer the bindings are built on. Two headers matter:
`include/elips/bls.h` for signatures, and `include/elips/pairing.h` with its
neighbours for the pairing itself.

## Requirements

- CMake 3.16 or newer
- A C11 compiler

That is all. **The library has no runtime dependencies**: nothing under `src/`
includes a third-party header, and the installed library links nothing but
libc.

To build and run the **test suite** you also need
[GMP](https://gmplib.org/) (`libgmp-dev` on Debian/Ubuntu, `brew install gmp`
on macOS). GMP is kept deliberately, as an independent oracle: the tests check
the Montgomery layer, the wide reduction and both inversions against it, and
agreement with a separate implementation is evidence in a way that agreement
with ourselves is not. Configure with `-DELIPS_BUILD_TESTS=OFF` to build
without it.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Signatures, in C

```c
#include "elips/bls.h"

uint8_t sk[ELIPS_BLS_SK_BYTES], pk[ELIPS_BLS_PK_BYTES];
uint8_t sig[ELIPS_BLS_SIG_BYTES];

elips_bls_keygen_random(sk);
elips_bls_sk_to_pk(pk, sk);
elips_bls_sign(sig, sk, msg, msg_len);

if (elips_bls_verify(pk, msg, msg_len, sig) == ELIPS_BLS_OK) {
    /* good */
}
```

Nothing returns a boolean; every function returns `ELIPS_BLS_OK`, which is
zero, or a negative code. `if (elips_bls_verify(...))` is therefore true on
*failure*, which is the safe way round. Link against `ELiPS::bls`.

`elips_bls` is a separate CMake target that consumes only the public headers of
the core, so a project that wants the pairing and not the protocol can ignore
it entirely.

## Choosing a curve

The arithmetic is compiled for one curve: the field width has to be a
compile-time constant for the loops to unroll and for the assembly to target a
known width.

```bash
cmake -B build -DELIPS_CURVE=BLS12_381    # default, and the only standardised one
cmake -B build -DELIPS_CURVE=BLS12_461
cmake -B build -DELIPS_CURVE=BN_462
```

**Use BLS12-381 unless you have a reason not to.** It is the curve with a
specification, so it is the one whose generators, encodings and hash-to-curve
outputs match other implementations byte for byte. The signature layer is
BLS12-381 only, because the other two have no registered ciphersuite and an
"IETF BLS signature on BN-462" would interoperate with nothing.

The test suite always builds all three, because the vectors have to cover all
three.

## Install and consume

```bash
cmake --install build --prefix /usr/local
```

```cmake
find_package(ELiPS REQUIRED)
target_link_libraries(app PRIVATE ELiPS::arith)   # or ELiPS::bls
```

`ELiPS_CURVE` is set in the package config so a consumer can read back which
curve was installed. No dependency is imposed on the consumer.

## Examples

Three runnable programs under `examples/`, built for the selected curve and run
by CTest so they cannot rot:

```bash
./build/examples/elips_example_pairing   # bilinearity, validation, timing
./build/examples/elips_example_hash      # hash-to-curve, encodings, what the
                                         # deserializer refuses and why
./build/examples/elips_example_bls       # the pairing identity behind BLS,
                                         # built from the low-level API
```

The BLS example predates `elips/bls.h` and builds a signature scheme out of the
primitives by hand. It is useful for seeing what the pairing is doing and is
explicitly not a scheme to deploy: it uses its own domain separation tag rather
than the IETF ciphersuite and has no proof of possession. **Use
`elips/bls.h`**, which has both.

---

# Contributing

**Status: under active modernization.** See
[`MODERNIZATION_PLAN.md`](MODERNIZATION_PLAN.md) for the roadmap,
[`PROGRESS.md`](PROGRESS.md) for what has landed and why, and
[`PRD-developer-api.md`](PRD-developer-api.md) for the bindings and demo.

## Sanitizer builds

```bash
cmake -B build-asan -DCMAKE_BUILD_TYPE=Asan     # address + undefined
cmake -B build-ubsan -DCMAKE_BUILD_TYPE=Ubsan   # undefined only
```

MemorySanitizer is not wired up: it needs every dependency instrumented, and
the GMP the tests link is not, so it would report noise rather than findings.

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
