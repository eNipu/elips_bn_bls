# ELiPS

Pairing-based cryptography over BLS12 and BN curves, in C, with no runtime
dependencies.

If you are here to **verify a BLS signature**, start below. If you are here for
the **pairing itself**, skip to [the C library](#the-c-library).

**[Try it in your browser](https://enipu.github.io/elips_bn_bls/demo/)** — many
signatures collapsing into one, computed in the page, no backend.

---

## Quick start

BLS signatures on BLS12-381: 32-byte secret keys, 48-byte public keys, 96-byte
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
than asserted: a statistical timing test and a scan of the generated object
code both run on every build. See
[Constant time](CONTRIBUTING.md#constant-time) for what each one catches and
why both are needed.

**It is not audited.** Implementing the draft and agreeing with another
implementation is a different claim from having been reviewed by a
cryptographer. No formal review has happened. Do not use it to protect
anything that matters yet.

## Performance

Measured on an Intel Xeon, same machine. The pairing figure is the one recorded
in [`bench/baseline.json`](bench/baseline.json); sign and verify were taken in
one run against each other, because on this machine the absolute numbers wander
by about 20% between runs while the ratios do not.

| | native | WebAssembly | |
|---|---|---|---|
| sign | 1.77 ms | 7.41 ms | 4.2x |
| verify | 3.43 ms | 18.56 ms | 5.4x |
| pairing | 1.39 ms | | |

The browser is slower for a structural reason rather than a missing
optimisation: `wasm32` has no 64×64 → 128 bit multiply, so every field
multiplication goes through a software helper and a pairing is tens of
thousands of them. Native x86-64 and AArch64 do it in one instruction.

---

# The C library

The layer the bindings are built on. Two headers matter:
`include/elips/bls.h` for signatures, and `include/elips/pairing.h` with its
neighbours for the pairing itself.

## Build

Requirements are CMake 3.16 or newer and a C11 compiler. That is all: **the
library has no runtime dependencies**. Nothing under `src/` includes a
third-party header, and the installed library links nothing but libc.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Building the test suite additionally needs GMP, which is used as an
independent oracle rather than by the library. See
[CONTRIBUTING.md](CONTRIBUTING.md).

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

## Curves

The library is compiled for one curve, chosen with `-DELIPS_CURVE`.
**Use BLS12-381 unless you have a reason not to**: it is the only one with a
specification, so it is the only one that interoperates. BLS12-461 and BN-462
are also supported, and the signature layer is BLS12-381 only.

Sizes, generators, encodings and the BLS12 `e³` convention: [`CURVES.md`](CURVES.md).

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

## Contributing

Build and test instructions, the sanitizer builds, the constant-time checks and
the Python reference oracle: [CONTRIBUTING.md](CONTRIBUTING.md).

## Licence

LGPL. See [`COPYING.LESSER`](COPYING.LESSER).

## Contact

khandaker@s.okayama-u.ac.jp
