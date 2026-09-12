# @elips/bls

BLS signatures on BLS12-381, compiled to WebAssembly from the
[ELiPS](https://github.com/eNipu/elips_bn_bls) pairing library.

```js
import { load } from "@elips/bls";
const elips = await load();

const sk  = elips.keygen();
const pk  = elips.skToPk(sk);
const msg = new TextEncoder().encode("hello");
const sig = elips.sign(sk, msg);

elips.verify(pk, msg, sig);        // returns undefined, or throws
```

One 68 KB file. No network fetch, no `.wasm` path to configure, no bundler
plugin: the module is embedded, so the same import works from a
`<script type="module">` tag, from a bundler and from node.

## Verification throws. It does not return a boolean.

```js
try { elips.verify(pk, msg, sig); }
catch (err) { reject(); }
```

A verify that returns a value invites `if (verify(...)) accept()`, which is
wrong in the direction that matters. Returning `undefined` makes that branch
never taken, so the mistake surfaces on the first *valid* signature rather than
on the first forged one.

Catch `ElipsError` to reject. A `TypeError` or `RangeError` means the calling
code passed the wrong type or length — a bug in the caller, not a bad
signature — so it is a different error, thrown before any cryptography runs.

## Aggregation

`n` signatures combine into one, still 96 bytes, verified against all `n`
public keys at once.

```js
const agg = elips.aggregate(sigs);        // 96 bytes, whatever n was
elips.aggregateVerify(pks, msgs, agg);
```

Messages must be distinct. For one shared message use `fastAggregateVerify`,
which takes the proofs of possession and checks them:

```js
const pops = sks.map((sk) => elips.popProve(sk));
elips.fastAggregateVerify(pks, pops, msg, agg);
```

There is no way to call it without the proofs. Aggregating keys over a shared
message is forgeable otherwise: an attacker who registers
`pk_evil = [t]G1 - sum(pk_honest)` can produce an aggregate for a group whose
keys they never held.

## Speed: slower than native, and here is how much

Measured on an Intel Xeon, node 22 against the native library with its
assembly backend, same machine, same run:

| | native | wasm | |
|---|---|---|---|
| sign | 1.77 ms | 7.41 ms | 4.2x |
| verify | 3.43 ms | 18.56 ms | 5.4x |

The reason is structural, not a missing optimisation. `wasm32` has no
64×64 → 128 bit multiply, so every field multiplication in the Montgomery
inner loop goes through a software helper, and a pairing is tens of thousands
of those. Native x86-64 and AArch64 use one instruction.

That is fine for verifying the occasional signature or for a demo. It is not a
basis for a throughput claim, and the numbers are published here rather than
omitted for exactly that reason.

## Reference

| | |
|---|---|
| `load()` | resolves to the API. Await once; cached after |
| `keygen(ikm?)` | 32-byte secret key. Deterministic if `ikm` is given |
| `skToPk(sk)` | 48-byte public key |
| `pkValidate(pk)` | throws unless `pk` is usable |
| `sign(sk, msg)` | 96-byte signature |
| `verify(pk, msg, sig)` | `undefined`, or throws |
| `aggregate(sigs)` | one 96-byte signature |
| `aggregateVerify(pks, msgs, sig)` | `undefined`, or throws. Distinct messages |
| `fastAggregateVerify(pks, pops, msg, sig)` | `undefined`, or throws. One message |
| `popProve(sk)` / `popVerify(pk, proof)` | proof of possession |

Everything takes and returns `Uint8Array`. Errors: `ElipsError` (base, catch
this to reject), `InvalidSignature`, `InvalidKey`.

Entropy comes from `crypto.getRandomValues` through emscripten's `getentropy`.

## What this is

`draft-irtf-cfrg-bls-signature`, proof-of-possession scheme,
minimal-pubkey-size: public keys in G1, signatures in G2. The variant Ethereum
uses.

Every operation is checked against [py_ecc](https://github.com/ethereum/py_ecc)
on every build, in both directions, by the *same* script that checks the C
library — `cli.mjs` presents the identical interface, so one script holds both
bindings to one bar rather than two that could drift apart.

**It is not audited.**

## Building

```bash
. /path/to/emsdk/emsdk_env.sh
./build.sh
node --test tests/*.test.mjs
```

`dist/` is generated and not tracked, for the same reason no other generated
file in this repository is: a committed copy is a second version of the
library that can disagree with the first.
