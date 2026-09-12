# elips

BLS signatures on BLS12-381, from the [ELiPS](https://github.com/eNipu/elips_bn_bls)
pairing library.

```bash
pip install elips
```

No system dependencies. The wheel links `libc` and nothing else.

```python
import elips

sk  = elips.keygen()
pk  = elips.sk_to_pk(sk)
sig = elips.sign(sk, b"hello")

elips.verify(pk, b"hello", sig)      # returns None, or raises
```

## Verification raises. It does not return a boolean.

This is the one thing to know before using the package.

```python
try:
    elips.verify(pk, msg, sig)
except elips.ElipsError:
    reject()
```

A verify that returns a value invites `if verify(...): accept()`, which is
wrong in the direction that matters. Here that branch simply never fires, so
the mistake shows up on the first *valid* signature rather than on the first
forged one.

`ValueError` and `TypeError` mean the calling code passed the wrong length or
the wrong type. That is a bug in the caller, not a bad signature, so it is a
different exception and is raised before any cryptography runs.

## Aggregation

`n` signatures combine into one, still 96 bytes, and verify against all `n`
public keys at once.

```python
sigs = [elips.sign(sk, msg) for sk, msg in zip(sks, msgs)]
agg  = elips.aggregate(sigs)          # 96 bytes, whatever n was
elips.aggregate_verify(pks, msgs, agg)
```

Messages must be distinct. A repeated message is refused, because without
proofs of possession an attacker can move a signature between signers. For one
shared message use `fast_aggregate_verify`, which takes the proofs and checks
them:

```python
pops = [elips.pop_prove(sk) for sk in sks]
agg  = elips.aggregate([elips.sign(sk, msg) for sk in sks])
elips.fast_aggregate_verify(pks, pops, msg, agg)
```

There is no way to call it without the proofs. Aggregating keys over a shared
message is forgeable otherwise: an attacker who registers
`pk_evil = [t]G1 - sum(pk_honest)` can produce an aggregate for a group whose
keys they never held.

## Reference

| | |
|---|---|
| `keygen(ikm=None)` | 32-byte secret key. Deterministic if `ikm` is given |
| `sk_to_pk(sk)` | 48-byte public key |
| `pk_validate(pk)` | raise unless `pk` is usable |
| `sign(sk, msg)` | 96-byte signature |
| `verify(pk, msg, sig)` | `None`, or raises |
| `aggregate(sigs)` | one 96-byte signature |
| `aggregate_verify(pks, msgs, sig)` | `None`, or raises. Distinct messages |
| `fast_aggregate_verify(pks, pops, msg, sig)` | `None`, or raises. One message |
| `pop_prove(sk)` / `pop_verify(pk, proof)` | proof of possession |

Exceptions: `ElipsError` (base, catch this to reject), `InvalidSignature`,
`InvalidKey`.

## What this is

`draft-irtf-cfrg-bls-signature`, proof-of-possession scheme,
minimal-pubkey-size: public keys in G1, signatures in G2. The variant Ethereum
uses.

Every operation is checked against
[py_ecc](https://github.com/ethereum/py_ecc), an independent implementation,
on every build — both directions, since an implementation can be
self-consistently wrong and pass a one-way check.

**It is not audited.** Implementing the draft and agreeing with another
implementation is a different claim from having been reviewed by a
cryptographer.

## Building from source

Needs a C compiler and nothing else:

```bash
pip install ./bindings/python
python -m pytest bindings/python/tests
```

`build_ffi.py` copies `include/` and `src/` from the repository root into
`vendor/` on every build, so the package is self-contained in an sdist and
cannot go stale against the C it wraps.
