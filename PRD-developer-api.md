# PRD — a developer-facing API, language bindings, and a public demo

Status: **proposed**, awaiting implementation.
Supersedes nothing. Amends `MODERNIZATION_PLAN.md` §10.8; see *Scope* below.

---

## 1. The problem, with evidence

ELiPS is a correct, fast, well-tested pairing library that a working developer
cannot use.

**The API is a toolkit, not a product.** 115 public functions across nine
headers, every one of them in terms of `ep_t`, `ep2_t`, `fp12_t` and
`limb_t[ELIPS_ORDER_LIMBS]`. To do the one thing most people come to a pairing
library for — a BLS signature — a caller must already understand Montgomery
form, cofactor clearing, subgroup membership, domain separation tags, and which
of two groups to hash into. `examples/03_bls_signature.c` is 217 lines and opens
with:

> READ THIS BEFORE COPYING ANY OF IT. This is a demonstration of the library's
> API, not a signature implementation you should deploy.

It is right to say so. That is the gap.

**Installation has a native dependency.** The README's first requirement is
GMP. That single line is why there is no `pip install`, no `cargo add`, no
`npm i`, and no browser build. It is also the reason the library cannot be
tried in under a minute, which is the only budget a developer evaluating a
library actually has.

**Nothing is packaged.** No presence on any package index. The only supported
consumption path is `cmake --build` followed by `find_package(ELiPS)`.

## 2. Who this is for

A developer who wants to verify a signature, is competent, and is not a
cryptographer. They should be able to get a correct result without reading a
specification, and should find it hard to get an incorrect one.

Explicitly **not** for: someone implementing a new pairing-based scheme. That
person wants the existing low-level API, which is not changing.

## 3. What has been verified

None of the plan below rests on recollection. Each of these was measured or
executed in this repository before the PRD was written.

**The library is already WASM-portable, except for GMP.** All nine library
sources were compiled for `wasm32-unknown-unknown` with the GMP header stubbed.
Nine of nine compile. The complete set of symbols a WASM link would have to
provide is: the GMP entry points below, `memcpy`, `memset`, `__multi3`, and an
entropy source. No inline assembly leaks, no platform lock-in; the `.S` files
correctly compile to empty objects off their target.

**The GMP dependency is three operations, not a library.**

| Site | Call | What it does |
|---|---|---|
| `src/util/random.c` | `mpn_sec_div_r` | wide random value mod `r` |
| `src/hash/hash_to_curve.c` | `mpn_sec_div_r` | wide hash value mod `p` |
| `src/arith/fp.c` `fp_inv_sec` | `mpn_sec_invert` | the inversion `fp_inv` **already replaced** |
| `src/arith/fp.c` `fp_inv_vartime` | `mpz_invert` | variable-time inversion |

The first two are the *same* operation at two widths. The primary `fp_inv` is
in-tree Bernstein–Yang divsteps and needs nothing.

**The IETF signature scheme interoperates today.** `py_ecc` 8.0.0 — already a
pinned test dependency in `tools/verify/requirements.txt` — implements
`draft-irtf-cfrg-bls-signature`. A signature it produced was verified by this
library, using only the existing public API:

```
decode pk  : ok            decode sig : ok            hash_to_g2 : ok
e(G1,sig) == e(pk,H(m)) : YES
wrong message rejected  : yes
```

So the compressed encodings, the hash-to-G2 ciphersuite DST and the pairing
check all already agree with an independent implementation, byte for byte. The
signature layer is a thin wrapper over verified behaviour, not new cryptography,
and it has an independent oracle that is already in the repository.

## 4. Scope, and the amendment to §10.8

`MODERNIZATION_PLAN.md` §8 and §10.8 both put a BLS signature layer out of
scope, on the reasoning that signing, verification, aggregation and the
ciphersuite registry are a protocol concern rather than a curve concern. That
reasoning stands and is not being overturned. §10.8 also says the layer "can be
built on top without touching this library", and that is exactly what happens
here:

- `elips_arith` — the core. **Unchanged in scope.** No protocol code enters it.
- `elips_bls` — a **new, separate CMake target** that consumes only the public
  headers of `elips_arith`. Deleting it would not affect the core.

§10.8 is amended to read: the layer lives in this repository, outside the core
library, as its own target.

## 5. Non-goals

- Rust bindings. Rust already has `blst` and `arkworks`, and the audience this
  PRD is written for is largely not writing Rust. Revisit on demand.
- Threshold signatures, identity-based encryption, zero-knowledge gadgets.
- Making the low-level API friendlier. It is aimed at a different reader.
- Competing with `blst` on performance in the browser. See *Risks*.
- A security audit. This layer will carry the same "not audited" statement the
  rest of the repository carries.

## 6. Design

### 6.1 One byte-oriented C API

The junior-facing surface is `include/elips/bls.h`. Everything crosses it as
bytes. No curve type, no limb, no Montgomery form, no DST argument.

```c
#define ELIPS_BLS_SK_BYTES   32
#define ELIPS_BLS_PK_BYTES   48
#define ELIPS_BLS_SIG_BYTES  96

int elips_bls_keygen(uint8_t sk[32], const uint8_t *ikm, size_t ikm_len);
int elips_bls_sk_to_pk(uint8_t pk[48], const uint8_t sk[32]);
int elips_bls_sign(uint8_t sig[96], const uint8_t sk[32],
                   const uint8_t *msg, size_t msg_len);
int elips_bls_verify(const uint8_t pk[48], const uint8_t *msg, size_t msg_len,
                     const uint8_t sig[96]);

int elips_bls_aggregate(uint8_t out[96], const uint8_t *sigs, size_t n);
int elips_bls_aggregate_verify(const uint8_t *pks, size_t n,
                               const uint8_t *msgs, const size_t *msg_lens,
                               const uint8_t sig[96]);
int elips_bls_fast_aggregate_verify(const uint8_t *pks, size_t n,
                                    const uint8_t *msg, size_t msg_len,
                                    const uint8_t sig[96]);

int elips_bls_pop_prove(uint8_t proof[96], const uint8_t sk[32]);
int elips_bls_pop_verify(const uint8_t pk[48], const uint8_t proof[96]);

const char *elips_bls_strerror(int code);
```

Every function returns `0` for success and a negative code otherwise. **There
is no boolean return anywhere**: a verification API that returns an `int` a
caller can accidentally treat as truthy is a footgun, and `if (verify(...))`
must not mean "valid".

**Variant: minimal-pubkey-size.** Public keys in G1 (48 bytes), signatures in
G2 (96 bytes). Chosen because it is what `py_ecc` implements, which is what
gives the layer an in-repo independent oracle, and because it is what Ethereum
uses, which makes the interoperability claim worth something. The 48-byte
aggregate of the other variant is not worth giving up the oracle.

**Ciphersuites**, exactly as registered:

- Basic: `BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_`
- Proof of possession: `BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_`

The DST is **not** a parameter. A caller who can pass the wrong DST will.

**Rogue-key safety is structural, not documented.** `fast_aggregate_verify`
over a shared message is unsafe without proof of possession, so the API does
not offer a way to call it without one. Where the current example warns in a
comment, the API refuses.

`elips_bls` is BLS12-381 only. The other two curves have no registered
ciphersuite, so an "IETF BLS signature on BN-462" would interoperate with
nothing and would be a claim the layer cannot make.

### 6.2 Python

`elips` on PyPI. CFFI over the C API above, shipped as a wheel with no native
dependency once GMP is gone.

```python
import elips
sk  = elips.keygen()
pk  = elips.sk_to_pk(sk)
sig = elips.sign(sk, b"hello")
elips.verify(pk, b"hello", sig)     # returns None, raises InvalidSignature
```

Verification **raises** rather than returning a bool, for the same reason the C
API returns a code: `if verify(...)` must not be a silent accept.

### 6.3 WASM and JavaScript

`@elips/bls` on npm. Emscripten build of `elips_arith` + `elips_bls`, wrapped in
a small typed module. Same shape as Python, `Uint8Array` throughout, verification
throws. Entropy from `crypto.getRandomValues`, injected at the `sysrand` seam.

### 6.4 The demo

A single static page: **aggregate signatures**, because it is the one property
of BLS that is immediately legible without any cryptography.

Sign a message with a few keys. Watch the signature count go up and the
*aggregate stay 96 bytes*. Verify the aggregate against all the public keys at
once. Then tamper with one message and watch it fail. The pairing is doing the
work and the user can see what it bought.

It runs entirely in the browser with no backend, which also makes it the
installation story: "here, it already works".

## 7. Milestones

Each one is independently shippable and independently verifiable.

| # | Deliverable | Gate |
|---|---|---|
| 1 | Remove the GMP runtime dependency | All 52 tests pass with GMP linked only into tests; dudect clean; `fp_difftest` still compares against GMP |
| 2 | `elips_bls` target and `bls.h` | Agrees with `py_ecc` on every operation, on every CI build; rejects each malformed input class |
| 3 | Python package | `pip install` from a built wheel on a clean container with no GMP; same cross-check from Python |
| 4 | WASM and JS package | Builds under emscripten; same cross-check under node; published size recorded |
| 5 | Browser demo | Loads and runs with no backend; a tampered message visibly fails |
| 6 | README rewrite | A developer who has not seen the repository reaches a verified signature in under five minutes |

Milestone 1 gates 3, 4 and 5. Milestone 2 gates 3, 4 and 5.

## 8. Success criteria

These are checks, not aspirations. Each is a CI job or it does not count.

1. `git grep -l gmp.h src/` returns nothing. GMP appears only under `test/` and
   `tools/`.
2. Every `elips_bls` operation agrees with `py_ecc` on randomized inputs, in CI,
   on every build.
3. Each rejection path has a test that *fails* if the check is removed. A
   verifier that cannot reject is not a verifier; this repository already
   applies that rule to its vector runners and it applies here.
4. `pip install` on a container with no GMP and no compiler, then a verified
   signature, in one script.
5. The demo page loads and verifies an aggregate with no network calls after
   load.
6. The low-level API is unchanged: all 52 existing tests still pass, untouched.

## 9. Risks, stated plainly

**WASM will be slow, and the README will say so.** `wasm32` has no 64x64→128
multiply, so `fp_mul` goes through `__multi3`. A pairing that takes ~1.5 ms
natively will be materially slower in the browser. This is fine for a demo and
must not be quietly presented as production performance.

**Removing GMP touches constant-time code.** The wide reduction is used on a
secret message in `hash_to_curve`. Mitigation: GMP stays as a *test* dependency
precisely so `fp_difftest` can keep checking the replacement against it, plus
dudect on the new routine, plus the existing RFC 9380 vectors, which would break
loudly if the reduction were wrong.

**A signature API invites production use.** The layer implements a real
ciphersuite and will interoperate, which makes "this is a demo" a weaker
defence than it is for the examples. Mitigation: say what is and is not true —
implemented to the draft, cross-checked against an independent implementation,
constant-time where it must be, **and not audited** — rather than discouraging
use with a warning nobody reads.

**Scope.** Six milestones is a lot. They are ordered so that stopping after any
one of them still leaves the repository better than it was.
