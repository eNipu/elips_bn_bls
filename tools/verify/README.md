# Independent verification

Everything in `test/` compares the C against `tools/reference/`, and both were
written by the same author. A shared misconception passes both in silence.

That is not a hypothetical. Issue #16 was a final exponentiation raising to the
wrong exponent. It survived for years because the result was still bilinear,
still non-degenerate, and agreed with itself perfectly.

This directory contains checks that share nothing with the library.

## Running it

```sh
tools/verify/run.sh              # everything
tools/verify/run.sh params       # PARI curve parameters
tools/verify/run.sh proofs       # symbolic proofs
tools/verify/run.sh crosscheck   # py_ecc pairing
```

A missing tool is reported as SKIPPED, not as a pass. The exit status is
non-zero if anything that ran, failed.

**In a clean container**, which is what CI uses and what removes any dependence
on what happens to be installed locally:

```sh
tools/verify/run-docker.sh
```

**Locally without Docker** you need PARI/GP, and `run.sh` will build a venv for
`py_ecc` on its own:

```sh
apt-get install --no-install-recommends pari-gp     # Debian/Ubuntu
brew install pari                                   # macOS
```

## What is here

### `curve_params.gp` — PARI counts the curves itself

PARI/GP is the number-theory engine SageMath itself calls. It does not take the
group order from the CM construction formula; it counts the points.

Covers all three curves. That matters because **BLS12-461 and BN-462 have no
published standard anywhere**: before this, their entire correctness rested on
two artefacts by one author agreeing with each other.

Verifies per curve: `p` and `r` prime; `#E(Fp)` equals `p + 1 - t` by PARI's own
count; the G1 cofactor; `#E'(Fp2)` for the twist and the G2 cofactor; the psi
and Frobenius multiplier `p mod r`; and for BN, the reduced GLV lattice basis
and its determinant.

The BN G2 cofactor line is load-bearing: `h2 = p + t - 1` is the identity the
fast cofactor chain in `hash_to_curve.c` depends on. If it were wrong,
`hash_to_g2` would land outside G2 and every downstream pairing would be wrong.

### `identities.gp` — proofs, not tests

`curve_params.gp` checks three specific curves. This proves the same claims as
identities in `Z[X]`, true for **every seed of each family**, including seeds
nobody has chosen yet.

A test says "it held for BN-462". A polynomial identity says "it holds".

Two things worth knowing about this file:

- An earlier draft verified the CM equation with a *power-series* square root.
  That checks an identity only to the truncation order and is not a proof,
  though it prints exactly like one. It uses `issquare()` over `Q[X]` now.
- Both `.gp` files assert the **number of checks that ran**. A statement that
  fails to parse never executes and never counts, so without that guard a
  syntax error reads as a clean pass. This guard has already caught exactly
  that here.

### `crosscheck_pyecc.py` — a pairing computed by someone else

`py_ecc` is the Ethereum Foundation's BLS12-381 implementation. It shares no
code, no author and no representation: its Fp12 is a direct degree-12
extension, not a tower.

The two are reconcilable because ELiPS's tower gives `w^6 = 1+u`, `u^2 = -1`,
hence `w^12 - 2w^6 + 2 = 0`, which is exactly py_ecc's modulus polynomial.
`identities.gp` proves that claim symbolically; the script does the basis change
and compares all twelve coefficients.

ELiPS's value comes out as py_ecc's **inverse**. That is a pure orientation
convention and both are valid non-degenerate pairings, so the script asserts the
inversion rather than tolerating a mismatch.

It also checks that it **rejects** `e^2`, `e^3` and `e^-3`. A cross-check that
cannot fail proves nothing, and `e^3` is precisely what issue #16 computed.

`dump_pairing.c` is the C side. It prints the *exact* pairing from
`pairing_final_exp_plain`, not `elips_pairing`, which returns `e^3` on BLS12 by
design. It prints the generators too: agreeing on an output while disagreeing on
the input proves nothing.

## What none of this does

It verifies the **mathematics the library is supposed to implement**, not the C
that implements it. A carry bug in `fp_mul` is invisible to every check here.

The layers catch different things, and do not substitute for each other:

| failure | caught by |
|---|---|
| wrong exponent, wrong cofactor, wrong lattice basis | this directory |
| wrong field arithmetic | `test/`, the KAT vectors, sanitizers |
| timing dependence on secrets | `dudect` targets |

Issue #16 was an algorithm-level bug: every field operation was correct. A
formally verified field library such as fiat-crypto would not have caught it.
This directory would have, on day one.

Only BLS12-381 gets a third-party *value* check, because py_ecc implements
nothing else. The other two curves get parameters and proofs, which is a real
anchor but a weaker one than an independent pairing value would be.
