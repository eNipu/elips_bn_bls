# ELiPS BN/BLS12 — Modernization Plan

**Status:** proposal for human review. No library code has been changed.
**Author:** analysis pass, 2026-09-05.
**Method:** Karpathy guidelines — state assumptions before coding, prefer the
simplest thing that works, keep changes surgical, and turn every task into a
goal with a mechanical pass/fail check.

---

## 0. Executive summary

The library computes correct, non-degenerate, bilinear pairings on a 462-bit BN
curve and a 461-bit BLS12 curve. It builds clean and its bilinearity tests pass.
That is the good news, and it is worth more than it sounds: the hard part
(tower construction, Frobenius constants, sparse multiplication, optimal-ate
loop parameters, final exponentiation addition chains) is *right*.

The problems are almost entirely below the mathematics:

1. **Representation.** `Fp` is a single heap-allocated `mpz_t`. An `Fp12`
   temporary is twelve separate heap allocations. The Miller loop allocates and
   frees thousands of times per pairing. This, not multiplication throughput, is
   where the time goes.
2. **Algorithmic choices.** Affine coordinates with a full modular inversion per
   point double/add; `mpz_mod` (generic division) after every multiply; no
   Montgomery form; no cyclotomic squaring in the final exponentiation.
3. **A short list of real defects** — memory errors, aliasing bugs, and two
   subroutines whose mathematics is simply wrong (Section 4).
4. **Repository state.** 1280 of 1624 tracked files are generated Doxygen HTML.
   Generated autotools output, a release tarball, built dylibs, dependency
   stamps and fontconfig caches are all committed.

Measured baseline on this machine (Apple Silicon arm64, `clang -O2`, GMP 6.3.0):

| Operation | BN-462 | BLS12-461 |
|---|---|---|
| Miller loop (opt-ate) | 3.85 ms | 2.49 ms |
| Final exponentiation | 5.07 ms | 5.92 ms |
| **Full pairing** | **~8.9 ms** | **~8.4 ms** |
| G1 scalar mult (plain / 2-split) | — | 1.42 / 0.82 ms |
| G2 scalar mult (plain / 2-split / 4-split) | — | 2.97 / 1.70 / 1.03 ms |
| G3 exponentiation (plain / 2-split / 4-split) | — | 7.68 / 5.99 / 3.49 ms |

For scale, a well-tuned BLS12-381 pairing is on the order of one millisecond on
comparable hardware. The gap is roughly an order of magnitude, and it is
addressable without touching the mathematics.

**Expected outcome of this plan:** 8–15x on the pairing, obtained mostly in
Phase 3 and Phase 4, with assembly (Phase 6) contributing the last ~1.3x and
deliberately deferred until measurement proves it is the bottleneck.

---

## 1. Assumptions, and where I need a decision

Karpathy rule #1: do not pick silently between interpretations. Six items need
a human answer. Everything else I will decide myself and record.

### 1.1 "Comparable to RELIC" needs a definition

RELIC does not ship the parameters this library uses. Its presets cover BN-254,
BN-382, BN-446, BLS12-381, BLS12-446, BLS12-455, BLS12-638. There is no BN-462
and no BLS12-461. So "comparable to RELIC" cannot be a direct A/B run as-is.

Three options:

| Option | Meaning | Cost |
|---|---|---|
| **A** | Keep BN-462/BLS12-461; compare against RELIC BLS12-455 and scale | cheap, approximate |
| **B** | Add BLS12-381 and BN-446 to *our* library, compare head-to-head | ~1 phase of work, exact |
| **C** | Add BN-462 to RELIC as a custom prime, compare head-to-head | fiddly, exact |

**Recommendation: B.** BLS12-381 is the parameter set the rest of the world
actually uses, it makes us comparable to RELIC, blst and arkworks at once, and
it gives us third-party known-answer vectors for free — which Phase 0 needs
anyway. BN-462 and BLS12-461 stay supported.

### 1.2 Which CPU are we optimizing for?

This machine is arm64. RELIC's fast paths are x86-64 (`x64-asm-4l` through
`x64-asm-12l`); its ARM assembly backend is 32-bit (`arm-asm-254`). So "match
RELIC's assembly" on this hardware means **writing AArch64 assembly that does
not exist upstream**, not porting theirs.

Decision needed: primary target is (a) AArch64, (b) x86-64, or (c) both. This
changes Phase 6 substantially. Default if unanswered: AArch64 first (it is the
development machine), x86-64 second, with a portable C fallback always present.

### 1.3 "Fast MPZ" — what I think you mean

There is no library called "fast MPZ". Reading it as *make the multiprecision
layer fast*, the actual wins, in descending order of payoff:

1. Stop using `mpz_t` for field elements. A 462-bit element is 8 limbs on a
   64-bit machine — put it in a fixed-size array inside the struct, on the
   stack, with no allocation.
2. Drop from `mpz_*` to GMP's `mpn_*` layer. `mpn_*` operates on caller-supplied
   limb arrays with no allocation, no size negotiation, and no sign handling.
3. Montgomery representation, so reduction is `mpn_*` shifts and multiply-adds
   instead of `mpz_mod`'s division.

This is what RELIC does (`dig_t fp_t[RLC_FP_DIGS]`, `INTEG` integrated
multiply-and-reduce). Confirm that this is the intent; it is the single largest
change in the plan and it rewrites every type in the library.

### 1.4 Bazel

You named Bazel. I will build Bazel. One honest note first, then I stop arguing:

GMP is autotools, RELIC is CMake, and every downstream consumer of a C crypto
library expects `pkg-config` or CMake. Bazel will need `rules_foreign_cc` to
build GMP from source, or a `cc_library` shim over a system GMP. It is perhaps
a day of extra plumbing versus CMake, and it makes the library slightly harder
for outsiders to consume.

Bazel does buy real things here: hermetic builds, a genuinely good test runner
with sharding, `--config=asan/ubsan/msan` as one-liners, and clean multi-toolchain
cross-compilation for the assembly phase. Given the assembly work coming in
Phase 6, that last point is worth the plumbing.

**Plan of record: Bazel primary, with a thin CMake file kept for consumers.**
Say the word if you want CMake-only.

### 1.5 Dead and near-dead code — delete or port?

| Component | State | Recommendation |
|---|---|---|
| KSS16 (`Fp4`, `Fp8`, `*_KSS16`) | incomplete, and mathematically broken (§4.2) | **delete** |
| `Fp16`, `EFp4`, `EFp8`, `EFp16` types | declared, no implementation file exists | **delete** |
| `EFp6` / `Fp6_sqrt` | `Fp6_sqrt` has its core commented out | **delete** |
| Tate, plain-ate, x-ate pairings | work, but superseded by opt-ate | **keep**, they are cross-checks |
| Timing/print scaffolding | welded into the library API | **extract** to a bench target |

Deleting the KSS16 path removes ~700 lines and one whole class of latent bugs.
It is not used by any pairing.

### 1.6 Constant-time: a requirement or not?

Currently every scalar multiplication branches on secret bits, and GMP's
`mpz_invert` and `mpz_mul` are variable-time. If this library is for research
and benchmarking, that is fine and Phase 5 shrinks a lot. If anything will ever
sign with a secret key, it is not fine.

Default if unanswered: implement constant-time scalar multiplication and
inversion, but keep the variable-time paths available behind an explicit
`*_vartime` suffix for public inputs, which is where the speed is.

---

## 2. What the code actually is

### 2.1 Shape

12,523 lines across 47 `.c` files and 50 headers. Two curve families (BN, BLS12)
plus an unfinished KSS16. A `src/` + `include/ELiPS_bn_bls/` split, one header
per source file, autotools build, one hand-driven `test/main.c`.

### 2.2 The tower

```
Fp                       p = 462 bits (BN) / 461 bits (BLS12)
Fp2  = Fp[u]/(u² + 1)                       u² = −1
Fp6  = Fp2[v]/(v³ − (1+u))                  cubic non-residue 1+u
Fp12 = Fp6[w]/(w² − v)
```

This is a standard, sound tower. `p ≡ 3 (mod 4)` holds for both curves, so
`u² = −1` is irreducible. `p ≡ 1 (mod 3)` holds, so the cube roots of unity used
for the Frobenius constants exist. I verified both.

Struct layout is nested by value — `Fp12` contains two `Fp6` by value, each
containing three `Fp2`, each containing two `Fp`, each containing one `mpz_t`.
So an `Fp12` is contiguous *except* that each of its 12 `mpz_t` points at a
separately malloc'd limb array. That is the core performance problem, and it is
also the reason the fix is mechanical: the nesting is already right, only the
leaf changes.

### 2.3 Verified-correct mathematics

I checked these against the defining equations and they are right. Listing them
matters as much as listing the bugs, because they constrain what the rewrite may
not disturb.

- **BN parameters.** `x = 2¹¹⁴ + 2¹⁰¹ − 2¹⁴ − 1`;
  `p = 36x⁴+36x³+24x²+6x+1` (462 bits), `r = 36x⁴+36x³+18x²+6x+1` (462 bits),
  `t = 6x²+1`. I confirmed `p + 1 − t = r` exactly, so the cofactor is 1.
- **BLS12 parameters.** `X = −2⁷⁷ + 2⁵⁰ + 2³³`;
  `p = ((X−1)²(X⁴−X²+1))/3 + X` (461 bits), `r = X⁴−X²+1` (308 bits).
  I confirmed the division by 3 is exact and that `r | p+1−t`, cofactor 153 bits.
- **`X_binary_opt` encodes exactly `6x+2`.** I summed the signed digits
  `{116:+1, 115:+1, 103:+1, 102:+1, 16:−1, 15:−1, 2:−1}` and got `6x+2` on the
  nose. This is the correct optimal-ate loop parameter for BN.
- **`weil()` computes `#E(F_p¹²) = p¹² − t₁₂ + 1` correctly** via the
  `t₂ → t₆ → t₁₂` Weil-sum recurrence. Only the *variable name* is wrong: it is
  stored in `EFpd_total`, which the header documents as the twisted curve's
  order over `F_p^{k/d}`. Misleading, not incorrect.

### 2.4 Two things that look like bugs and are not

I want these on the record, because a careless refactor will "fix" them and
break the library. I chased both and they are correct.

**The Miller loop bounds are right.** BN iterates `i = bn_X_length+1 = 115` down
to 0, while `6x+2` has its top digit at index **116**. That looks off by one. It
is not: `T` is initialized to `Q`, which *is* the leading `+1` digit — the
standard "seed the accumulator with the most significant digit, then loop over
the rest" pattern. Same for BLS12: it iterates from `bls12_X_length−1 = 76`
while `X`'s top digit sits at 77 with value `−1`, and `T` is initialized to
`−Q`. Both correct.

**`Fp2_set_ui` setting *both* coefficients is deliberate.** `Fp2_set_ui(&a, 1)`
yields `1 + u`, not the field element 1, which reads as a bug until you find
`set_basis()`, which relies on exactly that to construct the non-residue `1+u`.
It is a bad name for an intentional "broadcast to all coefficients" operation.
It stays a hazard (§4.4) but it is not a live defect.

---

## 3. Where the time goes

Ranked by expected payoff, from the measured 8.9 ms BN pairing.

**1. Heap allocation, ~40–50% of runtime.** `Fp12_init` is twelve `mpz_init`
calls. `ff_ltt` — called once per Miller iteration, 116 times per pairing —
constructs two `Fp12` and five `Fp2` temporaries, so roughly 34 mallocs and 34
frees per iteration, ~4000 per pairing, before any arithmetic happens.

**2. Modular inversion in the inner loop, ~20–25%.** The library uses affine
coordinates. `EFp2_ECD` and `ff_ltt` each perform a full `Fp2_inv`, which is a
`mpz_invert` plus multiplications. RELIC uses projective/Jacobian coordinates
and pays exactly one inversion per pairing, at the end. This is the single
biggest *algorithmic* difference between the two libraries.

**3. `mpz_mod` after every multiply, ~15–20%.** `Fp_mul` is `mpz_mul` then
`mpz_mod`. `mpz_mod` is truncating division. Montgomery reduction replaces it
with a multiply-and-shift; for a 462-bit modulus that is roughly a 2x
improvement on the reduction step alone, and it composes with lazy reduction
(reduce once per `Fp2`/`Fp6` operation instead of per `Fp` operation).

**4. No cyclotomic squaring, ~10% of the final exponentiation.** After the easy
part `f^(p⁶−1)(p²+1)`, the element lies in the cyclotomic subgroup, where
Granger–Scott compressed squaring costs substantially less than a generic
`Fp12_sqr`. The final exponentiation is 5–6 ms, the larger half of the pairing,
and it is currently all generic squarings. RELIC's `fp12_sqr_cyc` does this.

**5. Duplication.** `bn_p8sparse.c` and `bls12_p8sparse.c` are byte-identical
apart from function names. Same for `line_ate`, `frobenius`, `skew_frobenius`
and `twist` — the Frobenius pair differs in 21 lines, all of them names and
comments. Roughly 40% of `src/` is copy-paste. This is not a performance issue
but it doubles the cost of every subsequent fix and doubles the surface for
divergence bugs.

---

## 4. Correctness audit

This is the section you asked me to be hardest on. Each item was confirmed by
reading the code and, where the claim is mathematical, by independent
computation. Severity: **C**ritical (wrong results or UB on a live path),
**H**igh (UB or leak on a live path), **M**edium (latent).

### 4.1 Memory and undefined behaviour

| # | Sev | Location | Defect |
|---|---|---|---|
| M1 | **H** | `curve_settings.c:382,394,403` | `bls12_generate_prime()` calls `mpz_init` where it must call `mpz_clear` — on all three exit paths. Re-initializing a live `mpz_t` orphans its limb array. Four leaked integers per call. |
| M2 | **H** | `curve_settings.c:282,468` | `curve_parameters.EFpd_total` is **never `mpz_init`'d** — not in `init_bn_parameters()`, not in `init_bls12_parameters()` — yet `weil()` writes it and both G2 generators read it. It only works because the struct is a zero-initialized global and GMP tolerates `{alloc=0, size=0, d=NULL}`. Formally undefined. |
| M3 | **H** | `curve_settings.c:72,94` vs `bn_clears.c:54` | `curve_parameters.curve_a` is left uninitialized for BN and BLS12 (both `mpz_init` calls are commented out), then read by `mpz_cmp_ui` at `curve_settings.c:315` and passed to `mpz_clear` at teardown. `mpz_clear` on an uninitialized `mpz_t`. |
| M4 | **H** | `bn_line_ate.c:90`, `bls12_line_ate.c:90` | `f_ltq()` initializes `Fp2 E` and never clears it. This is inside the Miller loop: it leaks once per non-zero digit per pairing. |
| M5 | **M** | `bls12_finalexp.c:70,71` | `t4` and `t5` are initialized and never cleared; the teardown covers only `t0`–`t3`. |
| M6 | **M** | `bls12_scm.c:333` | Cleanup calls `EFp2_init(&twisted_Q_3x)` where it means `EFp2_clear`. Leak plus double-init. |
| M7 | **H** | 14 sites | `int length = mpz_sizeinbase(s,2); char binary[length]; mpz_get_str(binary,2,s);` — GMP writes `length` digits **plus a NUL terminator**, so every one of these overflows its stack buffer by one byte. Sites: `bn_fp.c:206`, `bn_fp2.c:416`, `bn_fp6.c:347`, `bn_fp12.c:275`, `bn_efp.c:223`, `bn_efp2.c:199`, `bn_efp6.c:195`, `bn_efp12.c:220`, `bn_miller_ate.c:47`, `bn_miller_tate.c:39`, `bls12_miller_ate.c:65`, `bls12_miller_tate.c:40`, `fp4.c:168`, `fp8.c:175`. Violates CERT **ARR30-C** and **STR31-C**; the VLA also violates **ARR32-C** since the size is derived from a runtime value. |

M7 is the one I would fix today regardless of the rest of the plan.

### 4.2 Wrong mathematics

| # | Sev | Location | Defect |
|---|---|---|---|
| X1 | **C** | `bn_fp.c:66` | `Fp_mul_basis_KSS16` reduces modulo `curve_parameters.curve_a` instead of `curve_parameters.prime`. `curve_a` is uninitialized for BN/BLS12 and set to **1** for KSS16 — so the result is identically zero. Every `Fp4`/`Fp8` operation built on it is wrong. |
| X2 | **C** | `fp4.c:213`, `fp8.c:221` | `Fp4_sqrt` and `Fp8_sqrt` both compute the group order as `p¹²−1`. `Fp4` has order `p⁴`, `Fp8` has order `p⁸`. Tonelli–Shanks over the wrong group does not return a square root. |
| X3 | **C** | `fp4.c:244`, `fp8.c:252` | Same two functions compute the Tonelli–Shanks step as `mpz_powm(t, 2, r−m−1, p)` — reducing an **exponent** modulo `p`. It must be the integer `2^(r−m−1)`. Two independent errors in the same routine. |
| X4 | **H** | `bn_fp6.c:315,316,319,320` | `Fp6_sqrt` has its Frobenius-constant multiplications commented out, so `tmp1` and `tmp2` are built from `x0` alone and the rest is garbage. The function returns a value that is not a square root. |
| X5 | **M** | `bn_fp2.c:163` vs `:218` | `Fp2_sqr_KSS16` calls `Fp_mul_basis` (the BN basis) while the otherwise-identical `Fp2_sqr_kss16` calls `Fp_mul_basis_KSS16`. Two functions differing only in case and in one call; at most one can be right. |

X1–X3 are confined to the KSS16/`Fp4`/`Fp8` path, which no pairing uses. That is
why the tests pass. Under §1.5 this code gets deleted rather than fixed, which
resolves all four at once — but if any of it is wanted, it needs rewriting from
the reference algorithm, not patching.

### 4.3 API contract violations

| # | Sev | Location | Defect |
|---|---|---|---|
| A1 | **C** | `bn_fp2.c:114,115` | `Fp2_mul_Fp(ANS, A, B)` computes `ANS = ANS * B`, ignoring `A` entirely. Correct only when the caller happens to alias `ANS == A`. |
| A2 | **H** | `bn_fp12.c:143,153,158` | `Fp12_add_mpz`, `Fp12_sub_ui` and `Fp12_sub_mpz` all read `ANS` where they should read `A`. Same defect class as A1. |
| A3 | **H** | `bn_final_exp.c:41,44,70,74`; `bls12_finalexp.c:77,81` | The final exponentiation **writes through its input pointer** (`Fp12_mul(A, &t0, &t1)`). `bn12_opt_ate` calls it as `bn_final_exp_optimal(ANS, ANS)` so the damage is invisible, but the published signature promises `A` is an input. Calling it twice on one value gives two different answers. |
| A4 | **H** | `bn_final_exp.c:87–89`; `bls12_finalexp.c:88–102` | Both final exponentiations **mutate global curve parameters mid-computation** and restore them afterwards — `X_binary[0]` for BN, and `bls12_X_length` plus three entries of `bls12_X_binary[]` for BLS12. The library is therefore not reentrant and not thread-safe, and any early return leaves the curve parameters corrupted for the whole process. |
| A5 | **M** | `bn_fp.c:99,100` | `Fp_inv` discards `mpz_invert`'s return value. On a non-invertible input the output is left undefined and the caller is never told. CERT **EXP12-C**. |
| A6 | **M** | `curve_settings.c` | `generate_bn_prime`/`generate_bn_order` return 1 for success; `bls12_generate_order` returns 0 for success. Every caller ignores the value, so `init_bls12_settings()` proceeds cheerfully when prime generation has failed. CERT **ERR33-C**. |
| A7 | **L** | `curve_settings.c:488` | `gmp_printf("E:y^2=x^3+4\n", curve_parameters.curve_b)` — an argument with no conversion specifier. CERT **FIO47-C**. |
| A8 | **L** | `curve_settings.h`, `bls12_timeprint.h` | `bls12_print_parameters` is declared in two headers. |

A3 and A4 together are why I rate the final exponentiation as the least
trustworthy code in the library despite producing correct answers today.

### 4.4 Design hazards

- **`Fp2_set_ui` / `Fp2_cmp_one` disagree.** `set_ui(&a,1)` produces `1+u`;
  `cmp_one(&a)` tests for `1+0u`. So `set_ui(&a,1); cmp_one(&a)` reports "not
  one". Both are individually defensible, together they are a trap. Rename to
  `Fp2_set_all_ui` and add a real `Fp2_set_one`.
- **`Fp_cmp` returns 0 or 1, never negative**, despite a name that implies
  `memcmp`-style ordering. `fp4.c` then writes `while (Fp4_cmp_mpz(&b,set_1)==1)`
  meaning "while not equal", which reads as its own opposite.
- **Global mutable state as algorithm input.** `curve_parameters`, `X_binary`,
  `bls12_X_binary`, `epsilon1/2`, `d12_frobenius_constant[][]`, the RNG `state`,
  and the timing globals `t0`/`t1` are all process-wide. Two curves cannot be
  used at once; `init_bn()` and `bls12_inits()` overwrite each other.
- **Benchmarking is welded into the API.** `bls12_plain_G1_scm` calls
  `gettimeofday` and writes the global `bls12_G1SCM_PLAIN`. Measurement belongs
  in a bench harness, not in the function under measurement.
- **RNG.** Every `*_sqrt` and `*_rational_point` builds a *local* `gmp_randstate_t`
  shadowing the global one and seeds it with `time(NULL)` — one-second
  granularity, entirely predictable. Fine for tests, unacceptable for key
  material. CERT **MSC32-C**.

### 4.5 What the current tests can and cannot catch

This matters more than any single bug above.

Every test has the same shape: draw random `s₁`, `s₂`, compute `e(P,Q)^(s₁s₂)`,
`e([s₁]P,[s₂]Q)` and `e([s₂]P,[s₁]Q)`, check all three agree and that the result
is neither 0 nor 1.

That verifies **bilinearity and non-degeneracy**. It does not verify that the
map is *the pairing it claims to be*. A truncated Miller loop, a wrong Frobenius
constant, or a botched final exponent can all yield a map that is still
bilinear, still non-degenerate, and still passes all three checks. The suite
would report success.

There are also no known-answer vectors, no `[r]P = O` subgroup checks, no
`assert`s, no exit code — `main` returns 0 unconditionally, so nothing is
CI-gradeable. Fixing this is Phase 0, and it gates everything else.

---

## 5. Target architecture

```
bazel/                    toolchains, GMP external, sanitizer configs
include/elips/
  config.h                generated: limb width, curve selection
  fp.h  fp2.h  fp6.h  fp12.h
  ep.h  ep2.h             elliptic curves over Fp and Fp2
  pairing.h               public API
  params.h                curve parameter descriptors
src/
  arith/
    fp_mont.c             Montgomery mul/sqr/reduce over mpn
    fp_inv.c              constant-time inversion
    fp_asm_a64.S          Phase 6, AArch64
    fp_asm_x64.S          Phase 6, x86-64
  tower/
    fp2.c fp6.c fp12.c    lazy reduction, cyclotomic ops
  curve/
    ep.c ep2.c            projective coordinates
    params_bn462.c params_bls12_461.c params_bls12_381.c
  pairing/
    miller.c              ONE implementation, parameterized
    line.c sparse.c
    final_exp.c           per-family addition chains
test/
  kat/                    known-answer vectors + the Sage script that made them
  unit/                   field-axiom and property tests
  vs_relic/               differential tests against RELIC
bench/                    all timing code lives here and nowhere else
```

Five principles, each answering a specific defect above:

1. **One implementation per algorithm.** The BN and BLS12 Miller loops differ
   only in the loop parameter and a final Frobenius correction. They become one
   function taking a `curve_params *`. This deletes ~40% of `src/`.
2. **No global state.** Curve parameters travel in a context struct passed
   explicitly. Fixes §4.4 and makes two curves usable in one process.
3. **Inputs are `const`.** Fixes A1, A2, A3 by making them compile errors.
4. **Stack allocation.** `typedef limb_t fp_t[FP_LIMBS]` — no malloc on any hot
   path. Fixes performance item 1 and, incidentally, M4/M5/M6, because there is
   nothing left to leak.
5. **Every function returns a status, or provably cannot fail.** Fixes A5, A6.

---

## 6. Phased plan

Each phase has a mechanical exit gate. If the gate does not pass, the phase is
not done. No phase begins before its predecessor's gate is green.

### Phase 0 — Build the oracle. Change nothing else.

The most important phase, and the one it is most tempting to skip. Right now
there is no way to tell a correct refactor from an incorrect one.

1. Tag the current commit as the reference implementation.
2. Write `test/kat/gen_vectors.sage` producing, for each curve:
   field arithmetic vectors (mul/sqr/inv/Frobenius at every tower level),
   curve vectors (double/add/scalar-mult with known scalars), and full pairing
   vectors — `e(P,Q)` for fixed generators, printed as 12 `Fp` coordinates.
   Sage is present at `/usr/local/bin/sage`; it computes these independently, so
   it is a genuine oracle and not a restatement of our own code.
3. For BLS12-381 (per §1.1 option B), cross-check the Sage vectors against
   published test vectors before trusting them.
4. Add real property tests: `[r]P = O` for G1 and G2; `e(P,Q)^r = 1`;
   `e(P₁+P₂,Q) = e(P₁,Q)·e(P₂,Q)`; the field axioms at every tower level.
5. Make the runner exit non-zero on failure.
6. Record the baseline timings in `bench/baseline.json`.

**Gate:** the *unmodified* library passes every generated vector, and the runner
returns non-zero when a vector is deliberately corrupted.

If the current library fails a Phase 0 vector, that is a finding, and it changes
the plan. I consider this maybe 20% likely — most plausibly in the final
exponentiation, given A3/A4.

### Phase 1 — Repository hygiene and Bazel

1. `git rm -r --cached` the 1280 generated Doxygen files, autotools output
   (`configure`, `Makefile`, `config.status`, `libtool`, `aclocal.m4`,
   `src/.deps/`), `src/.libs/`, `elipsbnbls12-1.0.1.tar.gz`, `fontconfig/`,
   `test/a.out`. A real `.gitignore` replaces the current two-line one.
   History is preserved; only the working tree stops carrying them.
2. `MODULE.bazel` with `rules_foreign_cc` building GMP hermetically.
3. `cc_library` for the library, `cc_test` per test group.
4. `.bazelrc` configs: `asan`, `ubsan`, `msan`, `opt`.
5. Keep a thin `CMakeLists.txt` for downstream consumers.
6. CI: build, test, and all three sanitizers on every push.

**Gate:** `bazel test //...` green from a clean checkout on a machine with no
GMP installed. Phase 0's vectors still pass, byte for byte.

Note: the sanitizer run will immediately surface M1–M7. That is intended — it
is Phase 1 proving the tooling works.

### Phase 2 — Fix the confirmed defects, under test, on the existing architecture

Deliberately *before* the rewrite, so each fix is a small reviewable diff
against a passing oracle rather than being smuggled into a 10,000-line change.

Order: M7 (stack overflows) → M1, M2, M3 (initialization) → A1, A2 (aliasing) →
M4, M5, M6 (leaks) → A5, A6, A7, A8 → delete the KSS16/`Fp4`/`Fp8`/`Fp16`/`EFp6`
tree, which retires X1–X5 wholesale.

Each fix is TDD: write the test that fails, then the fix, then watch it pass.
M7 gets a test compiled with ASan that currently reports a stack-buffer-overflow.

A3 and A4 are explicitly **not** fixed here — they are structural, and Phase 3
deletes the globals that cause them.

**Gate:** all Phase 0 vectors pass; ASan, UBSan and MSan all clean; measured
timings within noise of baseline (this phase is not supposed to change speed).

### Phase 3 — The representation change

The big one, and the source of most of the speedup.

1. `typedef uint64_t limb_t; typedef limb_t fp_t[FP_LIMBS];` sized from the
   prime at configure time. No allocation anywhere in `fp_t`.
2. Montgomery form: `to_mont`, `from_mont`, CIOS multiply-and-reduce over
   `mpn_*`. Conversion happens at the API boundary only.
3. Rewrite `fp2`/`fp6`/`fp12` on `fp_t`, with **lazy reduction** — accumulate
   double-width intermediates through an `Fp2` multiplication and reduce once at
   the end rather than after each `Fp` operation.
4. Curve parameters move into a context struct; the globals go away, taking A3
   and A4 with them.
5. Miller loops merge into one parameterized implementation.

Do this **one tower level at a time**, keeping the old implementation beside the
new one and differential-testing them against each other on random inputs until
the level is proven, then delete the old.

**Gate:** every Phase 0 vector passes; the differential test finds no
disagreement across 10⁶ random inputs per level; sanitizers clean; measured
speedup ≥ 3x on the full pairing. If it is under 3x, stop and re-measure before
continuing — the model of where time goes is wrong and the rest of the plan
needs revisiting.

### Phase 4 — Algorithmic improvements

1. **Projective coordinates** for `EFp` and `EFp2`, so the Miller loop performs
   zero inversions instead of one per iteration. One inversion at the end.
2. **Cyclotomic squaring** (Granger–Scott) in the final exponentiation and in
   `G3` exponentiation, with a comment stating the precondition — the element
   must lie in the cyclotomic subgroup, which holds after the easy part and
   nowhere else.
3. **Compressed squaring** for the `f^x` chains.
4. **Sparse multiplication** without the full `Fp12` temporary that
   `Pseudo_8_sparse_mul` currently allocates.
5. **Fixed windowing / wNAF** for scalar multiplication, replacing the current
   binary ladder.

Each is landed and measured independently. Anything that does not pay for itself
gets reverted rather than kept out of sentiment.

**Gate:** vectors pass; cumulative speedup ≥ 8x on the full pairing; every
individual optimization is shown to be a measured win.

### Phase 5 — Constant-time and API hardening

Scope depends on §1.6. Assuming yes:

1. Constant-time conditional select and swap; no branch on secret data.
2. Constant-time inversion — `x^(p−2)` via a fixed addition chain, or safegcd.
   GMP's `mpz_invert` is variable-time and must not touch secrets.
3. Montgomery ladder or fixed-window with full-table scan for scalar mult.
4. Subgroup membership checks on every point entering the API. Currently absent
   and required against small-subgroup attacks.
5. Full input validation at the trust boundary — points on curve, in subgroup,
   scalars in range.
6. `dudect` timing-leakage tests in CI.

**Gate:** `dudect` reports no leakage on secret-dependent paths; the variable-time
fast paths are reachable only through explicitly `_vartime`-suffixed functions.

### Phase 6 — Assembly

**Only if Phase 4's profile shows limb arithmetic is the bottleneck.** Ponytail
rule: do not write assembly to make malloc faster. My expectation is that after
Phases 3–4 the remaining hot spot genuinely is `mpn_mul`/Montgomery reduction,
in which case:

1. AArch64 first (§1.2): 8-limb Montgomery multiply using `MUL`/`UMULH` with
   `ADCS` carry chains. `fp_asm_a64.S`.
2. x86-64 with BMI2/ADX (`MULX`, `ADCX`, `ADOX`) — this is where RELIC's
   `x64-asm-*` backends get their advantage.
3. Every assembly routine ships with the C fallback it replaces, and both are
   run against the same vectors on every CI build.
4. Runtime CPU feature detection with a portable default.

**Gate:** assembly and C paths produce bit-identical output on all vectors;
measured win ≥ 1.2x, otherwise the assembly is deleted rather than carried.

---

## 7. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Phase 0 finds the current library is already wrong somewhere | ~20% | Better to know now. Fix before refactoring, not during. |
| Montgomery conversion introduces a subtle sign/carry bug | medium | Differential test against the mpz implementation, 10⁶ random inputs per level, before deleting the old path. |
| Lazy reduction overflows the accumulator | medium | Prove the bound explicitly for each modulus and assert it in debug builds. This is the classic lazy-reduction failure. |
| Bazel + GMP plumbing eats a week | medium | Timebox it. Fall back to system GMP via a `cc_library` shim; §1.4 is reversible. |
| Assembly wins nothing after Phases 3–4 | medium | Phase 6 is explicitly conditional on measurement. |
| Scope grows without end | **high** | This document is the scope. Anything not in it is a separate proposal. |

---

## 8. Explicitly not doing

- Adding curve families beyond BN, BLS12, and (per §1.1) BLS12-381.
- Reviving KSS16 — see §1.5 and §4.2.
- A protocol layer (BLS signatures, IBE). Different project.
- GPU or multi-threaded pairings.
- Keeping API compatibility with the current headers. The current API cannot
  express `const` inputs or a curve context, and preserving it would forfeit
  the fixes to A3 and A4.

---

## 9. What I need from you

1. §1.1 — comparison target. Recommendation: add BLS12-381.
2. §1.2 — AArch64, x86-64, or both.
3. §1.3 — confirm "fast MPZ" means fixed-width `mpn` plus Montgomery.
4. §1.4 — Bazel confirmed, or switch to CMake.
5. §1.5 — confirm deletion of KSS16 / `Fp4` / `Fp8` / `Fp16` / `EFp6`.
6. §1.6 — constant-time required, or research-grade acceptable.

With those six answered I will start at Phase 0, which changes no library code
and produces the test oracle everything after it depends on.

---

## Appendix A — Verification log

Claims in this document that were checked by computation rather than by reading:

| Claim | Method | Result |
|---|---|---|
| BN `p`, `r`, `t` polynomials, and `p+1−t = r` | Python, exact integers | confirmed, cofactor 1 |
| BN `p` is 462 bits, `r` is 462 bits | Python | confirmed |
| `X_binary_opt` = `6x+2` | summed signed digits | exact match |
| BN Miller loop bound is correct despite appearances | traced MSB-preload against digit index 116 | **not a bug** |
| BLS12 `p` is 461 bits, `r` is 308 bits, `3 \| (X−1)²(X⁴−X²+1)`, `r \| p+1−t` | Python | confirmed, cofactor 153 bits |
| BLS12 Miller loop bound is correct | traced `T = −Q` preload against digit index 77 | **not a bug** |
| `p ≡ 3 (mod 4)` and `p ≡ 1 (mod 3)` for both curves | Python | confirmed; tower is valid |
| `mpz_get_str` needs `sizeinbase+1` bytes for positive input | GMP documentation plus arithmetic | overflow of exactly 1 byte at 14 sites |
| `Fp2_set_ui` dual-set is relied on by `set_basis()` | read `set_basis`, checked `Fp2_mul_basis` implements `×(1+u)` | **not a bug** |
| Baseline timings | built with `clang -O2`, ran both suites | recorded in §0 |
| RELIC ships no BN-462 or BLS12-461 preset | enumerated `relic/preset/` and `relic_ep.h` | confirmed |

Two hypotheses I formed while reading were **disproved** by this process — the
Miller loop bounds and `Fp2_set_ui`. Both would have caused a regression had I
"fixed" them. This is the argument for Phase 0 in miniature.
