# ELiPS BN/BLS12 — Modernization Plan (rev. 2)

**Status:** scope settled, ready to start at Phase 0. No library code changed yet.
**Revised:** 2026-09-05, incorporating the six scope decisions.
**Method:** Karpathy guidelines — state assumptions before coding, prefer the
simplest thing that works, keep changes surgical, and turn every task into a
goal with a mechanical pass/fail check.

---

## 0. Executive summary

The library computes correct, non-degenerate, bilinear pairings on a 462-bit BN
curve and a 461-bit BLS12 curve. It builds clean and its bilinearity tests pass.
That matters more than it sounds: the hard part — tower construction, Frobenius
constants, sparse multiplication, optimal-ate loop parameters, final
exponentiation addition chains — is *right*.

The problems sit below the mathematics:

1. **Representation.** `Fp` is a single heap-allocated `mpz_t`. An `Fp12`
   temporary is twelve separate heap allocations. The Miller loop allocates and
   frees roughly 4000 times per pairing. This, not multiplication throughput, is
   where the time goes.
2. **Algorithmic choices.** Affine coordinates with a full modular inversion per
   point double/add; `mpz_mod` (generic division) after every multiply; no
   Montgomery form; no cyclotomic squaring in the final exponentiation.
3. **Twenty confirmed defects** — memory errors, aliasing bugs, and two
   subroutines whose mathematics is simply wrong (Section 4).
4. **Repository state.** 1280 of 1624 tracked files are generated Doxygen HTML.
   Generated autotools output, a release tarball, built dylibs, dependency
   stamps and fontconfig caches are all committed.

Measured baseline (Apple Silicon arm64, `clang -O2`, GMP 6.3.0):

| Operation | BN-462 | BLS12-461 |
|---|---|---|
| Miller loop (opt-ate) | 3.85 ms | 2.49 ms |
| Final exponentiation | 5.07 ms | 5.92 ms |
| **Full pairing** | **~8.9 ms** | **~8.4 ms** |
| G1 scalar mult (plain / 2-split) | — | 1.42 / 0.82 ms |
| G2 scalar mult (plain / 2-split / 4-split) | — | 2.97 / 1.70 / 1.03 ms |
| G3 exponentiation (plain / 2-split / 4-split) | — | 7.68 / 5.99 / 3.49 ms |

A well-tuned BLS12-381 pairing runs on the order of one millisecond on
comparable hardware. The gap is roughly an order of magnitude and it is
addressable without touching the mathematics.

**Target:** 8–15x on the pairing, obtained mostly in Phases 3 and 4. Assembly
(Phase 6) is expected to contribute the last ~1.3x and is deliberately
conditional on measurement.

---

## 1. Decisions (settled)

| # | Question | Decision |
|---|---|---|
| 1.1 | Comparison target | **Add BLS12-381** alongside BN-462 and BLS12-461 |
| 1.2 | Assembly architecture | **Both** AArch64 and x86-64 |
| 1.3 | Meaning of "fast MPZ" | Fixed-width `mpn` limb arrays + Montgomery (see below) |
| 1.4 | Build system | **CMake only** — Bazel dropped |
| 1.5 | KSS16 / Fp4 / Fp8 / Fp16 / EFp6 | **Delete now**, recorded in the backlog (§9) |
| 1.6 | Constant-time | **Required** |

### 1.1 BLS12-381 — better news than expected

Two facts discovered while checking this decision materially simplify the work.

**The tower is identical.** Standard BLS12-381 uses
`Fp2 = Fp[u]/(u²+1)`, `Fp6 = Fp2[v]/(v³−(u+1))`, `Fp12 = Fp6[w]/(w²−v)`.
That is character-for-character the tower this library already implements.
Adding BLS12-381 therefore needs **new parameters only** — no new extension
field code, no new Frobenius derivation, no new sparse multiplication.

I verified the parameters reproduce the published constants exactly:

| Property | Value | Checked |
|---|---|---|
| `X` | `−0xd201000000010000` | — |
| `p` | 381 bits | matches published prime |
| `r` | 255 bits | matches published order |
| `p mod 4` | 3 | `u²=−1` irreducible |
| `p mod 3` | 1 | cube roots of unity exist |
| G1 cofactor | 126 bits | — |

**Three curves need only two limb widths.** At 64-bit limbs: BLS12-381 is 6
limbs, BLS12-461 and BN-462 are both 8. So the assembly burden in Phase 6 is
2 widths × 2 architectures = **four routines**, not six. RELIC reaches the same
conclusion — its `x64-asm-6l` backend is the BLS12-381 one and `x64-asm-8l`
covers the 462-bit range.

BLS12-381 also brings a much cheaper Miller loop, which is why the rest of the
world standardized on it:

| Curve | Loop parameter | Bits | Non-zero digits |
|---|---|---|---|
| BN-462 | `\|6x+2\|` | 117 | 7 |
| BLS12-461 | `\|X\|` | 77 | 3 |
| BLS12-381 | `\|X\|` | 64 | 6 |

Finally, BLS12-381 has published third-party test vectors and a standardized
ciphersuite, which turns Phase 0's oracle from "trust our own Sage script" into
"cross-checked against an external authority."

### 1.2 Both architectures

AArch64 uses `MUL`/`UMULH` with `ADCS` carry chains. x86-64 uses BMI2/ADX
(`MULX`, `ADCX`, `ADOX`), which is where RELIC's advantage comes from. Both get
written, both ship with the C fallback they replace, and CI runs the fallback
and the assembly against the same vectors on every build.

### 1.3 "Fast MPZ" — recorded interpretation

You did not comment on this one, and every other answer presupposes it, so it
stands as plan of record. Say so if this is wrong, because it defines Phase 3.

There is no library called "fast MPZ". Read as *make the multiprecision layer
fast*, the wins in descending order of payoff:

1. Stop using `mpz_t` for field elements. Fixed-size limb array inside the
   struct, on the stack, no allocation.
2. Drop from `mpz_*` to GMP's `mpn_*` layer — caller-supplied limb arrays, no
   allocation, no size negotiation, no sign handling.
3. Montgomery representation, so reduction becomes shifts and multiply-adds
   instead of `mpz_mod`'s division.

### 1.4 CMake only

Bazel is dropped. This removes the `rules_foreign_cc` plumbing for GMP entirely
and makes the library consumable by anything that speaks `find_package` or
`pkg-config`. RELIC is also CMake, so the differential-test harness in Phase 0
can build both from one superbuild.

Sanitizers move from `.bazelrc` configs to `CMAKE_BUILD_TYPE` variants
(`Asan`, `Ubsan`, `Msan`), driven by CTest.

### 1.5 New decision that fell out of adding a third curve

Three curves with two different limb counts forces a choice the two-curve
version did not:

| Option | Behaviour | Cost |
|---|---|---|
| **A** | Compile-time curve selection, one curve per build (RELIC, blst) | fastest; separate binaries |
| B | Runtime context, arrays sized to the maximum | one binary; loops bounded by 8 limbs even for 381-bit |
| C | Compile every curve with specialized code, dispatch at runtime | best of both; largest binary and build |

**Plan of record: A.** It is both the fastest and the least code — the limb
count becomes a compile-time constant, loops unroll, and each assembly routine
is written against a known width. With assembly on two architectures, B would
force every 6-limb operation through 8-limb code paths and forfeit most of the
Phase 6 win. It also makes benchmarking against RELIC apples-to-apples, since
RELIC builds this way too.

Consequence: `-DELIPS_CURVE=BLS12_381|BLS12_461|BN_462` at configure time,
and the bench harness builds three times. This is reversible toward C later if
one binary supporting all three ever becomes a requirement.

### 1.6 Constant-time is required — and it changes earlier phases

This cannot be bolted on in Phase 5. It changes design decisions in Phases 3
and 4, so those phases now carry constant-time requirements directly. Appendix B
maps exactly where it bites; the summary is that **the pairing itself is nearly
constant-time already** and the real work is in scalar multiplication.

---

## 2. What the code actually is

### 2.1 Shape

12,523 lines across 47 `.c` files and 50 headers. Two curve families plus an
unfinished KSS16. A `src/` + `include/ELiPS_bn_bls/` split, one header per
source file, autotools build, one hand-driven `test/main.c`.

### 2.2 The tower

```
Fp                       p = 462 / 461 / 381 bits
Fp2  = Fp[u]/(u² + 1)                       u² = −1
Fp6  = Fp2[v]/(v³ − (1+u))                  cubic non-residue 1+u
Fp12 = Fp6[w]/(w² − v)
```

Sound, and — as §1.1 established — the same tower BLS12-381 needs. `p ≡ 3 (mod 4)`
and `p ≡ 1 (mod 3)` hold for all three curves; I verified each.

Struct layout is nested by value: `Fp12` holds two `Fp6`, each three `Fp2`, each
two `Fp`, each one `mpz_t`. So an `Fp12` is contiguous **except** that each of
its twelve `mpz_t` points at a separately malloc'd limb array. That is the core
performance problem, and it is also why the fix is mechanical: the nesting is
already correct, only the leaf changes.

### 2.3 Verified-correct mathematics

Checked against the defining equations and confirmed. Listing these matters as
much as listing the bugs, because they constrain what the rewrite may not
disturb.

- **BN-462.** `x = 2¹¹⁴ + 2¹⁰¹ − 2¹⁴ − 1`; `p = 36x⁴+36x³+24x²+6x+1` (462 bits),
  `r = 36x⁴+36x³+18x²+6x+1` (462 bits), `t = 6x²+1`. Confirmed `p + 1 − t = r`
  exactly, so the G1 cofactor is 1.
- **BLS12-461.** `X = −2⁷⁷ + 2⁵⁰ + 2³³`; `p = ((X−1)²(X⁴−X²+1))/3 + X`
  (461 bits), `r = X⁴−X²+1` (308 bits). Division by 3 is exact; `r | p+1−t`;
  cofactor 153 bits.
- **BLS12-381.** Verified against published constants — see §1.1.
- **`X_binary_opt` encodes exactly `6x+2`.** Summing the signed digits
  `{116:+1, 115:+1, 103:+1, 102:+1, 16:−1, 15:−1, 2:−1}` gives `6x+2` on the
  nose. Correct optimal-ate loop parameter for BN.
- **`weil()` computes `#E(F_p¹²) = p¹² − t₁₂ + 1` correctly** via the
  `t₂ → t₆ → t₁₂` recurrence. Only the *variable name* is wrong: it lands in
  `EFpd_total`, which the header documents as the twisted curve's order.

### 2.4 Two things that look like bugs and are not

On the record, because a careless refactor will "fix" them and break the library.
I chased both; both are correct.

**The Miller loop bounds are right.** BN iterates `i = bn_X_length+1 = 115` down
to 0 while `6x+2` has its top digit at index **116**. That reads as off-by-one.
It is not: `T` is initialized to `Q`, which *is* the leading `+1` digit — the
standard "seed the accumulator with the most significant digit, then loop the
rest" pattern. Same for BLS12: it runs from `bls12_X_length−1 = 76` while `X`'s
top digit sits at 77 with value `−1`, and `T` is initialized to `−Q`.

**`Fp2_set_ui` setting *both* coefficients is deliberate.** `Fp2_set_ui(&a,1)`
yields `1 + u`, not the field element 1 — which reads as a bug until you find
`set_basis()`, which relies on exactly that to construct the non-residue `1+u`.
A bad name for an intentional broadcast. It stays a hazard (§4.4), not a defect.

---

## 3. Where the time goes

Ranked by expected payoff against the measured 8.9 ms BN pairing.

**1. Heap allocation, ~40–50%.** `Fp12_init` is twelve `mpz_init` calls.
`ff_ltt` — called once per Miller iteration, 116 times per pairing — builds two
`Fp12` and five `Fp2` temporaries, so roughly 34 mallocs and 34 frees per
iteration, ~4000 per pairing, before any arithmetic happens.

**2. Modular inversion in the inner loop, ~20–25%.** Affine coordinates.
`EFp2_ECD` and `ff_ltt` each perform a full `Fp2_inv`. RELIC uses
projective/Jacobian coordinates and pays exactly one inversion per pairing, at
the end. This is the biggest *algorithmic* difference between the two libraries.

**3. `mpz_mod` after every multiply, ~15–20%.** `Fp_mul` is `mpz_mul` then
`mpz_mod`, and `mpz_mod` is truncating division. Montgomery reduction replaces
it with multiply-and-shift — roughly 2x on the reduction step alone — and it
composes with lazy reduction (reduce once per `Fp2`/`Fp6` operation rather than
per `Fp` operation).

**4. No cyclotomic squaring, ~10% of the final exponentiation.** After the easy
part `f^(p⁶−1)(p²+1)` the element lies in the cyclotomic subgroup, where
Granger–Scott compressed squaring is materially cheaper than a generic
`Fp12_sqr`. The final exponentiation is 5–6 ms — the larger half of the pairing
— and it is currently all generic squarings.

**5. `sprintf`/`strtol` inside scalar multiplication.** The split-scalar routines
build their window index one bit at a time with
`sprintf(str,"%c%c%c%c",...)` followed by `strtol(str,&e,2)` — a formatted
string round-trip per bit, ~154 times per scalar multiplication. It is setup
rather than per-iteration cost, and the `char str[5]` buffer does fit exactly,
so it is not a defect. It is merely one of the most expensive possible ways to
assemble four bits, and it is variable-time (§1.6).

**6. Duplication.** `bn_p8sparse.c` and `bls12_p8sparse.c` are byte-identical
apart from function names. Same for `line_ate`, `frobenius`, `skew_frobenius`
and `twist` — the Frobenius pair differs in 21 lines, all names and comments.
Roughly 40% of `src/` is copy-paste. Not a performance issue, but it doubles the
cost of every subsequent fix and doubles the surface for divergence bugs.

---

## 4. Correctness audit

The section to be hardest on. Each item was confirmed by reading the code and,
where the claim is mathematical, by independent computation. Severity:
**C**ritical (wrong results or UB on a live path), **H**igh (UB or leak on a
live path), **M**edium (latent).

### 4.1 Memory and undefined behaviour

| # | Sev | Location | Defect |
|---|---|---|---|
| M1 | **H** | `curve_settings.c:382,394,403` | `bls12_generate_prime()` calls `mpz_init` where it must call `mpz_clear`, on all three exit paths. Re-initializing a live `mpz_t` orphans its limb array. Four leaked integers per call. |
| M2 | **H** | `curve_settings.c:282,468` | `curve_parameters.EFpd_total` is **never `mpz_init`'d** — not in `init_bn_parameters()`, not in `init_bls12_parameters()` — yet `weil()` writes it and both G2 generators read it. It works only because the struct is a zero-initialized global and GMP tolerates `{alloc=0, size=0, d=NULL}`. Formally undefined. |
| M3 | **H** | `curve_settings.c:72,94` vs `bn_clears.c:54` | `curve_parameters.curve_a` is left uninitialized for BN and BLS12 (both `mpz_init` calls commented out), then read by `mpz_cmp_ui` at `curve_settings.c:315` and passed to `mpz_clear` at teardown. `mpz_clear` on an uninitialized `mpz_t`. |
| M4 | **H** | `bn_line_ate.c:90`, `bls12_line_ate.c:90` | `f_ltq()` initializes `Fp2 E` and never clears it. Inside the Miller loop: leaks once per non-zero digit per pairing. |
| M5 | **M** | `bls12_finalexp.c:70,71` | `t4` and `t5` initialized, never cleared; teardown covers only `t0`–`t3`. |
| M6 | **M** | `bls12_scm.c:333` | Cleanup calls `EFp2_init(&twisted_Q_3x)` where it means `EFp2_clear`. Leak plus double-init. |
| M7 | **H** | 14 sites | `int length = mpz_sizeinbase(s,2); char binary[length]; mpz_get_str(binary,2,s);` — GMP writes `length` digits **plus a NUL terminator**, so every one of these overflows its stack buffer by one byte. Sites: `bn_fp.c:206`, `bn_fp2.c:416`, `bn_fp6.c:347`, `bn_fp12.c:275`, `bn_efp.c:223`, `bn_efp2.c:199`, `bn_efp6.c:195`, `bn_efp12.c:220`, `bn_miller_ate.c:47`, `bn_miller_tate.c:39`, `bls12_miller_ate.c:65`, `bls12_miller_tate.c:40`, `fp4.c:168`, `fp8.c:175`. Violates CERT **ARR30-C** and **STR31-C**; the VLA also violates **ARR32-C**, its size deriving from a runtime value. |

M7 is the one worth fixing today regardless of the rest of the plan.

### 4.2 Wrong mathematics

| # | Sev | Location | Defect |
|---|---|---|---|
| X1 | **C** | `bn_fp.c:66` | `Fp_mul_basis_KSS16` reduces modulo `curve_parameters.curve_a` instead of `curve_parameters.prime`. `curve_a` is uninitialized for BN/BLS12 and set to **1** for KSS16, so the result is identically zero. Every `Fp4`/`Fp8` operation built on it is wrong. |
| X2 | **C** | `fp4.c:213`, `fp8.c:221` | `Fp4_sqrt` and `Fp8_sqrt` both compute the group order as `p¹²−1`. `Fp4` has order `p⁴`, `Fp8` has order `p⁸`. Tonelli–Shanks over the wrong group does not return a square root. |
| X3 | **C** | `fp4.c:244`, `fp8.c:252` | The same two functions compute a Tonelli–Shanks step as `mpz_powm(t, 2, r−m−1, p)` — reducing an **exponent** modulo `p`. It must be the integer `2^(r−m−1)`. Two independent errors in one routine. |
| X4 | **H** | `bn_fp6.c:315,316,319,320` | `Fp6_sqrt` has its Frobenius-constant multiplications commented out, so `tmp1` and `tmp2` are built from `x0` alone and the rest is garbage. Returns a value that is not a square root. |
| X5 | **M** | `bn_fp2.c:163` vs `:218` | `Fp2_sqr_KSS16` calls `Fp_mul_basis` (the BN basis) while the otherwise-identical `Fp2_sqr_kss16` calls `Fp_mul_basis_KSS16`. Two functions differing only in case and one call; at most one can be right. |

All five are confined to the KSS16 / `Fp4` / `Fp8` path, which no pairing uses.
That is why the tests pass. Under decision 1.5 this code is deleted, which
retires X1–X5 wholesale — see the backlog entry in §9 for what a future
implementation must not inherit.

### 4.3 API contract violations

| # | Sev | Location | Defect |
|---|---|---|---|
| A1 | **C** | `bn_fp2.c:114,115` | `Fp2_mul_Fp(ANS, A, B)` computes `ANS = ANS * B`, ignoring `A` entirely. Correct only when the caller happens to alias `ANS == A`. |
| A2 | **H** | `bn_fp12.c:143,153,158` | `Fp12_add_mpz`, `Fp12_sub_ui` and `Fp12_sub_mpz` all read `ANS` where they should read `A`. Same defect class as A1. |
| A3 | **H** | `bn_final_exp.c:41,44,70,74`; `bls12_finalexp.c:77,81` | The final exponentiation **writes through its input pointer** (`Fp12_mul(A, &t0, &t1)`). `bn12_opt_ate` calls it as `bn_final_exp_optimal(ANS, ANS)` so the damage is invisible today, but the signature promises `A` is an input. Calling it twice on one value gives two different answers. |
| A4 | **H** | `bn_final_exp.c:87–89`; `bls12_finalexp.c:88–102` | Both final exponentiations **mutate global curve parameters mid-computation** and restore them afterwards — `X_binary[0]` for BN, `bls12_X_length` plus six entries of `bls12_X_binary[]` for BLS12. The library is therefore not reentrant and not thread-safe, and any early return leaves the curve parameters corrupted process-wide. |
| A5 | **M** | `bn_fp.c:99,100` | `Fp_inv` discards `mpz_invert`'s return value. On a non-invertible input the output is undefined and the caller is never told. CERT **EXP12-C**. |
| A6 | **M** | `curve_settings.c` | `generate_bn_prime`/`generate_bn_order` return 1 for success; `bls12_generate_order` returns 0 for success. Every caller ignores the value, so `init_bls12_settings()` proceeds cheerfully after prime generation has failed. CERT **ERR33-C**. |
| A7 | **L** | `curve_settings.c:488` | `gmp_printf("E:y^2=x^3+4\n", curve_parameters.curve_b)` — an argument with no conversion specifier. CERT **FIO47-C**. |
| A8 | **L** | `curve_settings.h`, `bls12_timeprint.h` | `bls12_print_parameters` declared in two headers. |

A3 and A4 together are why the final exponentiation is the least trustworthy
code in the library despite producing correct answers today.

### 4.4 Design hazards

- **`Fp2_set_ui` and `Fp2_cmp_one` disagree.** `set_ui(&a,1)` produces `1+u`;
  `cmp_one(&a)` tests for `1+0u`. So `set_ui(&a,1); cmp_one(&a)` reports "not
  one". Individually defensible, together a trap. Rename to `Fp2_set_all_ui`
  and add a real `Fp2_set_one`.
- **`Fp_cmp` returns 0 or 1, never negative**, despite a name implying
  `memcmp`-style ordering. `fp4.c` then writes
  `while (Fp4_cmp_mpz(&b,set_1)==1)` meaning "while not equal", which reads as
  its own opposite.
- **Global mutable state as algorithm input.** `curve_parameters`, `X_binary`,
  `bls12_X_binary`, `epsilon1/2`, `d12_frobenius_constant[][]`, the RNG `state`
  and the timing globals `t0`/`t1` are all process-wide. Two curves cannot be
  used at once; `init_bn()` and `bls12_inits()` overwrite each other.
- **Benchmarking is welded into the API.** `bls12_plain_G1_scm` calls
  `gettimeofday` and writes the global `bls12_G1SCM_PLAIN`. Measurement belongs
  in a bench harness, not in the function under measurement.
- **RNG.** Every `*_sqrt` and `*_rational_point` builds a *local*
  `gmp_randstate_t` shadowing the global one and seeds it with `time(NULL)` —
  one-second granularity, entirely predictable. Fine for tests, unacceptable for
  key material. CERT **MSC32-C**. Given decision 1.6 this becomes a hard
  requirement, not a nicety.

### 4.5 What the current tests can and cannot catch

This matters more than any single defect above.

Every test has the same shape: draw random `s₁`, `s₂`, compute `e(P,Q)^(s₁s₂)`,
`e([s₁]P,[s₂]Q)` and `e([s₂]P,[s₁]Q)`, check all three agree and the result is
neither 0 nor 1.

That verifies **bilinearity and non-degeneracy**. It does not verify the map is
*the pairing it claims to be*. A truncated Miller loop, a wrong Frobenius
constant, or a botched final exponent can each yield a map that is still
bilinear, still non-degenerate, and still passes all three checks. The suite
would report success.

There are also no known-answer vectors, no `[r]P = O` subgroup checks, no
`assert`s and no exit code — `main` returns 0 unconditionally, so nothing is
CI-gradeable. Fixing this is Phase 0, and it gates everything else.

---

## 5. Target architecture

```
cmake/                    toolchain files, FindGMP, sanitizer presets
include/elips/
  config.h.in             generated: limb width, curve selection
  fp.h  fp2.h  fp6.h  fp12.h
  ep.h  ep2.h             elliptic curves over Fp and Fp2
  pairing.h               public API
  params.h                curve parameter descriptors
src/
  arith/
    fp_mont.c             Montgomery mul/sqr/reduce over mpn
    fp_inv.c              constant-time inversion
    asm/fp_mul_a64_6l.S   Phase 6 — AArch64, 6 limbs (BLS12-381)
    asm/fp_mul_a64_8l.S   Phase 6 — AArch64, 8 limbs (BN-462, BLS12-461)
    asm/fp_mul_x64_6l.S   Phase 6 — x86-64 BMI2/ADX, 6 limbs
    asm/fp_mul_x64_8l.S   Phase 6 — x86-64 BMI2/ADX, 8 limbs
  tower/
    fp2.c fp6.c fp12.c    lazy reduction, cyclotomic operations
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
  ct/                     dudect timing-leakage tests
  vs_relic/               differential tests against RELIC
bench/                    all timing code lives here and nowhere else
```

Six principles, each answering a specific defect above:

1. **One implementation per algorithm.** The BN and BLS12 Miller loops differ
   only in the loop parameter and a final Frobenius correction. They become one
   function taking a `curve_params *`. This deletes ~40% of `src/`.
2. **No global state.** Curve parameters travel in a context struct passed
   explicitly. Fixes §4.4 and A4.
3. **Inputs are `const`.** Makes A1, A2 and A3 compile errors.
4. **Stack allocation.** `typedef limb_t fp_t[FP_LIMBS]` — no malloc on any hot
   path. Fixes performance item 1 and, incidentally, M4, M5 and M6, because
   there is nothing left to leak.
5. **Every function returns a status, or provably cannot fail.** Fixes A5, A6.
6. **Secret-independent control flow by default.** Variable-time paths exist
   only behind an explicit `_vartime` suffix (§1.6, Appendix B).

---

## 6. Phased plan

Each phase has a mechanical exit gate. If the gate is not green, the phase is
not done, and no phase starts before its predecessor's gate passes.

### Phase 0 — Build the oracle. Change nothing else.

The most important phase and the easiest to skip. Today there is no way to tell
a correct refactor from an incorrect one.

1. Tag the current commit as the reference implementation.
2. Write `test/kat/gen_vectors.sage` producing, per curve: field arithmetic
   vectors (mul/sqr/inv/Frobenius at every tower level), curve vectors
   (double/add/scalar-mult with known scalars), and full pairing vectors —
   `e(P,Q)` for fixed generators, printed as twelve `Fp` coordinates. Sage is at
   `/usr/local/bin/sage` and computes these independently, so it is a genuine
   oracle rather than a restatement of our own code.
3. **Cross-check the BLS12-381 vectors against published third-party values**
   before trusting the generator. Then use the *validated generator* for BN-462
   and BLS12-461, where no third-party vectors exist. This is the concrete
   payoff of decision 1.1.
4. Add property tests: `[r]P = O` for G1 and G2; `e(P,Q)^r = 1`;
   `e(P₁+P₂,Q) = e(P₁,Q)·e(P₂,Q)`; field axioms at every tower level.
5. Make the runner exit non-zero on failure.
6. Record baseline timings in `bench/baseline.json`.

**Gate:** the *unmodified* library passes every vector for BN-462 and BLS12-461,
and the runner returns non-zero when a vector is deliberately corrupted.

If the current library fails a Phase 0 vector, that is a finding and it changes
the plan. I put this at roughly 20%, most plausibly in the final exponentiation
given A3 and A4.

### Phase 1 — Repository hygiene and CMake

1. `git rm -r --cached` the 1280 generated Doxygen files, autotools output
   (`configure`, `Makefile`, `config.status`, `libtool`, `aclocal.m4`,
   `src/.deps/`), `src/.libs/`, `elipsbnbls12-1.0.1.tar.gz`, `fontconfig/`,
   `test/a.out`. History is preserved; the working tree stops carrying them.
   (The `.gitignore` for this already landed with rev. 1 of this plan.)
2. `CMakeLists.txt` with `FindGMP`, an installable target, and a generated
   `config.h` carrying `FP_LIMBS` and the selected curve.
3. `-DELIPS_CURVE=BLS12_381|BLS12_461|BN_462`, defaulting to BLS12-381.
4. CTest targets per test group; `CMAKE_BUILD_TYPE` variants `Asan`, `Ubsan`,
   `Msan`, `Release`.
5. A CMake superbuild that also fetches and builds RELIC, for the Phase 4
   differential and benchmark comparison.
6. CI matrix: {AArch64, x86-64} × {3 curves} × {Release, Asan, Ubsan, Msan}.

**Gate:** `cmake --build . && ctest` green from a clean checkout, for all three
curves, on both architectures. Phase 0's vectors still pass byte for byte.

The sanitizer runs will immediately surface M1–M7. That is intended — it is
Phase 1 proving the tooling works before it is needed.

### Phase 2 — Fix the confirmed defects, on the existing architecture

Deliberately *before* the rewrite, so each fix is a small reviewable diff
against a passing oracle rather than being smuggled inside a 10,000-line change.

Order: M7 (stack overflows) → M1, M2, M3 (initialization) → A1, A2 (aliasing) →
M4, M5, M6 (leaks) → A5, A6, A7, A8 → delete the KSS16 / `Fp4` / `Fp8` /
`Fp16` / `EFp6` tree per decision 1.5, retiring X1–X5 and recording §9.

Each fix is TDD: write the failing test, then the fix, then watch it pass. M7
gets a test compiled with ASan that currently reports a stack-buffer-overflow.

A3 and A4 are explicitly **not** fixed here — they are structural, and Phase 3
deletes the globals that cause them.

**Gate:** all Phase 0 vectors pass; ASan, UBSan and MSan clean; timings within
noise of baseline, since this phase is not supposed to change speed.

### Phase 3 — The representation change

The big one, and the source of most of the speedup.

1. `typedef uint64_t limb_t; typedef limb_t fp_t[FP_LIMBS];` with `FP_LIMBS`
   fixed at configure time from the selected curve — 6 or 8. No allocation
   anywhere in `fp_t`.
2. Montgomery form: `to_mont`, `from_mont`, CIOS multiply-and-reduce over
   `mpn_*`. Conversion happens at the API boundary only.
3. Rewrite `fp2`/`fp6`/`fp12` on `fp_t` with **lazy reduction** — accumulate
   double-width intermediates through an `Fp2` multiplication and reduce once at
   the end rather than after each `Fp` operation.
4. Curve parameters move into a context struct; the globals go away, taking A3
   and A4 with them.
5. Miller loops merge into one parameterized implementation.
6. **Constant-time from the start** (decision 1.6): every `fp_t` operation is
   secret-independent. Reduction uses conditional subtraction via masked select,
   never a branch. `mpn_*` multiplication and addition are already
   data-independent in time; `mpz_mod` and `mpz_invert` are not, and both are
   gone by the end of this phase.
7. Add BLS12-381 parameters. Per §1.1 this is a parameter file, not new
   arithmetic — and it is the curve with third-party vectors, so it is the one
   to land and validate *first*, then BLS12-461 and BN-462 behind it.

Do this **one tower level at a time**, keeping the old implementation beside the
new and differential-testing them on random inputs until the level is proven,
then deleting the old.

**Gate:** every Phase 0 vector passes for all three curves; the differential
test finds no disagreement across 10⁶ random inputs per level; sanitizers clean;
measured speedup ≥ 3x on the full pairing. If it is under 3x, stop and
re-measure — the model of where time goes is wrong and the rest of the plan
needs revisiting before more work goes in.

### Phase 4 — Algorithmic improvements

1. **Projective coordinates** for `EFp` and `EFp2`, so the Miller loop performs
   zero inversions instead of one per iteration, and one at the end.
2. **Cyclotomic squaring** (Granger–Scott) in the final exponentiation and in
   `G3` exponentiation, with a comment stating the precondition — the element
   must lie in the cyclotomic subgroup, which holds after the easy part and
   nowhere else. Getting this wrong produces silently wrong results, so the
   precondition is asserted in debug builds.
3. **Compressed squaring** for the `f^x` chains.
4. **Sparse multiplication** without the full `Fp12` temporary that
   `Pseudo_8_sparse_mul` currently allocates.
5. **Constant-time fixed-window scalar multiplication**, replacing both the
   binary ladder and the `sprintf`/`strtol` window construction of §3.5. Table
   lookup is a full linear scan with masked select, not an index. This subsumes
   what would otherwise be Phase 5 work and is placed here because it replaces
   the same code.
6. **GLV / split-scalar decomposition** rewritten constant-time — the existing
   2-split and 4-split routines branch and use `mpz_tdiv_qr` on secret scalars.

Each item lands and is measured independently. Anything not paying for itself
is reverted rather than kept out of sentiment.

**Gate:** vectors pass; cumulative speedup ≥ 8x on the full pairing; every
individual optimization shown to be a measured win; head-to-head BLS12-381
benchmark against RELIC recorded, using the Phase 1 superbuild.

### Phase 5 — Remaining constant-time work and API hardening

Phases 3 and 4 carry most of the constant-time burden. What remains:

1. **Constant-time field inversion.** `x^(p−2)` via a fixed addition chain, or
   safegcd/Bernstein–Yang if measurement justifies the complexity. Note from
   Appendix B that inversion is needed only once per pairing after Phase 4, so
   the simple addition chain is likely sufficient — do not reach for safegcd
   without a measurement saying so.
2. **Subgroup membership checks** on every point entering the API. Currently
   absent and required against small-subgroup attacks.
3. **Full input validation** at the trust boundary: points on curve, points in
   subgroup, scalars in range.
4. **A real CSPRNG.** Replace every `time(NULL)`-seeded local `gmp_randstate_t`
   with `getrandom`/`arc4random_buf`, and remove the shadowing locals (§4.4).
5. **`dudect` timing-leakage tests in CI**, on both architectures.
6. Serialization and deserialization in the standard compressed point formats,
   which BLS12-381 interoperability requires.

**Gate:** `dudect` reports no leakage on any secret-dependent path on either
architecture; the variable-time fast paths are reachable only through
explicitly `_vartime`-suffixed functions; malformed and off-curve inputs are
rejected rather than processed.

### Phase 6 — Assembly

**Only if Phase 4's profile shows limb arithmetic is the bottleneck.** Do not
write assembly to make `malloc` faster. Expectation is that after Phases 3–4 the
remaining hot spot genuinely is Montgomery multiplication, in which case:

Four routines, per §1.1:

| Width | AArch64 | x86-64 |
|---|---|---|
| 6 limbs (BLS12-381) | `MUL`/`UMULH` + `ADCS` | `MULX`/`ADCX`/`ADOX` |
| 8 limbs (BN-462, BLS12-461) | `MUL`/`UMULH` + `ADCS` | `MULX`/`ADCX`/`ADOX` |

1. Start with 6-limb AArch64 — smallest routine, on the development machine, for
   the curve with external vectors. Prove the approach there before scaling.
2. Every assembly routine ships with the C fallback it replaces. Both run
   against the same vectors on every CI build, on every architecture.
3. Runtime CPU feature detection for BMI2/ADX with a portable default.
4. Assembly is constant-time by construction — no data-dependent branches, no
   data-dependent memory addressing. This is verified, not assumed.

**Gate:** assembly and C paths produce bit-identical output on all vectors, for
all four combinations; measured win ≥ 1.2x, otherwise that routine is deleted
rather than carried; `dudect` still clean.

---

## 7. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Phase 0 finds the current library is already wrong somewhere | ~20% | Better to know now. Fix before refactoring, not during. |
| Montgomery conversion introduces a subtle sign or carry bug | medium | Differential-test against the mpz implementation, 10⁶ random inputs per level, before deleting the old path. |
| Lazy reduction overflows the accumulator | medium | Prove the bound explicitly for each of the three moduli and assert it in debug builds. This is the classic lazy-reduction failure mode. |
| Constant-time requirement slows the library more than expected | medium | Appendix B shows the pairing is nearly constant-time already; the cost concentrates in scalar multiplication, where windowing recovers most of it. Keep `_vartime` variants for public inputs. |
| Three curves × two architectures multiplies CI and review load | medium | Compile-time curve selection (§1.5 option A) keeps each build simple; the matrix is CI's problem, not the code's. |
| Assembly wins nothing after Phases 3–4 | medium | Phase 6 is explicitly conditional on measurement, and each of the four routines is independently revertible. |
| Scope grows without end | **high** | This document is the scope. Anything not in it is a separate proposal, and §9 is where it waits. |

---

## 8. Explicitly not doing

- Curve families beyond BN-462, BLS12-461 and BLS12-381.
- A protocol layer (BLS signatures, IBE). Different project.
- GPU or multi-threaded pairings.
- One binary supporting all three curves at runtime — see §1.5, option A is the
  plan of record and option C remains available later.
- Keeping API compatibility with the current headers. The current API cannot
  express `const` inputs or a curve context, and preserving it would forfeit the
  fixes to A3 and A4.

---

## 9. Backlog

Work deliberately deferred, recorded so it is not lost. Nothing here is
scheduled; each item is a future proposal.

### 9.1 KSS16 (deleted in Phase 2, per decision 1.5)

The KSS16 support — `fp4.c`, `fp8.c`, the `Fp16`/`EFp4`/`EFp8`/`EFp16` type
declarations, and the `*_KSS16` / `*_kss16` variants in `bn_fp.c` and
`bn_fp2.c` — is removed in Phase 2 and can be recovered from git history at the
Phase 0 reference tag.

**Do not restore it from history.** A future implementation must be written
against a reference paper, because the deleted code carried five confirmed
mathematical defects (§4.2):

- `Fp_mul_basis_KSS16` reduced modulo the curve coefficient instead of the prime,
  making the result identically zero.
- `Fp4_sqrt` and `Fp8_sqrt` both used `p¹²−1` as the group order, when the
  correct orders are `p⁴−1` and `p⁸−1`.
- The same two functions reduced a Tonelli–Shanks *exponent* modulo `p`.
- `Fp2_sqr_KSS16` and `Fp2_sqr_kss16` disagreed on which basis multiplication
  to use; at most one was right.

What a future KSS16 effort needs, beyond fixing the above:

- A `k=16` tower: `Fp2 → Fp4 → Fp8 → Fp16`, with a correctly chosen non-residue
  and a derivation of the Frobenius constants, none of which existed.
- No `fp16.c` was ever written; `Fp16`, `EFp4`, `EFp8` and `EFp16` were type
  declarations with no implementation.
- The optimal-ate loop parameter, line functions and a `k=16` sparse
  multiplication.
- A `k=16` final exponentiation addition chain.
- KSS16 parameter generation, which existed but was never validated — the
  `generate_kss16_parameters` routine calls `exit(0)` on a non-prime result.

Realistically this is a project of comparable size to Phases 3 and 4 combined,
and it should not begin until the BN/BLS12 path is finished and proven.

### 9.2 Other deferred items

- **`Fp6_sqrt`** (X4) — deleted with the KSS16 tree since nothing else calls it.
  If a square root in `Fp6` is ever needed, write it fresh.
- **Tate, plain-ate and x-ate pairings** — kept through Phase 2 as cross-checks
  against opt-ate, and worth keeping in the test suite afterwards even though
  they are superseded for production use.
- **Option C from §1.5** — runtime dispatch across all three curves in one
  binary, if a consumer ever needs it.
- **Multi-pairing / products of pairings** — `e(P₁,Q₁)·e(P₂,Q₂)···` sharing one
  final exponentiation. A large win for verification-heavy protocols, and a
  natural follow-on once Phase 4 lands.

---

## Appendix A — Verification log

Claims checked by computation rather than by reading.

| Claim | Method | Result |
|---|---|---|
| BN-462 `p`, `r`, `t` polynomials, and `p+1−t = r` | Python, exact integers | confirmed, cofactor 1 |
| BN-462 `p` is 462 bits, `r` is 462 bits | Python | confirmed |
| `X_binary_opt` = `6x+2` | summed signed digits | exact match |
| BN Miller loop bound is correct despite appearances | traced MSB-preload against digit index 116 | **not a bug** |
| BLS12-461 `p` is 461 bits, `r` is 308 bits, `3 \| (X−1)²(X⁴−X²+1)`, `r \| p+1−t` | Python | confirmed, cofactor 153 bits |
| BLS12-461 Miller loop bound is correct | traced `T = −Q` preload against digit index 77 | **not a bug** |
| BLS12-381 `p` and `r` match published constants | Python, compared to standard hex values | exact match |
| BLS12-381 tower matches the library's existing tower | compared defining polynomials | identical; parameters only |
| `p ≡ 3 (mod 4)` and `p ≡ 1 (mod 3)` for all three curves | Python | confirmed; tower valid for each |
| Three curves need only two limb widths | `ceil(bits/64)` for 381, 461, 462 | 6, 8, 8 |
| Miller loop lengths and digit weights across curves | NAF weight computation | BN-462: 117 bits; BLS12-461: 77; BLS12-381: 64 |
| `mpz_get_str` needs `sizeinbase+1` bytes for positive input | GMP documentation plus arithmetic | overflow of exactly 1 byte at 14 sites |
| `Fp2_set_ui` dual-set is relied on by `set_basis()` | read `set_basis`, checked `Fp2_mul_basis` implements `×(1+u)` | **not a bug** |
| `char str[5]` in the split-scalar routines does not overflow | 4 chars + NUL = 5 | **not a bug**, but variable-time |
| Baseline timings | built with `clang -O2`, ran both suites | recorded in §0 |
| RELIC ships no BN-462 or BLS12-461 preset | enumerated `relic/preset/` and `relic_ep.h` | confirmed; BLS12-381 is present, which decision 1.1 exploits |

Three hypotheses formed while reading were **disproved** by this process — the
Miller loop bounds, `Fp2_set_ui`, and the `char str[5]` buffer. All three would
have caused a regression had they been "fixed". This is the argument for Phase 0
in miniature.

---

## Appendix B — Where constant-time actually bites

Decision 1.6 makes this required, so it is worth being precise about scope
rather than applying blanket paranoia. Blanket paranoia is slow and it hides the
places that genuinely matter.

**The pairing is nearly constant-time already.** The Miller loop's control flow
is driven by the *curve parameter* (`6x+2` for BN, `X` for BLS12), which is
public and fixed at compile time. The loop structure therefore leaks nothing.
The only secrets entering `e(P,Q)` are the point coordinates, and every
operation on them is straight-line. Once the underlying field operations are
constant-time — Phase 3, item 6 — the pairing is done.

**Scalar multiplication is where the work is.**

| Location | Leak | Fix | Phase |
|---|---|---|---|
| `EFp_SCM`, `EFp2_SCM`, `EFp12_SCM` | branch on each scalar bit | fixed-window, masked table scan | 4 |
| 2-split / 4-split SCM | `mpz_tdiv_qr` on secret scalars; `sprintf`/`strtol` window build | constant-time decomposition and window build | 4 |
| `Fp12_pow` (G3 exponentiation) | branch on each exponent bit | fixed-window over the cyclotomic subgroup | 4 |
| `Fp_inv` / `Fp2_inv` | `mpz_invert` is variable-time | addition chain `x^(p−2)`, or safegcd | 5 |
| All `*_sqrt`, `*_rational_point` | `time(NULL)`-seeded RNG | CSPRNG | 5 |
| `Fp_cmp` family | early-return comparison | constant-time compare returning a mask | 3 |

**What stays variable-time on purpose.** Operations on public data — verifying
a signature against a public key, hashing to a curve point from public input,
deserializing a public point — do not need the penalty. These keep fast paths
behind an explicit `_vartime` suffix so the choice is visible at every call site
rather than implied.

**Inversion deserves a note.** After Phase 4 introduces projective coordinates,
field inversion happens roughly *once per pairing* instead of once per Miller
iteration. A constant-time addition chain for `x^(p−2)` costs around 400
squarings and a handful of multiplications for a 381-bit prime — negligible
against a whole pairing. safegcd is several times faster but considerably more
intricate and easier to get subtly wrong. Start with the addition chain and only
reach for safegcd if a measurement says inversion is material.


## Gate revision (2026-09-06)

## Exit gate (revised 2026-09-06)

Vectors pass; every individual optimization shown to be a measured win; and the
**RELIC-relative gate**: BLS12-381 pairing within 25% of RELIC built with the
same class of arithmetic on the same machine, with scalar multiplication within
1.5x.

The original gate asked for a cumulative 8x against the legacy layer. It was
replaced because it was set before anyone knew where the time actually went, and
the Phase 0 model that produced it turned out to be wrong about inversion by a
factor of four or five. A ratio against the legacy library also measures the old
code's weaknesses rather than this one's quality. "Comparable to RELIC" was the
stated goal from the start, so the gate now says that directly.
