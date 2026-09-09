# BLS12-381 pairing research: subfield-scaled prepared Miller replay

**Status:** algebraically justified, tested local prototype, not production code.
**Novelty:** not established. The implementation improvement is real; a new
cryptographic result has not been established.

- Analysis: 2026-09-08 to 2026-09-09.
- Repository: `/Users/lamda/Documents/GitHub/elips_bn_bls`.
- Revision: `5c241fcc590b4381b9fcc35742e993b5907d900f`.
- Production repository files: unchanged.
- Prototype: [experiment.c](experiment.c).
- Timings: [first run](benchmark.log), [full-width test run](benchmark-fullwidth.log).

## 1. Result

The current prepared pairing multiplies a dense Fp12 accumulator by each
three-coefficient line using **15 Fp2 multiplications**. A normalized,
subfield-scaled kernel uses **8 Fp2 multiplications** instead.

The distinguishing implementation choice is to compute **twice the product**.
That avoids interpolation divisions, and final exponentiation removes the
extra factor. This deliberately optimizes the Miller value *up to an Fp2
factor*, not its exact representation.

Measured on an Apple M2 Max:

| Prepared terms | Current Miller | Candidate Miller | Less time | Current Miller + final exp | Candidate Miller + final exp | Less time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 330.167 us | 264.700 us | 19.83% | 791.533 us | 725.733 us | 8.31% |
| 2 | 515.033 us | 387.067 us | 24.85% | 975.167 us | 849.000 us | 12.94% |
| 8 | 1631.900 us | 1121.967 us | 31.25% | 2088.333 us | 1582.867 us | 24.20% |

These are **prepared arithmetic paths**, not complete checked public API calls.
They exclude subgroup validation, affine conversion, and offline table
construction. They include the candidate's additional online P normalization.
This is not a claim of these speedups for an arbitrary unprepared pairing.

## 2. What the implementation already does

Relevant source anchors:

- `src/pairing/miller.c:59`: 15M2 sparse line multiplication.
- `src/pairing/miller.c:110`: homogeneous tangent and point doubling.
- `src/pairing/miller.c:169`: homogeneous chord and point addition.
- `src/pairing/miller.c:196`: two Fp2-by-Fp line coefficient scalings.
- `src/pairing/miller.c:407`: prepared, shared Miller loop.
- `src/pairing/miller.c:547`: BLS12 fast final-exponentiation chain.
- `src/arith/fpx.c:216`: general Fp12 squaring, costing 12M2.
- `src/arith/fpx.c:322`: Granger–Scott cyclotomic squaring.
- `src/arith/fpx.c:432`: Karabina compressed squaring.
- `include/elips/fp_params.h:59-74`: 64 Miller iterations and 69 lines.
- `include/elips/pairing.h:41`: exact raw prepared-Miller output contract.

The main opportunities are not “add a Miller loop,” “use sparse arithmetic,”
“use cyclotomic squaring,” or “share final exponentiation.” Those are already
present. The remaining opportunity is to change the algebraic representation
the line-multiplication kernel is required to preserve.

The top-of-file Jacobian discussion in `miller.c` describes an older formula;
the implemented `dbl_line` and `add_line` use homogeneous coordinates.
The derivations below follow the actual code.

## 3. Normalize prepared lines

Use the repository's tower:

    Fp2 = Fp[u]/(u²+1)
    Fp6 = Fp2[v]/(v³-xi),  xi = 1+u
    Fp12 = Fp6[w]/(w²-v).

A prepared line is evaluated as

    L = a*yP + c*w³ + b*xP*w⁵.

Here a, b, c are in Fp2 and depend only on Q. For each fixed Q, store

    beta = c/a,    gamma = b/a.

For each online P compute

    iy = 1/yP,    xy = xP/yP.

The normalized line is

    Lhat = 1 + B*w³ + C*w⁵,
    B = beta*iy,    C = gamma*xy.

Thus L = (a*yP)*Lhat. Both paths still use two Fp2-by-Fp multiplications
per line.

### Costs that must not be hidden

- Online, per P: one Fp inversion and one Fp multiplication.
- Offline, per 69-line Q table: batch-invert 69 nonzero a coefficients and
  multiply each inverse into b and c. This costs **1I2 + 342M2** beyond
  generating the original lines:
  `3*(69-1) + 2*69 = 342`.
- Stored table: 19,872 bytes becomes 13,248 bytes, a one-third reduction.
- Batch-normalizing several online P values can amortize their inversions;
  the measured prototype does not implement this further improvement.

The setup cost matters for one-use Q values. The proposal primarily targets
reused G2 arguments, such as a fixed generator and recurring public keys.

### Why division is permitted

For a nonidentity P of odd prime order r, yP is nonzero. For a nonidentity
order-r Q, every doubling state in this loop is finite and has nonzero Y.
At the addition positions the running point is not +/-Q, so H is nonzero.

These facts were checked using the exact signed scalar prefixes for all 69
lines. They hold for every nonidentity Q in this prime-order subgroup, not
just the tested generator multiples. Inputs still need validation, and the
prototype checks each a coefficient before batch inversion.

This is not permission to normalize arbitrary unvalidated point tables.

## 4. Derive the eight-multiplication kernel

Write the accumulator as f = A + D*w, with A,D in Fp6. Since

    Lhat = 1 + w³*(B+C*v),

we obtain

    f*Lhat = [A + v²*D*(B+C*v)] + [D + v*A*(B+C*v)]*w.

We therefore need two products of a quadratic polynomial by the same
linear polynomial.

For T = t0 + t1*v + t2*v², form four Fp2 products:

    m0   = t0*B
    m3   = t2*C
    mplus  = (t0+t1+t2)*(B+C)
    mminus = (t0-t1+t2)*(B-C).

Polynomial evaluation at 0, 1, -1, and infinity gives:

    H0 = 2*(m0 + xi*m3)
    H1 = mplus - mminus - 2*m3
    H2 = mplus + mminus - 2*m0.

Then

    H(T) = H0 + H1*v + H2*v² = 2*T*(B+C*v).

Unlike ordinary interpolation, these formulas require no halving.
Compute H(A) and H(D), sharing B+C and B-C, then set

    Anew = 2*A + v²*H(D)
    Dnew = 2*D + v*H(A).

The output is **2*f*Lhat**. Its cost is **8M2** plus additions,
subtractions, and cheap multiplications by xi/v. This is not an optimality
proof: a lower-rank or better scheduled kernel has not been ruled out.

## 5. Why the pairing stays exactly the same

Let E = (p¹²-1)/r. In BLS12-381,

    (p²-1) divides E.

Consequently s^E = 1 for every nonzero s in Fp2.

Suppose the candidate accumulator is s times the original at an iteration
boundary. Squaring changes the discrepancy to s². Replacing a line L by
2*L/(a*yP) introduces another nonzero Fp2 factor. Thus the discrepancy stays
in Fp2 throughout the whole loop, including multi-pairing interleaving.

After final exponentiation, both results are identical.

The current fast BLS12 chain computes exponent 3E rather than E; the argument
also holds for 3E. The prototype tested both conventions. This does not
silently change the library's documented cubed-pairing convention.

### API consequence

The candidate **does not preserve the raw Miller value**. The existing
`pairing_miller_prec` promises bit-identical output to `pairing_miller`.

Do not silently replace it. Introduce a distinct normalized prepared-table
type and a separate final-exponentiated path, or explicitly version the
low-level contract. Preserve the existing exact-raw path.

## 6. Field-operation accounting

M2 denotes a general Fp2 multiplication; Mp an Fp multiplication.
Here M2 = 3Mp and S2 = 2Mp in the actual implementation.
Counts below do not treat additions as free in runtime measurements.

| Per evaluated line | Current code | Normalized Karatsuba control | Scaled interpolation candidate |
| --- | ---: | ---: | ---: |
| Dense accumulator times line | 15M2 | 10M2 | 8M2 |
| Two Fp2-by-Fp coefficient scalings | 4Mp | 4Mp | 4Mp |
| Combined multiplication count | 49Mp | 34Mp | 28Mp |

The 10M2 control is also implemented and tested. It separates the benefit of
normalization from the benefit of the eight-product kernel.

For 69 lines, the candidate saves 483M2 = 1,449Mp in line multiplication.
Subtract the extra online 1Ip + 1Mp per P when comparing whole loops.

With the existing 64 general Fp12 squarings included:

    current:   64*36 + 69*49     = 5685 Mp
    candidate: 64*36 + 69*28 + 1 = 4237 Mp + 1 Ip.

These counts exclude all additive work and offline preparation, and apply
to a single prepared Miller loop before final exponentiation.

## 7. Another concrete reduction: tangent construction

The implemented homogeneous tangent has

    a = 2*xi*Y*Z²
    b = -3*X²*Z
    c = 3*X³ - 2*Y²*Z.

Using Y²*Z = X³ + b_twist*Z³, it is Z times

    a' = 2*xi*Y*Z
    b' = -3*X²
    c' = Y² - 3*b_twist*Z².

Discarding the common nonzero Z is legal for the reduced pairing. The
surrounding point-doubling code already computes YZ, X², Y², Z², and
3*b_twist*Z². Reusing these removes **four M2 per doubling** from line
formation. For BLS12-381, b_twist = 4*xi, so computing 12*xi*Z² by
additions and `mul_xi` removes the existing general multiplication by
3*b_twist as well.

The existing combined tangent/doubling costs 11M2 + 3S2, counting its
general constant multiplication. The derived form can use 6M2 + 3S2,
with additional cheap constant operations. Over 64 doublings that is
320M2 removed.

This identity passed 32 projectively rescaled point checks. A complete
optimized doubling implementation and its timing were **not** produced;
do not attribute the prepared-path benchmark gains to this reduction.
This is a conventional curve-equation/denominator-elimination improvement,
not a novelty claim.

## 8. Validation actually performed

All following checks passed:

1. 2,000 full-width Fp12 kernel tests against dense multiplication.
   Cases include zero B/C, B=C, B=-C, zero accumulator, and identity.
   Both the exact 10M2 control and doubled-output 8M2 candidate are checked.
2. 32 homogeneous tangent-line identities on projectively rescaled points.
3. Four groups of random small nonzero generator multiples, eight P/Q pairs
   per group. For n=0 through n=8, both candidate kernels match the current
   complete reduced pairing: 36 comparisons per candidate.
4. An additional comparison using the exact plain final exponent E.
5. The same prototype tests under AddressSanitizer and UBSan, with the
   repository's arithmetic also instrumented.
6. Existing tests `pairing.BLS12_381`, `edge.BLS12_381`, and
   `pairing.detects_corruption`: 3/3 passed in the isolated baseline build.
7. Exact integer interpolation identity on all six bilinear basis pairs.
8. Exact scalar-prefix checks of normalization preconditions.
9. Exact divisibility and BLS12 hard-exponent identity checks.

The initial standalone compile could not locate Homebrew's GMP header;
adding `/opt/homebrew/include` and `/opt/homebrew/lib` fixed that environment
issue. Subsequent builds passed with warnings treated as errors.

Not performed: a complete repository test run, other curves, assembly review,
constant-time leakage measurements, multi-architecture benchmarks, production
integration, or an independent external implementation comparison of the
prototype itself. The existing pairing KAT suite supplies the baseline
independent-oracle connection. Tests are evidence, not a formal verification
of the C implementation.

## 9. Benchmark method and limits

- CPU: Apple M2 Max, arm64, macOS Darwin 25.6.0.
- Compiler: Apple Clang 21.0.0, `-O3`, C11.
- GMP: 6.3.0, Homebrew installation.
- 15 repetitions of 30 calls for every method and workload.
- Warmup pass; methods interleaved with rotating order.
- Two runs yielded similar ratios. Tables in section 1 use the second run.
- Current and candidate kernels were built into the same test executable.
- The eight-product method also beat the 10M2 normalized control.
- No secret or production inputs were used.
- Timings are local measurements, not confidence intervals or portable claims.

The candidate includes one inversion per P inside each timed call. It does
not hide this cost in offline setup. Both sides exclude subgroup validation,
point-to-affine conversion, prepared-table construction, and the public API's
outer chunk-accumulator multiply. Final exponentiation is included only in
the explicitly labeled Miller + final measurements.

## 10. Prior work and novelty assessment

Relevant prior work located during targeted public searches:

1. Costello and Stebila, **Fixed Argument Pairings**,
   [ePrint 2010/342](https://eprint.iacr.org/2010/342).
   Precomputing the G2 work is established and already present here.
2. Mori, Akagi, Nogami, and Shirase, **Pseudo 8–Sparse Multiplication for
   Efficient Ate–Based Pairing on Barreto–Naehrig Curve**,
   [DOI 10.1007/978-3-319-04873-4_11](https://doi.org/10.1007/978-3-319-04873-4_11).
   This is directly relevant prior art for normalized/pseudo-sparse arithmetic.
   “8-sparse” is a sparsity label, not evidence that the paper uses 8M2.
3. Devegili et al., **Multiplication and Squaring on Pairing-Friendly Fields**,
   [ePrint 2006/471](https://eprint.iacr.org/2006/471).
   Extension-field multiplication optimization is an established subject.
4. **Fast AVX-512 Implementation of the Optimal Ate Pairing on BLS12-381**,
   [ePrint 2025/1283](https://eprint.iacr.org/2025/1283).
   Relevant modern implementation baseline, but a different hardware target.
5. **Simpler and Faster Pairings from the Montgomery Ladder**,
   [ePrint 2025/672](https://eprint.iacr.org/2025/672).
   Relevant alternative algorithmic direction; its title alone does not
   establish an applicable speedup for this BLS12-381 implementation.

These searches identified prior work; they were not an exhaustive,
full-text comparative novelty review. Search absence is not proof of novelty.
No claim is made that the precise scaled 8M2 construction is absent from
existing papers or libraries.

The defensible result is:

> A tested arithmetic optimization absent from the current code, assembled
> from established ideas, with a potentially research-worthy combination
> whose publication-level novelty remains unresolved.

The research question worth pursuing is whether optimizing an entire
square-and-line step *modulo a nonzero subfield scalar* admits a better
arithmetic circuit than separately optimizing exact field operations.
The prototype demonstrates a useful instance of this approach. A lower
operation count for the fused step has not been found or proved here.

## 11. Recommendation

1. Use the normalized 8M2 prepared path as the first implementation candidate,
   with an explicit final-exponentiation-only contract.
2. Retain subgroup checks and validate normalized-table construction.
3. Keep the current raw Miller API and its existing equality tests.
4. Add proof-oriented tests for equality after both final-exponent conventions.
5. Benchmark public checked APIs, realistic key reuse, cold/warm tables, and
   more than eight terms before selecting a default.
6. Separately implement and measure the homogeneous tangent reduction.
7. Compare the precise kernel with the pseudo-sparse and extension-field
   literature before labeling any part a novel algorithm.

## 12. Replay

Existing local binaries:

```sh
/tmp/elips-bls381-research.aoHAvY/experiment
/tmp/elips-bls381-research.aoHAvY/experiment bench
/tmp/elips-bls381-research.aoHAvY/experiment-asan
ctest --test-dir /tmp/elips-bls381-research.aoHAvY/build \
  --output-on-failure \
  -R '^(pairing\.BLS12_381|edge\.BLS12_381|pairing\.detects_corruption)$'
```

The prototype includes the original `miller.c` to obtain access to its
static reference kernel, and links against the baseline static library.
It is deliberately a research harness, not a proposed production file.
Paths are temporary and may be cleared by the operating system.
