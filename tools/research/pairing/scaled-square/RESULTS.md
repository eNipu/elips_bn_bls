# Result: scaled Fourier squaring in the prepared Miller loop

**Status:** working research prototype, not production integration.
**Novelty:** not established. This is an implementation improvement assembled
from established interpolation and denominator-elimination ideas.

## Measured result

Relative to the previously committed **8M2 normalized-line path**, keeping
the same line kernel and changing only the general squaring:

| Workload | Run 1: less median time | Run 2: less median time |
| --- | ---: | ---: |
| One general Fp12 squaring | 10.15% | 10.18% |
| Prepared Miller loop, 1 term | 5.71% | 5.37% |
| Prepared Miller + final exponent, 1 term | 1.88% | 1.61% |
| Prepared Miller loop, 2 terms | 3.24% | 3.66% |
| Prepared Miller + final exponent, 2 terms | 2.05% | 1.19% |
| Prepared Miller loop, 8 terms | 0.97% | 0.79% |
| Prepared Miller + final exponent, 8 terms | 1.07% | 0.62% |

These are incremental gains, **not** another comparison against the original
15M2 line kernel. The new one-term Miller medians were 250.267 and 255.600 us,
against 265.433 and 270.100 us for the 8M2 control.

A Chung–Hasan-style 4S4+M4 scaled squaring was included as a stronger control.
Fourier squaring took about 4.4% less time than this control's standalone
square in both runs. Its advantage over that control's complete pairing was
small and noisy; at eight terms it did not consistently win the Miller loop.
This is not a claim to beat every optimized implementation of that formula.

Raw data:

- `benchmark-20260909-run1.log`
- `benchmark-20260909-run2.log`

`analyze_bench.py` recomputes these tables from the `SAMPLE` records and also
reports median within-block ratios. No confidence interval, multi-machine
conclusion, or production API speedup is asserted.

## Construction

Use the existing field tower, just regroup its coordinates:

```text
u^2 = -1, xi = 1+u
s = w^3, s^2 = xi
Fp4 = Fp2[s]/(s^2-xi)
Fp12 = Fp4[w]/(w^3-s)
f = a + b*w + c*w^2
```

There is no expensive basis conversion. If the stored halves are
`(g0,g2,g4)` and `(g1,g3,g5)`, then the Fp4 coordinates are
`a=(g0,g3)`, `b=(g1,g4)`, `c=(g2,g5)`.

Evaluate the quadratic at the four fourth roots of unity and infinity:

```text
P = (a+b+c)^2
M = (a-b+c)^2
I = (a+u*b-c)^2
J = (a-u*b-c)^2
K = 4*c^2

S = P+M
T = I+J
U = P-M
V = -u*(I-J)
```

The output coordinates are

```text
h0 = S+T-K + s*(U-V)
h1 = U+V   + s*K
h2 = S-T

h0 + h1*w + h2*w^2 = 4*f^2.
```

This is length-four Fourier interpolation with the degree-four coefficient
supplied separately. Multiplication by `u` is a coefficient swap and sign
change. Multiplication by `s` swaps Fp2 coordinates and applies `mul_xi`.
No division, inversion, or general constant multiplication is required.

The earlier Gaussian candidate at `{0,1,-1,u,infinity}` is retained in both
the Python and C files. It has the same multiplication count but more
interpolation additions and measured slower.

## Why the pairing is unchanged

For `E=(p^12-1)/r`, `p-1` divides E on these curves. Thus every nonzero Fp
scale vanishes after exponentiation by E, and also after 3E for ELiPS's
fast BLS12 convention.

Relative to the same normalized-line loop, the scalar discrepancy obeys

```text
scale_0 = 1
scale_(k+1) = 4 * scale_k^2 mod p.
```

Each subsequent line multiplication is linear in the accumulator, so it
does not change that discrepancy. The identity works for shared Miller
loops as well. The Python tests verify this exact raw-value relation at
the end of full loops, not just equality after final exponentiation.

This does **not** preserve the raw square or raw Miller value. Do not use it
inside the current bit-identical `pairing_miller_prec` API. Nor can it
replace Granger–Scott squaring in the final exponentiation: a general
Fp-scaled element need not remain cyclotomic.

## Operation accounting

Here M2=3Mp, S2=2Mp, and the implemented S4=3S2=6Mp.
“Linear work” includes modular additions/subtractions and coefficient
permutations; it is included in the measured execution times.

| Square | Multiplicative work | Leading Mp count |
| --- | --- | ---: |
| Production general square | 2M6 = 12M2 | 36 |
| Scaled asymmetric cubic control | 4S4 + M4 | 33 |
| Gaussian scaled square | 5S4 = 15S2 | 30 |
| Fourier scaled square | 5S4 = 15S2 | 30 |

Outside the five S4 calls, Fourier uses 19 Fp4 additions/subtractions
(including doublings), versus 23 for the Gaussian version, plus the cheap
`u` and `s` operations. It reduces the multiplication count, not every
category of operation. This is not a lower-bound claim.

With 64 general squares and 69 evaluated lines in the current schedule, a
single prepared replay has leading counts:

```text
8M2 control: 64*36 + 69*28 + 1 = 4237 Mp + 1 Ip
Fourier:     64*30 + 69*28 + 1 = 3853 Mp + 1 Ip
```

The extra `1 Mp + 1 Ip` is the online P normalization, present in both.
The 384Mp saving per shared loop explains why the percentage falls as the
number of terms increases. No final-exponentiation operations were removed.

## Validation

Performed successfully on the preserved code:

- Exact Gaussian and Fourier polynomial identities over `Z[i,a,b,c,s]`,
  using SymPy in `../inline_checks.py`.
- 1,009 Python field cases per curve on BLS12-381, BLS12-461, and BN-462,
  for each candidate. Includes zero, one, minus one, basis vectors, and
  full-width random values.
- Instrumented Python counts: exactly 15 S2 and zero general M2 per
  candidate square; 23 external A4 for Gaussian and 19 for Fourier.
- All four committed BLS12-381 pairing KATs, for both candidates, using the
  exact final exponent E.
- Shared products with 0, 2, 4, and 5 terms and explicit raw scale checks.
- Fourier replay cancellation check, e(P,Q)*e(-P,Q)=1.
- Python interpolation corruption control; checks also pass under `-O`.
- C: 4,000 full-width/edge cases for all four square implementations,
  checked against dense multiplication, with separate-output and in-place
  calls.
- C: three sets of eight points with full-width public scalar multiples;
  reduced-pairing comparison for n=0..8, all four methods.
- C: an additional exact final-exponent comparison for all four methods.
- Fresh standalone Release and ASan+UBSan builds with `-Wall -Wextra
  -Werror`; both research CTest tests pass in each build.
- Relocated legacy norm and primitive-ratio programs execute successfully.
- Reconstructed lattice unit tests and symbolic searches pass. The slow
  full-size lattice searches were **not** repeated this round.

Not performed: production integration, a complete repository test run,
constant-time leakage tests, C pairing validation on BN/BLS12-461, or
multi-architecture benchmarks. Python field checks on the other curves
are not a claim that this C replay, which lacks BN's final extra lines,
implements their pairings.

## Benchmark conditions

Apple M2 Max, Apple Clang 21.0.0, GMP 6.3.0; a fresh `-O3 -DNDEBUG`
research build of this repository. Two independent runs, each with:

- 31 samples per method/workload after two warmup blocks;
- 4,000 dependent calls per square sample;
- 30 complete replays per Miller or Miller-plus-final sample;
- all four methods interleaved in rotating order;
- no other research test/build running concurrently with timed runs.

Each replay includes one online inversion and one multiplication per P.
Timings exclude Q table construction/normalization, point-to-affine
conversion, subgroup checks, and the public API's chunk-accumulator wrapper.
Some samples contain noticeable OS noise. The small whole-pairing
differences, especially at eight terms, should be treated cautiously.

## Prior art and novelty boundary

- Devegili, Ó hÉigeartaigh, Scott, and Dahab,
  [Multiplication and Squaring on Pairing-Friendly Fields,
  ePrint 2006/471](https://eprint.iacr.org/2006/471):
  section 2 explicitly removes Toom-Cook denominators through subfield
  scaling. Section 4 gives both integer-point Toom-Cook-3x and CH-SQR3x.
  The implemented cubic control is from that established family.
- Chung and Hasan,
  [Asymmetric Squaring Formulae](https://www.lirmm.fr/arith18/papers/Chung-Squaring.pdf):
  established alternatives to symmetric polynomial squaring.
- Fourier/Toom-Cook interpolation and choosing cheap roots of unity are
  standard polynomial-arithmetic techniques.

The particular roots-of-unity schedule above was derived and tested here,
and is absent from the current production squaring. The targeted searches
did not establish whether this exact specialization has appeared elsewhere.
**Do not call it a new cryptographic algorithm or publication-level novel
result without a substantially deeper prior-art review.**

The useful finding is narrower: the earlier “squaring floor” conclusion was
false, and a concrete, reduced-pairing-compatible formula gives a repeatable
incremental improvement on this implementation.
