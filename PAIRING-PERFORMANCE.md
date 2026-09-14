# Why the pairing is slower than blst

Measured, not argued. Everything here was taken on one machine in one sitting:
an Intel Xeon at 2.80 GHz, blst `de54cd4` built with `-D__ADX__` (its fast
assembly path), ELiPS at `a206f81`, both `-O3`. Operation counts come from
callgrind, which counts calls exactly and is independent of which backend runs.

## The gap

| | ELiPS | blst | ratio |
|---|---:|---:|---:|
| pairing | 1060 µs | 403 | 2.6x |
| miller loop | 400 | 180 | 2.2x |
| final exponentiation | 425 | 239 | 1.8x |

## It is not the base field

This is the first thing to rule out, and it rules out cleanly.

| | ELiPS | blst |
|---|---:|---:|
| Fp multiply | **100.6 cycles** | 98.6 |
| Fp square | 102.4 | 101.0 |

Two percent apart. The `mulx`/`adcx`/`adox` assembly in
`src/arith/fp_mul_x86_64.S` is competitive with blst's. Nothing below the
extension tower explains a 2.2x gap, and no amount of further work on
`fp_mul` will close it.

One level up, the picture changes:

| | ELiPS | blst |
|---|---:|---:|
| Fp2 multiply | **372.8 cycles** | 322.3 |

16% slower, and that difference is structural rather than incidental.

## Where the work actually goes

Calls in one miller loop:

| | ELiPS | blst | ratio |
|---|---:|---:|---:|
| 384-bit multiplications | 8,550 | 6,867 | **1.25x** |
| Montgomery reductions | 8,550 | ~3,336 | **2.56x** |
| Fp-level add/sub | 27,929 | ~17,700 | **1.6x** |

Three separate problems. The reduction count is the largest and the most
interesting.

### Why ELiPS reduces 2.56x as often

`fp2_mul` in `src/arith/fpx.c` is Karatsuba, which is the right algorithm:

```c
fp_mul(t0, a[0], b[0]);
fp_mul(t1, a[1], b[1]);
fp_mul(s,  a0+a1, b0+b1);
```

Three multiplies rather than four. But each `fp_mul` is a *fused* multiply and
Montgomery reduction, so three multiplies means **three reductions**.

blst does not do that. Its top miller-loop function by instruction count is
`__mulx_384` — a raw 384×384 multiply with **no** reduction — and its Fp6
routines are named `mul_fp6x2` and `mul_by_xy0_fp6x2`, where the `x2` is
double width. It multiplies into 768-bit accumulators, adds and subtracts at
768 bits (`add_mod_384x384`, `sub_mod_384x384`), and reduces **once per output
coefficient**. That is lazy reduction, and it is the single biggest structural
difference between the two libraries.

The reduction is worth about half of a Montgomery multiply. Two independent
measurements agree:

- Counting instructions in `ROUND6` of `fp_mul_x86_64.S`: 29 in the multiply
  half, 25 in the reduction half — **46%**.
- Timing the halves separately in portable C: `mul_raw` 126.9 cycles, `redc`
  132.2, against `fp_mul_portable` at 267.7. The parts sum to the whole, and
  the reduction is **50%**.

So 5,200 surplus reductions is real work, not bookkeeping.

## What lazy reduction is actually worth — measured, not projected

A prototype Fp2 multiply with lazy reduction (three raw multiplies, two
reductions, 768-bit combination) was written and checked against the shipping
`fp2_mul` over 200,000 random inputs: **0 mismatches**.

Compared like with like — both sides built from the *same* portable
primitives, so the comparison measures the algorithm and not anyone's
assembly:

| | cycles |
|---|---:|
| eager (3 multiplies, 3 reductions) | 933.8 |
| lazy (3 multiplies, 2 reductions) | **862.1** |

**7.7% saved.** Not the 14% that removing one reduction in three suggests: the
768-bit adds and subtracts eat half the gain. That is the honest number, and
it is why lazy reduction at the Fp2 level alone is not the answer.

The saving grows at Fp6 and Fp12, where more products accumulate before a
reduction is forced — blst reduces once per Fp6 coefficient, amortising one
reduction over six multiplications instead of one over one. That is where its
2.56x comes from. But the 768-bit overhead grows too, and the Fp2 measurement
says to expect roughly half of the naive projection.

**One bound worth recording.** The first version of the prototype failed 4.3%
of random inputs. Both Karatsuba combinations can go negative, so a multiple
of p² is added to keep them non-negative — p² ≡ 0 mod p, so it changes
nothing. But the bias must also leave the accumulator **below p·R**, or the
single conditional subtraction at the end of REDC cannot land in [0,p). The
first attempt used p·2³⁸⁴, which *is* p·R, so the sum sat exactly at the limit
and the reduction came out short on the inputs that pushed it over. The
correct biases are p² for the real part (keeping it in [0, 2p²)) and 2p² for
the imaginary part ([0, 3p²)), both comfortably under p·R.

## The other lever, already prototyped in this repository

`tools/reference/normalized_miller_ref.py` implements a normalized,
subfield-scaled line multiplication. Running it:

```
PASS 1000 kernel checks; counted 8M2 and 1I2+342M2 normalization
C cost model (1M2=3Mp): 15M2 -> 8M2 per line; 1449Mp saved over 69 lines.
```

**15 Fp2 multiplications per line become 8.** Over the 69 lines in a
BLS12-381 miller loop that is 1,449 Fp multiplications saved out of ELiPS's
8,550 — **17% of the miller loop's multiplications**, from an algorithmic
change rather than an implementation one.

The catch is in the same output: the table normalization costs `1I2 + 342M2`,
which is 1,026 Mp plus an Fp2 inversion. For a single pairing against a fresh
G2 point that eats most of the 1,449 Mp saving.

**And this library's usage does not amortise it.** In `pairing_check` the G2
arguments are the signature and the message hash, both fresh per verification;
the fixed argument is `-g1`, which is in G1. So table reuse, which is what
makes prepared-line techniques pay elsewhere, does not apply to BLS
verification as this library performs it.

Two ways that could still be made to work, neither yet attempted:

1. The prototype's own header says preparation "deliberately uses slow affine
   arithmetic for readability; its timing is NOT a model of the C
   implementation's projective preparation." The 342M2 is an artifact of the
   Python, not a floor. A projective preparation could be much cheaper, and
   establishing what it actually costs is the first thing to do.
2. `elips_pairing_prec` already exists for callers who hold a fixed G2 point.
   Where a verifier checks many signatures against one public key, the table
   *is* reusable, and there the 8M2 kernel is close to free.

## What this does and does not add up to

The three counts — 1.25x multiplications, 2.56x reductions, 1.6x additions —
are together worth roughly 1.7 to 2x, which with blst's better instruction-level
parallelism accounts for the measured 2.2x. The map is complete enough to work
from.

But neither lever alone closes it:

- Lazy reduction: measured 7.7% at Fp2, plausibly 15–20% through Fp6/Fp12.
- The 8M2 kernel: 17% of miller multiplications, currently offset by a
  normalization cost this library's call pattern cannot amortise.

Both together, if the normalization can be made cheap, land somewhere near
30–35% off the miller loop. That takes 400 µs to about 270 — still 1.5x behind
blst's 180.

## The additions, traced

Every `fp_add` and `fp_sub` in the miller loop, attributed to its caller by
parsing callgrind's `cfn=`/`calls=` pairs inside each `fn=` block:

| `fp_add` 14,613 | | | `fp_sub` 13,316 | |
|---|---:|---|---|---:|
| from `fp2_add` | 8,148 | | from `fp2_mul` | 7,821 |
| from `fp2_mul` | 5,214 | | from `fp2_sub` | 4,296 |
| from `fp2_mul_xi` | 867 | | from `fp2_mul_xi` | 867 |
| from `fp2_sqr` | 384 | | from `fp2_sqr` | 192 |

**Every leaf routine is exactly at its formula minimum.** Karatsuba `fp2_mul`
needs 2 additions and 3 subtractions: 2,607 calls give 5,214 and 7,821.
`fp2_add` needs 2, `fp2_mul_xi` needs 1 of each for ξ = 1+u. There is no waste
at the Fp level and nothing to reclaim there.

So the additions are not a defect in their own right. They are a *consequence*
of how many Fp2-level operations the tower performs — `fp2_mul` alone creates
13,035 of the 27,929, **47%**, purely as Karatsuba overhead. Reduce the number
of Fp2 multiplications and the additions fall with them.

That reframes the question, and the answer is one function.

## The doubling step

`dbl_line` in `src/pairing/miller.c` costs, per call: **11 Fp2 multiplications
and 3 Fp2 squarings**. blst's `line_dbl` costs **3 multiplications and 8
squarings**. Both confirmed by attributing each library's Fp2 routines to
their callers; 64 doubling steps in ELiPS against 63 in blst.

An Fp2 multiplication is 3 Fp multiplications and an Fp2 squaring is 2, so:

| | per doubling step |
|---|---:|
| ELiPS `dbl_line` — 11M + 3S | **39 Mp** |
| blst `line_dbl` — 3M + 8S | **25 Mp** |

**1.56x, on the step the Miller loop runs 64 times.** Over the loop that is
2,496 Mp against 1,575 — a difference of **921 Mp, which is 55% of the entire
1,683 multiplication gap.**

The cause is visible in the source. `dbl_line` uses the Renes–Costello–Batina
*complete* addition formulas, and forms X³, X²Z and YZ² explicitly for the
line. Complete formulas handle every exceptional case, which is a real virtue
in general and an unnecessary one here: inside a Miller loop T is a point of
odd prime order r being doubled fewer than r times, so it is never the
identity and never equals its own negative. The dedicated formula is safe in
this context, and it trades multiplications for squarings — which is the
direction that pays, because a squaring is two thirds the cost.

Adopting a 3M + 8S doubling-and-line formula would be worth, arithmetically:

- **−896 Fp multiplications**, −10.5% of 8,550
- **−1,600 Fp additions and subtractions**, −5.7% of 27,929

roughly **8% off the miller loop**, from one function, with no change to the
field layer and no new representation. It attacks both gaps at once, which
neither of the other two levers does.

## Revised ordering

| | worth | risk |
|---|---|---|
| `dbl_line` formula, 11M+3S → 3M+8S | ~8% of the miller loop | contained: one function, existing KATs, complete-vs-dedicated needs the subgroup argument written down |
| 8M2 line kernel (`normalized_miller_ref.py`) | 17% of multiplications | normalization cost this call pattern cannot amortise; needs a projective preparation first |
| lazy reduction through Fp6/Fp12 | 7.7% measured at Fp2 | touches the whole tower; adds additive work to remove reductions |

The doubling step is now the first thing to do: it is the largest single
identified item, the most contained, and the only one that reduces
multiplications and additions together.

**It still does not close 2.2x.** All three together land somewhere near 30%,
taking 400 µs to roughly 280 against blst's 180. The remaining 762 Mp of the
multiplication gap is not yet attributed — `fp12_mul_sparse035` at 9 Fp2
multiplications per call and `fp6_mul` at 6 are the places left to look.

## Reproducing any of this

```bash
git clone --depth 1 https://github.com/supranational/blst && (cd blst && ./build.sh)
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/bench_BLS12_381 --reps 21          # ELiPS side
```

Operation counts come from `valgrind --tool=callgrind` on a program that
performs exactly one miller loop, parsing `calls=` out of the output. Call
counts are backend-independent, which matters because valgrind hides ADX from
CPUID and both libraries fall back to their portable paths under it.
