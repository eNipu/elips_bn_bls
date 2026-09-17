# Why the pairing is slower than blst

Measured, not argued.

**The timing table in this document was wrong and has been replaced.** The
earlier figures (pairing 1060 us against 403, miller 400 against 180) were not
reproducible. Re-measuring both libraries in one sitting through one harness
put blst's miller loop at 304 us, not 180, while ELiPS moved only from 400 to
468. A uniform host change cannot scale one library by 1.67x and the other by
1.15x, so the two halves of the old table were not taken under the same
conditions, whatever the document claimed. The ratios below replace them.

What survives unchanged: every operation count. Those come from callgrind,
which counts calls exactly and does not care which backend runs or how fast
the host is.

## How the numbers below were taken

Intel Xeon at 2.80 GHz. blst `8921f76` built with `-D__ADX__`, ELiPS at
`4f52583`, both `-O3`, both linked into one binary (`bench/cross_blst.c`) that
times every case through the same function-pointer table and the same
interleaved loop. Medians over 15 to 31 reps.

Three things were checked before any number here was believed:

- **The harness detects a regression.** Five extra `fp2_mul` were planted in
  `dbl_line`, which is arithmetically +11.2%. Alternating clean and planted
  binaries three times reported +9.3%, +9.1%, +8.9%. The plant was reverted.
- **blst is on its fast path.** The `-D__ADX__` build selects the ADX assembly
  at compile time and would fault on a CPU without it; it ran. `objdump` finds
  1,017 `mulx`/`adcx`/`adox` sites in the archive.
- **The ratios are stable.** Three independent runs gave miller 1.518x,
  1.504x, 1.530x. Medians agree to 2% even where a single descheduled sample
  pushes the reported spread above 40%.

## The gap

| | ELiPS | blst | ratio |
|---|---:|---:|---:|
| miller loop | **433 us** | 310 | **1.40x** |
| final exponentiation | 522 | 409 | 1.22x |
| pairing, arithmetic only | 963 | 716 | 1.34x |
| pairing, validation included | 1232 | 837 | 1.47x |
| G1 subgroup check | **61** | 50 | **1.24x** |
| G2 subgroup check | **103** | 61 | **1.69x** |

The miller row is after the doubling-step change below, and the two subgroup
rows are after the scalar-multiplication change after it. Before those the
loop was 468 us and 1.52x, and the checks were 96 and 161 us, 1.85x and
2.49x. The pairing row carries both.

Two rows, not one, because `elips_pairing` validates both input points and
`blst_miller_loop` plus `blst_final_exp` do not. Comparing the first against
the second is comparing an API to an algorithm. Validation is **22.9%** of the
ELiPS pairing and 14.4% of blst's.

That also answers the question of whether the subgroup checks are merely an
API difference. They are not. blst performs the same two checks and is 2.49x
faster at the G2 one, which is the single largest ratio in this table.

## The counts over-predict the gap

This is the most useful thing the re-measurement changed.

The operation counts below say ELiPS performs 1.25x the multiplications,
2.56x the reductions and 1.6x the additions. The old document reasoned that
those were "together worth roughly 1.7 to 2x" and that blst's better
instruction-level parallelism explained the rest of a 2.2x gap.

The measured miller gap is 1.52x, which is **less** than the operation counts
predict. ELiPS is not losing time per operation; it is losing operations. Its
per-operation efficiency is at least blst's, which is consistent with the Fp
multiply being at parity.

That makes the levers below more valuable than they looked, not less. Work
that removes operations should convert to time roughly one for one, and the
starting point is 1.52x rather than 2.2x.

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

## Lazy reduction, measured and rejected

**Do not build this.** The section below it records a 7.7% prototype gain at
Fp2. That number was measured with portable primitives on both sides, and it
does not survive contact with the shipping assembly. Reproduce any of this
with `bench/lazy_probe.c`.

The two primitives a lazy tower needs were written in assembly, as the halves
of `ROUND6` separated, and checked against `fp_mul` over random inputs:

| | ns | share of a fused multiply |
|---|---:|---:|
| fused `fp_mul` (asm, multiply and reduce) | 36.05 | 100% |
| `mulw` raw 384x384 product | 17.74 | **49%** |
| `redcw` 768 -> 384 reduction | 22.63 | **63%** |

**The halves sum to 112% of the whole they replace.** The fused CIOS routine
interleaves both carry chains across its multiply and reduction halves; split
apart, each has a tighter dependency chain and less to overlap. That 12% is
paid on every multiplication and it is what the reduction saving has to beat.

At Fp2 it does not even start: three fused multiplies cost 106.1 ns, and three
raw products plus two reductions plus the wide combination cost 106.0 ns.
**0.1%.**

### Fp6 is where the amortising happens, and it still fails

Reductions amortise 3:1 at Fp6, one per output coefficient against one per
multiplication. A complete lazy Fp6 multiply was written and checked against
`fp6_mul` over 20,000 random inputs with 0 mismatches. Every wide value is kept
in `[0, p*R)` by conditionally adding or subtracting `p*R`, which is what
blst's `add_mod_384x384` does, and which is exact because
`REDC(w + p*R) = w*R^-1 + p == w*R^-1 (mod p)`. `p*R` has six zero low limbs,
so that conditional only ever touches the top half: a 768-bit modular add
costs **1.68 ns**.

| `fp6_mul` | ns | |
|---|---:|---:|
| eager, shipping | 1030.5 | |
| lazy, ideal: 18 `mulw` + 6 `redcw` + 40 wide adds | **522.2** | **-49.3%** |
| lazy, measured in C | **1801.5** | **+74.8%** |

**The algorithm is worth 49%. The C implementation gives away 1,279 ns, more
than the entire eager multiply.**

The overhead is memory. A 768-bit intermediate is 96 bytes and does not live in
registers across a C function boundary: `fp2_mulw` measures 153.4 ns where its
own parts sum to 58.3, the difference being the spill and reload of values
`mulw` has just written. One `fp6_mul` moves about two kilobytes of stack that
the eager version never touches, because 384-bit values fit in registers.

blst does not pay this because `mul_fp6x2` **is assembly**: its double-width
intermediates never leave registers inside the routine.

### What this actually costs, and why it stops here

Lazy reduction in this library is not "carry the double-width representation up
through Fp6 and Fp12". It is **write the Fp6 and Fp12 multiplication in
assembly**, for every curve, with a C fallback that must stay eager because in
C the transformation is a large loss. That forks the tower rather than
extending it, and it is a different and much larger project than the one that
was scoped.

The gate was set in advance at 6% and the measurement is -74.8%. Stopped.

## What the earlier prototype measured — superseded, kept for the record

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

1. ~~The 342M2 is an artifact of the Python, not a floor.~~ **Wrong, and the
   arithmetic says so exactly.** Montgomery batch inversion of 69 values costs
   3(69-1) = 204 M2, and beta = c/a with gamma = b/a costs 2 per line = 138
   M2. 204 + 138 = 342. The figure is the cost of normalizing 69 lines and is
   independent of language and of coordinate system. The prototype's caveat
   about slow affine arithmetic applies to the point arithmetic that produces
   the lines, which is separate and genuinely improvable. The two were
   conflated. Net saving on a single pairing is about 423 Fp multiplications,
   4.9%, with an Fp2 inversion still owed.
2. `elips_pairing_prec` already exists for callers who hold a fixed G2 point.
   Where a verifier checks many signatures against one public key, the table
   *is* reusable, and there the 8M2 kernel is close to free.

## What this does and does not add up to

The three counts, 1.25x multiplications, 2.56x reductions and 1.6x additions,
over-predict the measured 1.52x rather than under-predicting it. See "The
counts over-predict the gap" above. The map is complete enough to work from,
and it is a map of operations to remove rather than of time to claw back.

But neither lever alone closes it:

- Lazy reduction: measured 7.7% at Fp2, plausibly 15–20% through Fp6/Fp12.
- The 8M2 kernel: 17% of miller multiplications, currently offset by a
  normalization cost this library's call pattern cannot amortise.

Both together land somewhere near 25% off the miller loop, and the
normalization cannot be made cheap: see the correction below. That takes 468
us to about 350 against blst's 304, which is 1.15x.

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

## The doubling step, and what replaced it

**Done.** `dbl_line` cost 11 Fp2 multiplications and 3 squarings per call, 39
Fp multiplications. blst's `line_dbl` costs 3 and 8, which is 25. It now costs
**4 and 6, which is 24**, and the Miller loop is below blst on this step.

Two changes, and the first is the one that matters.

**The line is scaled by one more factor of Z, and the curve equation pays for
the constant term.** Every scaling in this loop is free because any Fp2 factor
dies in the final exponentiation. Taking one more factor of Z turns the
constant coefficient into `(3X^3 - 2Y^2 Z)/Z`, and `Y^2 Z = X^3 + b'Z^3` gives

```
(3X^3 - 2Y^2 Z)/Z = (3(Y^2 Z - b'Z^3) - 2Y^2 Z)/Z = Y^2 - 3b'Z^2
```

The cubing is gone. What remains is one squaring the doubling already needs
and one multiplication by the constant 3b'. That was the whole cost of the old
line: `X^3`, `X^2 Z` and `Y Z^2` were three multiplications spent on terms the
curve equation makes unnecessary.

**The doubling drops the complete formulas for the dedicated ones.** The
Renes-Costello-Batina formulas handle every exceptional case, which is a real
virtue in general and an unnecessary one here. The dedicated formulas for
`Y^2 Z = X^3 + b'Z^3` are undefined only at `Z = 0` and `Y = 0`, and T reaches
neither: it is always `[v]Q` for Q of odd prime order r, so it is never
2-torsion, and `ELIPS_LOOP` is a NAF whose prefix values satisfy
`|v| > 2^(k-1)` and stay far below r, so it is never the identity. The
argument is written out in the source, because a wrong one here fails silently
rather than loudly.

The usual halvings in these formulas are avoided by scaling all three output
coordinates by 4, which is free in projective coordinates.

Measured by callgrind over one Miller loop:

| | before | after | delta |
|---|---:|---:|---:|
| `fp_mul` | 8,566 | 7,606 | **-960, -11.2%** |
| `fp_add` | 14,618 | 13,978 | -640 |
| `fp_sub` | 13,325 | 12,557 | -768 |
| `fp2_mul` | 2,609 | 2,161 | -448 |
| `fp2_sqr` | 192 | 384 | +192 |

-960 is exactly 15 Fp multiplications times 64 doubling steps, and the
squaring count rises by exactly 3 per step. The counts are the formula, with
nothing unaccounted.

Measured by alternating the two binaries, which is the only comparison this
repository trusts:

| | before | after | |
|---|---:|---:|---:|
| BLS12-381 miller | 470.1 us | 431.5 | **-8.5%, 100% confidence** |
| BLS12-461 miller | 891.9 | 809.4 | **-9.5%** |
| BN-462 miller | 1392.0 | 1263.4 | **-9.0%** |
| BLS12-381 pairing | 1266.4 | 1232.2 | -2.7% |
| BLS12-461 pairing | 2236.2 | 2159.7 | -3.8% |
| BN-462 pairing | 2969.0 | 2842.1 | -4.4% |

All three curves improve and not one row on any of them reports `slower`,
which matters as much as the miller rows: `hash_to_g1` moves 201.9 to 201.8
and `final_exp` 522.4 to 522.3, so nothing outside the doubling step was
disturbed.

**11.2% of the multiplications bought 8.5% of the time.** The shortfall is
worth naming rather than rounding away: the formula trades multiplications for
squarings and additions, and the additions are not free. Against blst the
miller loop goes from 1.52x to **1.40x**.

One warning about method. Comparing the new numbers against the recorded
`bench/baseline.json` instead would have reported `hash_to_g1` 25% slower and
the BLS12-461 scalar multiplications 13 to 15% slower, none of which this
change can touch. That is host drift between two recorded runs, and it is
exactly why `compare.py --ab` alternates binaries rather than comparing files.

## The subgroup checks

**Done.** The G2 check was the largest single ratio in the table at 2.49x. It
is now 1.69x, and the G1 check 1.24x from 1.85x.

The test itself did not change, and that is the point. Both libraries use the
same one, Scott, https://eprint.iacr.org/2021/1130:

```
G1   sigma^2(P) == [-z^2]P
G2   psi(P)     == [z]P
```

blst's `POINTonE2_in_G2` and `ep2_in_subgroup` are the same three lines. The
entire 2.49x was in how `[z]P` gets computed.

ELiPS used `ep2_mul`, the general constant-time fixed-window ladder: a
16-entry table, and every window scans the whole table under a mask so that no
memory address depends on the scalar. That is the right routine for a caller's
secret scalar and the wrong one here, because **z is a curve parameter
compiled into `fp_params.h`**. Nothing about it is secret, and it is chosen to
be sparse: |x| on BLS12-381 is `0xd201000000010000`, six set bits in
sixty-four.

`ep_mul_pubconst` and `ep2_mul_pubconst` expand that constant into
non-adjacent form and skip the zero digits. The control flow depends on the
constant alone, so the routines remain **constant time in the point**, which
is the property that matters when the point is what an attacker supplies. They
keep the complete addition formulas, which unlike the Miller loop is not
something to give up here: the caller is validating a point it does not trust,
so every exceptional case has to work rather than be argued away.

**NAF and not plain binary, and that distinction was not free.** Binary was
written first and measured first, and it made two of the three curves slower:
BLS12-461's `ep2_in_subgroup` by 6.9% and BN-462's by 11.7%. The seeds are
sparse in signed-digit form and dense in bits:

| constant | bits | popcount | NAF |
|---|---:|---:|---:|
| BLS12-381 seed | 64 | 6 | 6 |
| BLS12-461 seed | 77 | **43** | **3** |
| BN-462 seed | 115 | **101** | **4** |
| BN-462 6x^2 | 231 | **106** | **18** |

BLS12-461 pays 42 additions in binary and 2 in NAF. These seeds are chosen
near powers of two, which is precisely the shape plain binary represents
worst, and BLS12-381 was the one curve sparse enough in bits to hide it.

`ep2_in_subgroup`, counted on all three curves:

| | `fp_mul` | `ep2_add` | `fp2_cselect` |
|---|---|---|---|
| BLS12-381 | 3,017 -> **2,047** (-32%) | 23 -> **5** | 768 -> **0** |
| BLS12-461 | 3,637 -> **2,282** (-37%) | 27 -> **2** | 960 -> **0** |
| BN-462 | 9,527 -> **7,239** (-24%) | 65 -> **17** | 2,784 -> **0** |

Every addition count is exactly the NAF weight less the leading digit: 6-1,
3-1, 18-1. The counts are the representation, with nothing unaccounted.
Measured on BLS12-381:

| | before | after | |
|---|---:|---:|---:|
| `ep2_in_subgroup` | 158.7 us | 102.4 | **-35.8%** |
| `ep_in_subgroup` | 96.2 | 61.5 | **-36.3%** |
| `pairing` | 1239.8 | 1143.2 | **-7.8%** |
| `bls_verify` | 3227.4 | 2933.5 | **-9.1%** |
| `miller` | 434.2 | 434.1 | unchanged, as it should be |

and on the other two, which is where the first attempt at this failed:

| | before | after | |
|---|---:|---:|---:|
| BLS12-461 `ep2_in_subgroup` | 308.5 us | 186.4 | **-39.3%** |
| BLS12-461 `pairing` | 2138.1 | 1963.4 | **-8.8%** |
| BN-462 `ep2_in_subgroup` | 813.9 | 596.3 | **-26.6%** |
| BN-462 `pairing` | 2851.3 | 2623.9 | **-8.0%** |

No row on any curve reports `slower`.

Nothing here required a new exactness argument, and that is worth stating
plainly rather than glossing: the test is untouched, so the gcd conditions
`tools/reference/subgroup_ref.py` derives and `gen_params.py` re-asserts still
cover it. Only the route to `[z]P` changed, and that was checked against
`ep_mul` over 200,000 random point-and-scalar pairs plus the cases a
double-and-add loop gets wrong: k = 0, k = 1, a declared bit length longer
than the value, and the point at infinity.

**What is left.** 1.69x on G2 is no longer the biggest ratio but it is not
parity either. The remainder is formula cost: ELiPS doubles and adds with the
complete formulas throughout, blst uses dedicated ones. The same trade as the
Miller loop, but the safety argument is harder because the input is untrusted
rather than a known multiple of a prime-order point, so it is not taken here.

## The final exponentiation, traced

It had never been traced. It is 42% of the pairing and 1.22x off blst, and the
ranking of everything else was being made without it. The result is the
opposite of what was expected.

**ELiPS already does less multiplication work than blst here.** The entire gap
is reductions.

| one final exponentiation | ELiPS | blst | |
|---|---:|---:|---|
| cyclotomic squarings | 321 | 315 | parity |
| Fp2 squarings | **2,234** | 2,841 | **ELiPS 21% fewer** |
| Fp2 multiplications | 678 | 646 | parity |
| 384-bit multiplications | **6,568** | 7,707 | **ELiPS 15% fewer** |
| Montgomery reductions | 6,568 | **4,225** | **ELiPS 55% more** |

The chain, attributed:

```
pairing_final_exp_fast    1        fp12_exp_param     5
fp12_sqr_cyc             81        fp12_sqr_cyc_run  10
fp12_mul                 34        fp6_mul          104
fp12_frobenius            3        fp12_conj          9
fp2_sqr               2,234        fp2_mul          678
```

`2,234 x 2 + 678 x 3 = 6,502` against an actual `fp_mul` of 6,568. The 66
remaining are the 11 `fp_inv` and 11 `fp2_inv` of Karabina decompression and
their helpers. **The map is complete**, not complete to within 10%.

### Why ELiPS is ahead on multiplications

Karabina. `fp12_exp_param` walks the exponent in runs of consecutive zero
digits and squares a compressed representation across any run of 16 or more,
which is cheaper per squaring but costs an Fp2 inversion to undo. Per
exponentiation the runs are 2, 2, 3, 9, 32, 16: the 32 and the 16 go through
`fp12_sqr_cyc_run`, the rest through plain Granger-Scott squaring. Over five
exponentiations that is 240 compressed and 80 uncompressed, plus one
standalone, so 321 cyclotomic squarings against blst's 315 with the same
algebra but no compression. Same squarings, 21% fewer Fp2 squarings to do
them.

**A candidate checked and rejected.** Lowering `ELIPS_KARABINA_MIN_RUN` below
16 to catch the run of 9 does not pay: nine compressed squarings save roughly
45 Fp2 squarings, and the decompression they would need costs an Fp2 inversion,
which is the divstep at `3 * FP_BITS` iterations. The threshold is where it
should be.

### What the gap actually is

Both halves of a Montgomery multiply cost about the same, measured earlier in
this document at 46% and 50% by two independent methods. Counting each as one
unit:

| | raw multiplies | reductions | total |
|---|---:|---:|---:|
| ELiPS | 6,568 | 6,568 | **13,136** |
| blst | 7,707 | 4,225 | **11,932** |

1.10x against a measured 1.22x, the remainder being additions and blst's
instruction-level parallelism. And the counter-factual is the useful part:

**If ELiPS reduced once per output coefficient as blst does, its final
exponentiation would do about 10,768 units against blst's 11,932. It would be
ahead by 10%.**

### The hard-part chain is already blst's, checked and closed

The obvious follow-up to the above is that ELiPS's hard-part chain might simply
be longer than blst's. It is not.

| | exp by x | fp12 multiplies | cyclotomic squarings |
|---|---:|---:|---:|
| ELiPS | 5 | **34** | 321 |
| blst | 5 | 35 | 315 |

**One fewer multiplication than blst, six more squarings.** The six are the
3*lambda convention this library computes by design, worth about 1.6%. Both
invert once in the easy part. There is no shorter chain to adopt here, so the
1.22x is reductions and only reductions, which is what the section above
already concluded from the operation counts. Recorded so the chain is not
re-examined a third time.

### The ranked list

1. **Lazy reduction, and nothing else.** 2,343 surplus reductions is the whole
   of it. This is the same lever as the Miller loop, where the surplus ratio is
   larger still at 2.56x.
2. There is no second item. Karabina already beats what blst does on the
   squaring side, the squaring counts are at parity, and the multiplication
   count is in ELiPS's favour.

**This closes the question the trace was run to answer.** The final
exponentiation was the last place large enough to hide a novel algorithmic
result, and it does not contain one: ELiPS is already ahead of blst on
algorithm here and behind only on representation. Every remaining microsecond
in the pairing, in both halves, now runs through lazy reduction.

## hash_to_g2, the largest ratio in the library

Issue #42 measured `hash_to_g2` at 6.7x off blst, on 18,982 `fp_mul`. Traced,
it splits cleanly:

| | `fp_mul` | share |
|---|---:|---:|
| `fp2_exp`, inside `fp2_sqrt` | 11,580 | **61%** |
| cofactor clearing, the GLV ladder | 6,924 | 36% |
| the SSWU map itself | ~480 | 3% |

**Cofactor clearing is already at parity** and is not the problem: 132
`ep2_dbl` against blst's 127. The whole gap is the square root.

### What was fixed

`fp2_exp` was plain binary square-and-multiply, branching on every exponent
bit. Its two call sites are both in `fp2_sqrt`, with exponents `(p-3)/4` and
`(p-1)/2` derived from the modulus and public, and those exponents are dense:

| | bits | popcount | non-zero 4-bit windows |
|---|---:|---:|---:|
| BLS12-381 `(p-3)/4` | 379 | **228** | 92 |
| BLS12-461 `(p-3)/4` | 459 | 201 | 98 |
| BN-462 `(p-3)/4` | 460 | 208 | 79 |

**The same shape as the subgroup checks**: a dense exponent walked one bit at a
time. A 4-bit window turns 228 multiplications into 92 plus 7 to build the
table. The squarings are untouched and they are two thirds of the cost, so this
is worth about a quarter of an exponentiation, not half.

It changes no security property. The routine already branched on the exponent;
indexing a table by a window of a public exponent leaks nothing the branch did
not.

| `hash_to_g2` | before | after | |
|---|---:|---:|---:|
| `fp_mul` | 18,982 | 15,870 | **-16.4%** |
| `fp2_mul` | 4,198 | 3,150 | -25.0% |
| `fp_sub` | 16,973 | 13,845 | -18.4% |

Measured, median of three alternating runs:

| | before | after | |
|---|---:|---:|---:|
| `hash_to_g2` | 716.2 us | 633.5 | **-16.8%** |
| `bls_sign` | 1111.5 | 916.1 | **-17.1%** |
| `bls_verify` | 2474.8 | 2213.7 | **-10.9%** |
| `miller`, `final_exp`, `pairing` | | | unchanged |

and on the other two curves, exact counts rather than timings because the
timing noise on that run exceeded the effect:

| `hash_to_g2` | `fp_mul` before | after | |
|---|---:|---:|---:|
| BLS12-381 | 18,982 | 15,870 | **-16.4%** |
| BLS12-461 | 26,965 | 23,521 | **-12.8%** |
| BN-462 | 28,471 | 24,115 | **-15.3%** |

No row on any curve reports `slower`. The A/B on the two larger curves put
`hash_to_g2` at -17.2% and -13.6%, matching the counts, but at 62% and 75%
confidence: within each of those runs every other row drifted together, by -1
to -5% on BLS12-461 and +0.4 to +3.2% on BN-462, and `hash_to_g2` stands
clearly outside that drift in both. The counts are the stronger evidence here
and they are backend independent.

`bls_sign` moves almost one for one with `hash_to_g2`, which is the check that
the attribution was right: signing is a hash to G2 and a scalar multiplication.

### What remains, and it is most of it

**blst performs ONE exponentiation per map. ELiPS performs four.**

`map_sswu` computes `sqrt(gx1)` and `sqrt(gx2)` and discards one, deliberately:
choosing between them would branch on the message, which RFC 9380 treats as
secret. Each `fp2_sqrt` is then two exponentiations, because the complex method
needs `a^((p-3)/4)` and `(1+alpha)^((p-1)/2)`. Two maps, two roots, two
exponentiations: eight.

blst uses `sswu_opt` from the CFRG draft, where one call returns both the root
and the square/non-square flag, and the other branch is `y2 = y1 * u^3`. It
also carries `(xn, xd)` as a fraction and never inverts; ELiPS inverts three
times per map.

The obvious shortcut does not apply here, and it is worth recording why rather
than leaving someone to rediscover it. The identity

```
g(x2) == (Z u^2)^3 g(x1)
```

does hold, checked on 400 random inputs. But deriving the second root from the
first needs `sqrt(-Z^3)` as a constant, and on this curve **neither `Z^3` nor
`-Z^3` is a square**, because `p^2 == 9 mod 16` rather than `3 mod 4`. That is
exactly why blst passes `recip_ZZZ` and `magic_ZZZ` tables into
`recip_sqrt_fp2`. Adopting it means implementing the 9-mod-16 square root, not
adding a constant.

### One exponentiation instead of four: done

The map now computes **one** square root, not two, and that root costs one
exponentiation rather than two. `fp2_sqrt` is no longer on the G2 path at all.

`q = p^2` is 9 mod 16 on all three curves, so the 2-Sylow of Fp2* has order 8.
For `y = g(x1)^((q+7)/16)`, the value `t = y^2/g(x1)` is an 8th root of unity
and its index `k` decides everything:

```
k even   g(x1) is the square, and (y * zeta^(-k/2))^2    == g(x1)
k odd    g(x2) is the square, and (y * FIX[k] * u^3)^2   == g(x2)
```

The odd case works because `g(x2) = (Z u^2)^3 g(x1)`, so the missing factor is
`sqrt(Z^3 zeta^-k)`. It exists precisely because `Z^3` and an odd power of
`zeta` are both non-squares and a product of two non-squares is a square. That
is the shortcut the paragraph above says does not exist in the simple form: it
does exist, once the correction is allowed to depend on `k`.

`t` is never formed, because that would need an inversion. The loop compares
`y^2` against `g(x1) * zeta^k` for the eight `k`, which is eight
multiplications against a divstep, and selects under a mask so `k` stays
hidden. Nothing branches on the message.

`gen_h2c_params.py` derives `zeta`, `ZETA[8]` and `FIX[8]` from the curve
parameters and **re-checks the whole table against the map over 32 random
inputs before emitting the header**, so a wrong constant fails generation
rather than silently producing points that are on the curve, in the group, and
disagree with every other implementation.

| `hash_to_g2` | before P8 | after P8 | after P9 |
|---|---:|---:|---:|
| `fp_mul` | 18,982 | 15,870 | **11,558** |
| `fp2_sqr` | 3,064 | 3,080 | **1,536** |
| `fp2_exp` | 8 | 8 | **2** |
| `fp2_sqrt` | 4 | 4 | **0** |

**-27.2% against P8, -39.1% against where issue #42 started.** Against blst's
4,919 multiplications the ratio is now 2.35x, from 3.86x.

Measured:

| | before | after | |
|---|---:|---:|---:|
| `hash_to_g2` | 750.7 us | 580.3 | **-23.3%** |
| `bls_sign` | 1023.1 | 847.6 | **-16.7%** |
| `bls_verify` | 2430.3 | 2257.0 | -6.4% |
| `hash_to_g1`, `miller`, `pairing` | | | unchanged |

Median of three alternating runs, every `hash_to_g2` and `bls_sign` row at
100% confidence and `miller` and `pairing` never leaving the noise.

R11 is exact here rather than statistical: BLS12-461 and BN-462 use
Shallue-van de Woestijne, not SSWU, and their `hash_to_g2` operation counts are
**byte-identical** before and after, 23,521 and 24,115. BLS12-381's G1 also
keeps the two-root form, because Fp is 3 mod 4 and its 2-Sylow has order 2, so
the eight-element table does not apply.

### What is left in hash_to_g2

Cofactor clearing is now the largest part at 6,924 of 11,558, and it is already
at parity with blst in point operations. The rest of the gap to blst is the
same lazy reduction that P6 rejected: blst's Fp2 operations are cheaper, not
fewer. The remaining SSWU item is blst's fraction form, which carries
`(xn, xd)` and never inverts, against the three Fp2 inversions still here.

**The SvdW curves have the bigger untouched item.** BLS12-461 and BN-462
compute *three* square roots per map, six exponentiations against BLS12-381's
one, which is why their counts are twice BLS12-381's. Nothing here addresses
that.

### A note on the machine, and why the README was not updated

Everything above is an A/B on one host in one sitting, which is what makes the
percentages trustworthy. The absolute microseconds are not comparable across
the other tables in this document: the container moved from a 2.80 GHz Xeon to
a 2.10 GHz one partway through, so `hash_to_g2` reads 716 us here and 1010 us
in `bench/baseline.json` for code that differs by this one commit.

`bench/baseline.json` and the README table were deliberately **not**
regenerated. The README says "both columns on that one machine", and its
WebAssembly column was measured on the 2.80 GHz host. There is no `emcc` in
this container and the checked-in `bindings/js/dist/elips_wasm.mjs` predates
this change, so a regenerated native column would have been paired with a
WebAssembly column from a different machine and a different commit, under a
sentence claiming otherwise.

The record is therefore a consistent snapshot at `dd63ab1` on the 2.80 GHz
host, and it understates current signing by about 17%. Re-measuring **both**
columns together, on one machine, is what fixes it.

## The line multiply: both avenues closed

`fp12_mul_sparse035` applies the line to the accumulator 69 times per Miller
loop and is the largest single item in it. Two ways to make it cheaper were
checked. Neither survived.

**Karatsuba on the middle term: 207 fewer multiplications, no time.** The term
`f1 * (0, c3, c5)` was schoolbook, six Fp2 multiplications where Karatsuba
needs five. It was implemented, and the identity is exact with no scaling
freedom used, so the pairing output stayed byte-identical: digest
`76da830b565dbafb` over 3,000 inputs.

| miller loop | before | after |
|---|---:|---:|
| `fp_mul` | 7,606 | 7,399 (-207) |
| `fp2_mul` | 2,161 | 2,092 (-69) |
| `fp_sub` | 12,557 | 12,626 (+69) |

2.7% of the loop's multiplications, and the measured change was **-0.1% on the
assembly path and -0.1% on the portable one**, both reported `unchanged` across
repeated runs. Not the noise floor hiding a win: two tight runs agreeing near
zero.

The six schoolbook products are mutually independent and fill issue slots that
the Karatsuba chain then serialises, since `p4` waits on two sums and the last
coefficient waits on `p4`. It was reverted. **This is the same lesson as the
lazy reduction: an operation count is not a cost**, and it is the second time
in this document that removing arithmetic bought nothing because of what it did
to the dependency graph. The portable path was measured precisely because
multiplication is roughly twice as expensive there relative to addition, which
is where the trade should have paid if it paid anywhere.

**Rescaling so that yP becomes 1 is the 8M2 kernel, already rejected.**
`miller.c` carried a note that the legacy code reached two non-trivial
coefficients this way, "worth copying if the final measurement asks for it".
Making the w^0 coefficient 1 means dividing each line by its own w^0
coefficient, which varies per line, so it is a per-line inversion. That is
exactly the batch-normalised kernel whose `1I2 + 342M2` normalisation costs
more than it saves on the single-pairing path. The note now says so rather than
inviting a third attempt.

## The assembly path, as far as it got

P6 concluded that lazy reduction could not be reached from C. That was
confirmed end to end here, in the shipping code path rather than in a
microbenchmark: wiring P6's verified lazy Fp6 into `fpx.c` and running the
normal benchmark gives

```
miller +22.5%   final_exp +10.5%   pairing +13.6%   gt_exp +43.0%    all SLOWER
```

**But the obstacle is one routine, and writing it in assembly recovers most of
the loss.** The Karatsuba for a wide Fp2 product, done as three calls to the
raw-product routine plus a C combination, measures 128 ns -- against 137 for
the eager `fp2_mul` it is supposed to beat, even though the three raw products
alone are 52. The rest is call overhead (five callee-saved registers saved and
restored three times) and three 96-byte buffers stored and immediately
reloaded.

`bench/fp2_mulx2_x86_64.S` does the whole thing in one routine, which is what
blst's `mulx_382x` is for:

| | ns |
|---|---:|
| wide Fp2, C calling the raw product | 128.5 |
| wide Fp2, **one assembly routine** | **83.0** |
| shipping eager `fp2_mul` (reduced output) | 137.1 |

**0.61x the eager multiply, and its output is not even reduced yet.** Checked
against the C version over random inputs and the zero and one edge cases, and
all 59 tests pass with it driving a lazy Fp6 inside the library.

End to end that turns the loss into a much smaller one:

| | lazy Fp6 in C | with the assembly multiply |
|---|---:|---:|
| miller | +22.5% | **+7.2%** |
| pairing | +13.6% | **+3.9%** |
| final exponentiation | +10.5% | +3.0%, `unchanged` |
| `gt_exp` | +43.0% | +12.8% |

**Fifteen points recovered on the miller loop from one routine.** What remains
is arithmetic, not mystery: six wide multiplies at 83 ns plus six reductions at
23 is 636 ns against the eager Fp6 multiply's 989, which should be 36% faster.
The Fp6-level combination -- twenty-two 768-bit conditional add and subtract
operations over eight 192-byte temporaries, still in C -- is costing more than
350 ns, more than the reductions save.

### Then the combination, and the surprise in it

With the assembly multiply the lazy Fp6 was still 7.2% behind. Instruction
counts located the rest exactly: `fp6_mul`'s own cost, excluding everything it
calls, was **1,030,704 instructions per miller loop** against 950,328 for all
six wide multiplies together. The 768-bit combination was costing more than the
multiplications it existed to enable.

It was not the algorithm. It was `unsigned __int128` borrow extraction. Writing
the same twenty-two wide add and subtract operations with `_addcarry_u64` and
`_subborrow_u64`, which map to `adc` and `sbb` directly, cut that self cost to
**431,233** and moved the whole tower:

| | lazy Fp6 in C | + assembly multiply | + adc/sbb intrinsics |
|---|---:|---:|---:|
| miller | +22.5% | +7.2% | **+0.2%, unchanged** |
| pairing | +13.6% | +3.9% | +0.4%, unchanged |
| final exponentiation | +10.5% | +3.0% | +0.6%, unchanged |
| `gt_exp` | +43.0% | +12.8% | **-1.4%** |

**From 22.5% behind to parity, with `gt_exp` slightly ahead.** All 59 tests
pass with the lazy tower driving the real library.

### Finished, and it ships

The combination is now assembly too: `src/arith/fp6_comb_x2_x86_64.S` does all
twenty-two 768-bit operations in one body with one prologue.

| | lazy Fp6 in C | + asm multiply | + adc/sbb | **+ asm combination** |
|---|---:|---:|---:|---:|
| miller | +22.5% | +7.2% | +0.2% | **-6.3%** |
| pairing | +13.6% | +3.9% | +0.4% | **-3.4%** |
| final exponentiation | +10.5% | +3.0% | +0.6% | **-2.9%** |
| `gt_exp` | +43.0% | +12.8% | -1.4% | **-12.8%** |
| `gt_exp_ct` | | | | **-8.3%** |

**Eighteen fused multiply-and-reduce become eighteen raw products and six
reductions, one per output coefficient.** That is what blst does and what this
library did not.

Against blst, over three runs:

| | before | after |
|---|---|---|
| miller loop | 1.40x | **1.29x** |
| final exponentiation | 1.22x | **1.18x** |
| pairing, validation included | 1.36x | **1.28x** |

The result is bit-for-bit what it was: the pairing digest over 3,000 inputs is
`76da830b565dbafb` on the lazy path, on the eager path forced with
`ELIPS_NO_ASM`, and on the commit before any of this. Nothing about the value
changed, only how it is computed.

**Three routines, and each one was necessary.**

| | what it removes |
|---|---|
| `fp2_mulx2_x86_64.S` | three prologues and three 96-byte round trips per Fp2 product |
| `adc`/`sbb` intrinsics | `__int128` borrow extraction, 5,000 instructions per `fp6_mul` |
| `fp6_comb_x2_x86_64.S` | twenty-two more prologues and round trips |

Written in plain C the same algorithm is **22.5% slower** than the eager tower.
Every one of those three steps was needed to get from there to -6.3%, and the
middle one is not assembly at all.

**Scope.** The lazy path needs BMI2 and ADX, the same requirement as the fused
multiply, so it reuses the choice `fp_mul` already makes at startup. The
8-limb curves, non-x86-64 targets, and any machine without ADX take the eager
path, which is unchanged and still covered by the full suite under
`ELIPS_NO_ASM`. That is a fork in the tower, and it is the price of this.

## The sparse line multiply, lazy

The lazy Fp6 left one dense-reduction site inside the miller loop, and it was
the biggest: `fp12_mul_sparse035` is 42% of the loop.

It spent **33 Montgomery reductions per call**. Nine in `t0 = f0*c0`, eighteen
in the six-product `t1 = f1*(0,c3,c5)` block, and six inside the `fp6_mul` that
builds `s` -- that one already took the lazy path, but it reduced its result
only for the caller to subtract it narrow. An Fp12 has twelve Fp coefficients,
so **twelve is the floor**, and holding all fifteen products at 768 bits pays
exactly that.

The function moved from `src/pairing/miller.c` to `src/arith/fpx.c`, where the
wide machinery lives, and is declared in `elips/fpx.h`. The eager body is
unchanged and is still what everything without ADX runs.

### The combination in C, again

Written first with the `adc`/`sbb` intrinsics, the thirty wide operations the
combination needs gave the miller loop only **-2.2%** -- 21 reductions saved
per call and most of the saving spent getting 768-bit values in and out of
thirty function calls. Exactly the shape the Fp6 combination had.

### The combination in assembly

`src/arith/fp12_sparse_comb_x2_x86_64.S` does all thirty in one body with one
prologue: the `t1` combination, `f0 = t0 + v*t1`, and `f1 = s - t0 - t1`,
leaving twelve wide values for twelve `redcw`.

| | C combination | **assembly combination** |
|---|---:|---:|
| miller | -2.2% | **-5.2%, -5.6%** |
| pairing | -1.1% | **-2.0%, -2.4%** |
| `bls_verify` | -0.8% | **-1.6%, -2.4%** |

Two `--ab` runs, both at 100% agreement on those rows.

The wide add and subtract are now in `src/arith/wide_x86_64.h`, included by
both combination routines. Checked: the instruction stream `fp6_comb_x2` emits
is byte-identical to what it emitted with its own copy of the macros.

**`gt_exp` moved -2.9% and -3.2%, and this change is not on its path.** Nothing
in a GT exponentiation calls the sparse multiply. What moved it is code layout
in `fpx.o` from the function being added to that translation unit. It is real
and it reproduces, but it is not this change working, and it would be dishonest
to bank it.

### Against blst

| | before the lazy tower | after Fp6 | **after the line multiply** |
|---|---:|---:|---:|
| miller loop | 1.40x | 1.29x | **1.23x** |
| final exponentiation | 1.22x | 1.18x | 1.18x |
| pairing, validation included | 1.36x | 1.28x | **1.26x** |

### Verification

This is an exact algebraic identity, so the output must not move, and it does
not. The pairing digest over 3,000 inputs is identical on the lazy path, on the
eager path forced with `ELIPS_NO_ASM`, and on the commit before the change, on
**all three curves**: `b83e6405b41a18f2`, `91ac62b9ac181589`, `fe33dae4426b9a0a`.

The digest was shown able to fail: mis-indexing one of the thirty operands in
the lazy path moved it to `abe544634aa1e2ef` while the eager path stayed
correct, which also confirms `ELIPS_NO_ASM` really does take the eager branch.
`bench/fp12_sparse_comb_x2_test.c` checks the assembly against the C
combination over 200,000 random wide inputs and was likewise shown to fail on a
one-operand change. 59 tests pass on both paths, `ct_branch_scan` is clean over
51 objects, and both py_ecc cross-checks pass: the exact 12-coefficient pairing
comparison and 320 BLS signature checks.

## Where the gap actually was, measured

After the lazy tower the obvious assumption was that the remaining 1.23x on the
miller loop was per-multiply speed. It was not. Head to head on the same
machine, same harness, with the outputs checked bit-identical:

| ns | ELiPS | blst |
|---|---:|---:|
| `fp_mul`, fused Montgomery | **34.6** | 36.1 |
| raw 384x384 | 16.7 | **15.8** |
| `fp2_sqr` | 82.8 | 83.7 |
| `fp2_mul` | 124.6 | **119.0** |

**The kernels are at parity.** And so are the reductions: ~3,294 per miller
loop against blst's ~3,336. That lever is spent.

What is not at parity is the multiplication count. Callgrind on one miller loop
each: **ELiPS 7,602 Fp multiplications, blst 6,867.** Call counts are
backend-independent, which is what makes this comparison legal under valgrind.

Most of the 735 sat in one routine.

## Rotating the line by w: 15 Fp2 multiplications become 13

Our line is non-zero at w^0, w^3 and w^5. Storage puts w^0, w^2, w^4 in the
first Fp6 half and w^1, w^3, w^5 in the second, so the line splits **1+2** --
`L0 = (c0,0,0)`, `L1 = (0,c3,c5)` -- and Karatsuba over Fp12 costs
3 + 6 + 6 = 15. blst's line splits 2+1 and costs 5 + 3 + 5 = 13.

Scaling the line by one power of w moves it to (1,4,0), which in our storage is

    L0 = (xi c5, 0, c3),   L1 = (c0, 0, 0)

a **2+1** split. The scaling is a relabelling plus one `fp2_mul_xi`.

The two-non-zero Fp6 multiply is five products, not six. For `f = (a0,a1,a2)`
and `L = A + B v^2` with `v^3 = xi`:

    r0 = a0A + xi(a1B)
    r1 = a1A + xi(a2B)
    r2 = a2A + a0B  =  (a0+a2)(A+B) - a0A - a2B

so `a0A, a1B, a1A, a2B, (a0+a2)(A+B)` are enough, and all five are independent
of one another.

### The probe came first, and it should have

This project had already refuted a multiplication-count win in this exact
routine: Karatsuba on the six-product middle term removed 207 `fp_mul` per
miller loop and bought **nothing**, because those six products were independent
and filled issue slots the Karatsuba chain then serialised.

So before writing anything, two of the fifteen wide products were deleted from
the shipped routine, aliasing their outputs so the answer was wrong and the
shape was right, and the miller loop was measured against master, ABBA, 16
pairs:

    shipped  369.0 us      two products fewer  360.0 us      -2.4%

Every pair agreed. 138 fewer Fp2 products saved 9.0 us against a theoretical
11.4, which is coherent, so the machine is partly issue-bound but not entirely.
Worth building. The probe cost one edit.

### The fixup that turned out not to exist

The first version tracked the accumulated power of w through the loop and
divided it out at the end. It was written, it was correct, and then sabotaging
it changed no answer.

**Every power of w dies in the final exponentiation**, not just the Fp2 factors
the file header already relied on. `w^6 = xi` lies in Fp2, so the order of w
divides `6(p^2-1)`; r divides `p^12-1` and no smaller `p^k-1` because the
embedding degree is 12, and r is far larger than 6; so r divides neither
factor, r cannot divide the order of w, and `w^((p^12-1)/r) = 1`.

The exponent tracking, the per-curve fixup and `fp12_mul_w_pow` were all
deleted. The rotation is free. `pairing_test.c` now checks the property
directly on every vector, scaling a miller value by `w^1` through `w^5` and
confirming the final exponentiation is unmoved; sabotaged with `1 + w`, which
is not a power of w, it fails on all four vectors.

### Result

`src/arith/fp12_line_comb_x2_x86_64.S` does the new combination: forty-four
wide operations under one prologue, against the fifty-two the unrotated form
needed across two routines.

| | change | agreement |
|---|---:|---:|
| miller | **-3.8%, -3.4%** | 100%, 100% |
| pairing | **-1.5%, -1.8%** | 88%, 100% |
| `bls_verify` | -0.9%, -1.2% | 100%, 88% |

Fp multiplications per miller loop: **7,602 to 7,188**, exactly the 414
predicted (69 lines x 2 Fp2 products x 3). The gap against blst's 6,867 is now
321.

Against blst:

| | after Fp6 | after the lazy line | **after the rotation** |
|---|---:|---:|---:|
| miller loop | 1.29x | 1.23x | **1.18x** |
| final exponentiation | 1.18x | 1.18x | 1.18x |
| pairing, validation included | 1.28x | 1.26x | **1.24x** |

Over the two commits, the miller loop is **8.4% faster** and the pairing 3.1%.
`gt_exp` is where it was, which is correct: nothing in a GT exponentiation
touches the line multiply, and the -2.9% the previous commit showed on it was
code layout, declined at the time and duly given back.

### Verification

The pairing output is unchanged on all three curves, on both the lazy and the
`ELIPS_NO_ASM` eager paths: `b83e6405b41a18f2`, `91ac62b9ac181589`,
`fe33dae4426b9a0a`. 59 tests pass on both paths.
`bench/fp12_line_comb_x2_test.c` checks the assembly against the C combination
over 200,000 random wide inputs, and was shown to fail on a one-operand change.
`ct_branch_scan` is clean over 54 objects. Both py_ecc cross-checks pass: the
exact 12-coefficient pairing comparison and 320 BLS signature checks.

The rotation exponent differs per curve -- `w^3` on BLS12-381, `w^5` on
BLS12-461, `w^1` on BN-462 -- so if the freedom had been narrower than claimed,
the KAT vectors would have caught it on at least one of them.

## fp12_sqr_cyc, lazy: probed and refuted

The final exponentiation is the largest remaining gap at 1.18x, and
`fp12_sqr_cyc` is 87% of it, about 321 calls each. Callgrind gives its exact
shape: **18 `fp_mul` and 98 `fp_add`/`fp_sub`** per call. Twelve reductions is
the floor for a twelve-coefficient output, so a lazy version saves 6 of 18,
where the Fp6 work saved 12 of 18.

### The kernel costs, which is where this was decided

| ns | |
|---|---:|
| `fp_mul`, product and reduction fused | 34.55 |
| raw 384x384 product alone | 16.67 |
| `redcw`, reduction alone | 23.08 |
| `fp_add` / `fp_sub` | 3.95 |
| one wide modular add/sub | ~7.0 |

The wide figure comes from the two shipped combination routines, whose
operation counts are known exactly: `fp6_comb_x2` is 22 operations at 141.1 ns
and `fp12_line_comb_x2` is 44 at 308.5 ns, a slope of 7.6 ns and an average of
6.4 to 7.0.

**Read the first three rows again: 16.67 + 23.08 = 39.75 against 34.55.
Splitting a multiply into a product and a separate reduction costs 15% MORE
than the fused form.** Lazy reduction is never free; it only ever wins by
amortisation, and the amortisation has to be steep.

For `fp12_sqr_cyc` it is not:

| | now | lazy |
|---|---:|---:|
| 18 fused multiplies | 621.9 ns | |
| 18 raw + 12 reductions | | 577.1 ns |
| 98 narrow add/sub | 387.1 ns | |
| ~71 wide + ~18 narrow | | 568.1 ns |
| **total** | **1,009 ns** | **1,145 ns** |

A predicted **13% regression**, and the multiply side wins only 44.8 ns.

### The probe, run anyway

Three of the nine `fp2_sqr` calls inside `fp12_sqr_cyc` were replaced with two
raw products and no reduction, and six with two raw products and two `redcw`:
eighteen products and twelve reductions, exactly the lazy budget, with the
answer wrong and everything else unchanged. ABBA, ten rounds:

| | shipped | 18 raw + 12 redc | change | agreement |
|---|---:|---:|---:|---:|
| `final_exp` | 478.7 | 476.5 | -0.4% | 40% |
| `gt_exp` | 860.2 | 851.8 | -1.0% | 80% |
| `gt_exp_ct` | 331.4 | 329.9 | -0.5% | 70% |

The model predicted -3.0% on `final_exp` from the multiply side alone.
**Only about a seventh of it appeared, at agreement no better than a coin.**
The standalone kernel timings are throughput numbers from tight loops; inside
the real routine the multiplies interleave with dependent additions and the
saving hides. Same lesson as the refuted Karatsuba on the sparse middle term.

So the measured win is under 1% and the modelled cost is over 12%. Refuted.

A narrower variant was checked on paper before giving up: keep the outer
`3a +- 2b` combination narrow and make only the three `fp4_sqr` lazy. That
still reaches twelve reductions, because three `fp4_sqr` produce exactly twelve
coefficients, and it needs only 33 wide operations instead of 71. It comes to
1,041 ns against 1,009. Also a regression, for a new wide Fp2 squaring kernel
and a combination routine.

### The rule this gives the project

Every lazy-reduction change so far, scored by reductions saved per wide
operation added:

| | saved / wide | outcome |
|---|---:|---|
| `fp6_mul` | 12 / 22 = 0.55 | -6.3% miller |
| sparse line multiply | 21 / 52 = 0.40 | -5.4% miller |
| `fp12_sqr_cyc` | 6 / 71 = 0.085 | refuted |

Below roughly 0.3, do not write it. Above it, probe first anyway.

## ep2_in_subgroup, traced

The subgroup-check section above ends by asserting the remaining gap is
"formula cost: ELiPS doubles and adds with the complete formulas throughout,
blst uses dedicated ones". That was reasoning. Here is the trace, and it turned
out the assertion was right about the residual and wrong about all of it.

### The decomposition

| | us | share |
|---|---:|---:|
| `ep2_in_subgroup`, whole | 92.39 | |
| `ep2_mul_pubconst` | 89.50 | **96.9%** |
| `ep2_on_curve` | 0.93 | 1.0% |
| `ep2_psi` | 0.25 | 0.3% |
| `ep2_eq` | 0.25 | 0.3% |
| one `ep2_dbl` | 1.255 | |
| one `ep2_add` | 1.910 | |

64 doublings and 5 additions is 64 x 1.255 + 5 x 1.910 = 89.9, which is the
whole ladder. **Nothing outside the doubling matters**, and the doubling alone
is 89% of the check.

### It is not the chain, and it is not the test

blst's `POINTonE2_in_G2` is Scott eprint 2021/1130, `psi(P) == [z]P`, and so is
ours. blst's `POINTonE2_times_minus_z` is a hardcoded chain: one double, then
`add_n_dbl` by 2, 3, 9, 32, 16. That is **63 doublings and 5 additions**. Our
NAF of the BLS12-381 seed has weight 6, so 63 doublings and 5 additions.
Identical.

Fp multiplications per check, callgrind, both libraries:

| | ELiPS | blst |
|---|---:|---:|
| before | 2,047 | 1,281 |

One difference worth stating because it makes the row not quite like for like:
`blst_p2_affine_in_g2` does **no on-curve check**, and `ep2_in_subgroup` does.
That is 1.0% of ours, so it is not the cause, but ours is doing strictly more.

### What the trace actually found

`PT(dbl)` opens with

```c
F(mul)(t0, p->y, p->y);     /* Y^2  */
F(mul)(t2, p->z, p->z);     /* Z^2  */
```

**Two of its nine products are squarings written as multiplies**, and
`PT(on_curve)` has three more. `fp2_sqr` is two Fp products against
`fp2_mul`'s three: 82.5 ns against 132.6.

The formula is untouched and so is the exception-freeness argument. These are
the same values, computed the cheaper way.

### One thing that had to be measured, not assumed

Routing the Fp instantiation through `F(sqr)` too made G1 **slower**: `ep_mul`
+1.3% at 100% agreement. `fp_sqr` IS `fp_mul(a, a)` -- the note on it in
`fp.c` says why a dedicated one was measured and not written -- so going
through it buys nothing and adds a call. The template now asks the
instantiation, via `EC_CHEAP_SQR`, which the Fp2 one sets and the Fp one does
not.

### Result

| | change | agreement |
|---|---:|---:|
| `ep2_in_subgroup` | **-6.6%, -7.3%** | 100%, 100% |
| `ep2_mul` | **-4.8%, -5.2%** | 100%, 100% |
| `ep2_mul_glv` | -2.0%, -2.3% | 75%, 88% |
| `bls_sign` | -1.5%, -2.0% | 88%, 100% |
| `hash_to_g2` | -1.0%, -3.5% | 62%, 88% |

Fp multiplications per check: **2,047 to 1,916**, which is exactly the 131
squarings each saving one. Against blst on the same host, `in_g2` goes
**1.73x to 1.61x**.

One run showed `ep_mul_glv` at +3.0% and 100% agreement. It is not real: all
fifteen G1 functions are instruction-identical before and after, and
`ep_mul_glv` differs only in branch target addresses because it moved in the
object file. A second run put it at +0.7% at 75%.

### The residual, and why it stays

1,916 against blst's 1,281. We double with the RCB complete formulas, now
7 Fp2 multiplications and 2 squarings; blst uses dedicated ones. That is the
remaining 1.61x, and it is a deliberate choice rather than an oversight: the
caller is validating a point it does not trust, so every exceptional case has
to work rather than be argued away. The Miller loop could give up the complete
formulas because T there is a known multiple of a prime-order point. Here there
is no such argument available.

## Revised ordering

| | worth | risk |
|---|---|---|
| `dbl_line` formula, 11M+3S → 3M+8S | ~8% of the miller loop | contained: one function, existing KATs, complete-vs-dedicated needs the subgroup argument written down |
| 8M2 line kernel (`normalized_miller_ref.py`) | 17% of multiplications | normalization cost this call pattern cannot amortise; needs a projective preparation first |
| lazy reduction through Fp6/Fp12 | 7.7% measured at Fp2 | touches the whole tower; adds additive work to remove reductions |

The doubling step is now the first thing to do: it is the largest single
identified item, the most contained, and the only one that reduces
multiplications and additions together.

**The target is now 1.52x, not 2.2x, and parity is in reach.** The doubling
formula and lazy reduction together take 468 us to roughly 366, which is 1.20x
against blst's 304. The final exponentiation is already at 1.22x. Neither
number was reachable under the old table, and both follow from the same
levers. The 762 Mp recorded here as unattributed was an artifact:
`fp12_mul_sparse035` is 15 Fp2 multiplications per call, not the 9 counted.
The trace followed direct calls and not the six inside `fp6_mul`, which over
69 lines is 1,242 Fp multiplications.

## Reproducing any of this

```bash
git clone --depth 1 https://github.com/supranational/blst
(cd blst && CFLAGS="-O3 -fno-builtin -fPIC -D__ADX__" ./build.sh)
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG"
cmake --build build -j
gcc -O3 -std=c11 -Iinclude -Iblst/bindings -DELIPS_CURVE_BLS12_381 \
    bench/cross_blst.c -o cross \
    build/libelips_arith_BLS12_381.a blst/libblst.a -lm
./cross --reps 31                          # both libraries, one harness
```

Operation counts come from `valgrind --tool=callgrind` on a program that
performs exactly one miller loop, parsing `calls=` out of the output. Call
counts are backend-independent, which matters because valgrind hides ADX from
CPUID and both libraries fall back to their portable paths under it.
