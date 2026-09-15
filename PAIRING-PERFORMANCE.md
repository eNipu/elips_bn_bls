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
