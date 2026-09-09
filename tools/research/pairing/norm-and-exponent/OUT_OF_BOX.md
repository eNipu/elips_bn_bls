# Out-of-box pairing optimization sweep (BLS12-381 / BN)

Research note, local only. Baseline: branch `research/normalized-miller-prototype`
(8M2 normalized prepared-line kernel committed, measured −19..−31% Miller loop).

## Bottleneck inventory (prepared path, per Miller iteration, Mp = Fp mult)

- General Fp12 squaring: 36 Mp (2 Karatsuba Fp6 muls; fp12_sqr in src/arith/fpx.c)
- Scaled sparse line multiply (committed kernel): 8 M2 + 2 Fp2·Fp scalings = 28 Mp
- Per iteration: 64 Mp; loop total ≈ 64·36 + 69·28 = 4236 Mp
- Squaring is now the dominant term (2304 Mp, 54% of the loop) and is sequential
  (64 doublings forced by the BLS12 seed; already minimal loop length for optimal ate).

## Verified new invariant (research lead, not a speedup)

Prescribed norms of the complete (denominator-retaining) Miller function:
for scalar m with ladder point T = [m]Q,

- N_{12/6}(f) = (xP − xQ)^m / (xP − xT)   in Fp6
- N_{12/4}(f) = (yP − yQ)^m / (yP − yT)   in Fp4
- Fp2 compatibility: N_{6/2}(N_{12/4}(f)) = (xP³ − xQ³)^m / (xP³ − xT³)

Verified numerically on BLS12-381 and BN-462, m ∈ {2,3,5} (norm_state.py, all PASS).
Also verified: the inclusion-exclusion projector g = f⁶·C/(A³·B²) equals
f^((p⁶−1)(p²+1)(p²−2)) and lands in the Φ12 cyclotomic subgroup.

Potential (unproven): track (compressed f, T) as Miller state; norms update by
multiplying known per-line norms. Blockers: norm updates per step may cost as much
as the squaring saved; recovery needs decompression inversions (Karabina-like).

## Prior art boundary

- Naehrig–Barreto–Schwabe, "On compressible pairings and their computation" (SAC 2008):
  raise each line to q³−1 inside the loop, work in T2(Fq3), compress to 2 Fq
  elements using the norm relation −3b1b2ξ + ξ + 3b0² = 0 (their Prop. 1). This is
  the same relation family as the independently derived 3ab − 3c²ξ = 1. Covers
  generalized Eta/Ate pairings on 6|k curves. => the "compressed Miller state via
  norm structure" direction is published; any novelty claim must be narrow.
- Karabina ePrint 2010/525: torus compression factor 4/6 (already used in this repo's
  final exponentiation).
- Mori et al. pseudo-sparse multiplication; Costello–Stebila fixed-argument pairings:
  adjacent to the committed 8M2 kernel; novelty of that kernel still unresolved.

## Quantitatively rejected candidates (op counts verified)

1. Early easy-part injection (line → line^(p⁶−1) per step, loop in cyclotomic
   subgroup with Granger–Scott squarings): 72 vs 64 Mp/iter; loop 4878 vs 4236 Mp.
   Lines become dense; loss ≈ 15%. (NBS variant works only with their torus product
   trick and their curve/orientation conditions.)
2. Window-2 line merging for fixed Q: support grows 3 → 5 (L0·L1) / 6 (L0²·L1),
   verified numerically; merged multiply ≈ 18 M2 vs 2×8 = 16 M2. Wash/loss.
3. Fused square-then-line trilinear kernels: no shared products; ≈ 123 vs 64 Mp.
4. Shorter BLS12 final-exponent relations: bounded symbolic nullspace search
   (x-degree ≤ 4, p-degree ≤ 3, λ·Φ12/3 multiplier) found zero relations.

## Round 2: scaling freedom, measured ratios, final-exp audit

**Verified theorem (new):** (p^6−1) | (p^12−1)/r on BLS12-381, BLS12-461 and
BN-462, so the Miller accumulator may be scaled by arbitrary Fp6* factors with
zero effect on the pairing value. The committed kernel uses only Fp2* freedom.
(Stronger than needed for any kernel found so far; recorded as design space.)

**Measured primitive ratios (M2 Max, this library):** fp_mul ~42 ns,
fp_add ~7-11 ns (M/A ~ 4-6), fp2_mul ~3.4-3.6 Mp, fp2_sqr ~2.1-2.2 Mp.

**Rejected with measurements:**
- Fp4-cubic tower + scaled Toom-3 squaring (constants absorbed into the free
  Fp6* scalar): 5 S4 = 30 Mp but ~50 interpolation Fp-adds at M/A ~5 cost
  ~10 Mp => ~42 Mp > 36 Mp floor. Line side is a wash (3 M4 + scalings ~29 Mp
  vs 28 Mp). Dead on this hardware.
- Cyclotomic cubing chains (ePrint 2025/958, 2026/885): BLS12-381 seed x has
  ternary weight 28 (dense); binary weight is 6 (sparse). Cubing wins only on
  ternary-sparse seeds (SG54, BLS15/27). Not applicable.
- ePrint 2025/1387 multi-exp chain: 10% claim is for a different BLS12 curve
  (44-bit seed t = -2^44+2^17+1). On BLS12-381's 64-bit HW-6 seed the repo's
  nested chain (5 exp-by-x, ~11 muls) already beats their generic multiexp
  (5 exps, ~25-40 muls).

**Final-exp depth search (validated LLL + Babai, exact arithmetic):**
- Sanity: depth-5 basis (x^5 p^3) recovers the known max|coeff|=3
  representation == the current chain. Method trustworthy.
- Depth-4 bases (x^4 p^3, x^4 p^4): nearest representations need coefficients
  ~1.5e19 (the size of x itself) => each coefficient would cost a full
  64-bit exponentiation. No useful 4-exp-by-x chain exists in this family.
- Conclusion: the repo's final exponentiation is at the known optimum; the
  bounded search rules out the one remaining structural improvement.

## Net

- Strongest result remains the committed 8M2 kernel (+ separate tangent reduction,
  11M2+3S2 → 6M2+3S2, already tested).
- The remaining bottleneck (general Fp12 squaring in the loop) sits at the known
  lower bound for this tower; the only routes around it (compressed states) are
  published and lose on BLS12-381's tower.
- The prescribed-norm invariant is exact, verified, and may be worth keeping as a
  runtime self-check / test oracle; turning it into a speed mechanism is an open
  problem with published neighbors.
- After round 2, every remaining avenue was closed with either a measurement
  (Toom adds), a validated lattice search (no 4-exp chain), or a parameter
  property (ternary-dense seed). The practical frontier is engineering, not
  algebra: lazy reduction in fp2/fp6 arithmetic and vectorization (cf. ePrint
  2025/1283) are the known-open wins for this codebase.
