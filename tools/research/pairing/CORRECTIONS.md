# Corrections to the earlier research notes

The archived notes are retained as requested, not silently rewritten.
Their conclusions are historical hypotheses, not an optimality certificate.

1. **No hardware or mathematical floor was proved.** The existing 36Mp
   general squaring is one implementation. A primitive cost estimate did not
   benchmark a complete competing squaring. The new scaled-square experiments
   implement and test counterexamples to the claim that this implementation
   cannot be improved.
2. **LLL plus Babai is heuristic, not exact closest-vector enumeration.**
   Recovering the known degree-five representation is a positive control. A
   large-coefficient degree-four result does not rule out a smaller one, and
   says still less about all exponentiation circuits. Degree in x is not in
   one-to-one correspondence with the number of exponentiations in a DAG.
3. **Subfield scaling is established prior art, not a new theorem here.**
   Devegili et al., [ePrint 2006/471](https://eprint.iacr.org/2006/471),
   section 2, explicitly uses final exponentiation to discard nonzero
   proper-subfield factors and clear Toom-Cook interpolation denominators.
   The divisibility checks remain valid; the novelty label does not.
4. **Partial easy exponentiation is not enough for Granger–Scott squaring.**
   Raising to p^6-1 gives the quadratic norm-one subgroup, not in general the
   Phi12 cyclotomic subgroup. The additional p^2+1 step is necessary.
   The old illustrative operation comparison is not a complete implementation
   of an early-projection algorithm.
5. **The remaining optimization space was not exhausted.** Dense support alone
   is not a lower bound on multiplication cost. “No shared products exist”
   was not proved. Seed Hamming weights do not exclude every mixed-radix
   chain. Claims of universal impossibility or optimality are retracted.
6. **The norm experiment is evidence, not a cheap runtime check.** It tests
   one public P/Q pair at m=2,3,5 on each of two curves. It retains vertical
   denominators and uses costly full-field operations. No compressed
   recurrence, general exceptional-case handling, or runtime saving was
   implemented.
7. **Norm notation in the old memo was inconsistent.** The common Fp2 norm
   of an Fp4 value B is N_(4/2)(B), computed by `norm6(B)` in the Fp12
   embedding, not N_(6/2)(B).
8. **The 64 squarings are specific to the current signed-digit schedule.**
   BLS12-381's seed has a 64-bit absolute binary expansion, which can use 63
   doubling steps. Neither schedule's operation count proves a universal
   lower bound for pairing computation.
9. **The two Karabina topics should not be conflated.** Torus compression
   (ePrint 2010/525) is related background, not the same construction as the
   compressed cyclotomic squaring used in the production final exponentiation.

The valid earlier outcomes are the explicit identities, the tests actually
run, the bounded searches' actual outputs, and the workload-specific timing
measurements. Publication-level novelty and optimality remain unestablished.
