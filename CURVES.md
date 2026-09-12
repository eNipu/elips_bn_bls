# Curves

| Curve | p | r | Generators | Serialized G1 / G2 | hash-to-curve |
|---|---|---|---|---|---|
| BLS12-381 | 381 bits | 255 bits | from the specification | 48 / 96 bytes | RFC 9380 `SSWU_RO_`, byte-exact |
| BLS12-461 | 461 bits | 308 bits | generated | 58 / 116 bytes | RFC 9380 `SVDW_RO_`, no registered suite |
| BN-462 | 462 bits | 462 bits | generated | 59 / 118 bytes | RFC 9380 `SVDW_RO_`, no registered suite |

Sizes are the compressed encodings; uncompressed is twice each. BN-462 needs one
byte more than its field width because the three flag bits do not fit otherwise.

Only BLS12-381 has a specification to conform to. For it, the generators, the
point encodings and the hash-to-curve outputs are pinned against the published
values and match any conforming implementation. The other two curves have no
standard, so their vectors pin self-consistency rather than interoperability.

**On BLS12 the pairing returns `e³`, not `e`.** That is a property of the
standard final-exponentiation chain, which RELIC and the original library also
use, and not a defect: since `gcd(3, r) = 1` it is still bilinear and
non-degenerate, so any protocol that only compares pairings is unaffected. Raw
values will not match an implementation that outputs `e`; use
`pairing_final_exp_plain` when the exact value is needed. On BN the value is
`e` exactly.

## Choosing one

The arithmetic is compiled for one curve: the field width has to be a
compile-time constant for the loops to unroll and for the assembly to target
a known width.

```bash
cmake -B build -DELIPS_CURVE=BLS12_381    # default, and the only standardised one
cmake -B build -DELIPS_CURVE=BLS12_461
cmake -B build -DELIPS_CURVE=BN_462
```

**Use BLS12-381 unless you have a reason not to.** It is the curve with a
specification, so it is the one whose generators, encodings and hash-to-curve
outputs match other implementations byte for byte. The signature layer is
BLS12-381 only, because the other two have no registered ciphersuite and an
"IETF BLS signature on BN-462" would interoperate with nothing.
