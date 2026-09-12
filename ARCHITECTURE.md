# Architecture

How the pieces fit, and why the boundaries are where they are. Every edge
below was read out of the includes and the build files rather than recalled;
the layering is enforced by what each header is allowed to include.

## The shape of it

```
                     ┌──────────────────────────────────────────┐
   consumers         │  Python (cffi)    JavaScript (wasm)   C  │
                     └────────┬──────────────┬──────────────┬───┘
                              │              │              │
                              └──────────────┴──────────────┘
                                             │
                                    the byte-oriented API
                                             │
  ┌──────────────────────────────────────────▼──────────────────────────────┐
  │  elips/bls.h            sign · verify · aggregate · proof of possession  │
  │  src/bls/               draft-irtf-cfrg-bls-signature, BLS12-381 only    │
  │                         hkdf.c  RFC 5869, for key generation            │
  └─────┬───────────┬───────────────┬──────────────┬────────────┬───────────┘
        │           │               │              │            │
  ┌─────▼─────┐ ┌───▼──────────┐ ┌──▼──────────┐ ┌─▼────────┐ ┌─▼─────────┐
  │ pairing.h │ │hash_to_curve │ │ serialize.h │ │ random.h │ │ sha256.h  │
  │ miller.c  │ │ RFC 9380     │ │ compressed  │ │ + sysrand│ │ RFC 6234  │
  │ optimal   │ │ SSWU / SVDW  │ │ encodings,  │ │ CSPRNG   │ │           │
  │ ate, final│ │ + isogenies  │ │ validation  │ │          │ │           │
  │ exp       │ │              │ │             │ │          │ │           │
  └─────┬─────┘ └───┬──────────┘ └──┬──────────┘ └─┬────────┘ └───────────┘
        │           │               │              │
        └───────────┴───────┬───────┴──────────────┘
                            │
                    ┌───────▼────────┐
                    │  elips/ec.h    │   G1 over Fp, G2 over Fp2
                    │  src/arith/    │   ladders, GLV, subgroup tests
                    │  ec.c          │   ec_tmpl.h instantiated twice
                    └───────┬────────┘
                            │
                    ┌───────▼────────┐
                    │  elips/fpx.h   │   Fp2 → Fp6 → Fp12
                    │  fpx.c         │   Karabina compressed squaring
                    └───────┬────────┘
                            │
                    ┌───────▼────────┐
                    │  elips/fp.h    │   Fp, Montgomery form
                    │  fp.c          │   divstep inversion
                    │  wide.c        │   wide reduction (replaced GMP)
                    │  ct.h          │   the optimisation barrier
                    └───────┬────────┘
                            │
                 ┌──────────▼───────────┐
                 │  fp_mul_x86_64.S     │  mulx/adcx/adox   (CPUID-gated)
                 │  fp_addsub_x86_64.S  │  adc/sbb/cmov     (always on)
                 │  fp_mul_aarch64.S    │  mul/umulh/adcs   (baseline)
                 └──────────────────────┘
```

Arrows point from a module to what it depends on. There are no cycles and no
upward edges: `fp` knows nothing about curves, `ec` knows nothing about
pairings, and nothing below `bls` knows a signature scheme exists.

## Why the layers are where they are

**One curve per build.** `ELIPS_CURVE` selects BLS12-381, BLS12-461 or BN-462
at configure time, and the field width becomes a compile-time constant. That is
not a packaging convenience: the loops unroll and the assembly targets a known
limb count only because `FP_LIMBS` is fixed. The test suite builds all three,
because vectors have to cover all three. `CURVES.md` has the details.

**`elips_bls` is a separate target**, not a layer inside the core. It consumes
only the public headers, so a project that wants the pairing and not the
protocol links `ELiPS::arith` and never compiles it. This is the line
`MODERNIZATION_PLAN.md` §10.8 drew and it still holds: hash-to-curve is a curve
operation with a curve specification and belongs below; ciphersuites and
proof-of-possession do not.

**The signature layer is BLS12-381 only.** The other two curves have no
registered ciphersuite, and an "IETF BLS signature on BN-462" would
interoperate with nothing.

**GMP is a test dependency, not a runtime one.** Nothing under `src/` includes
a third-party header and the installed library links nothing but libc. GMP is
kept deliberately as an *independent oracle*: the tests check the Montgomery
layer, the wide reduction and both inversions against it, because agreement
with a separate implementation is evidence in a way agreement with ourselves is
not. `wide.c` exists because that reduction used to be `mpn_sec_div_r`.

## The two things that cross every layer

These are the reason several boundaries look the way they do.

**Constant time.** Every routine that touches a secret key is branch-free, and
the masks that make it so go through `ct_mask` in `src/arith/ct.h`. That helper
exists because a masked select is a value the compiler can prove is one of two
things, and both gcc and clang have turned one back into a branch here — in
`wide.c`, then in `glv_divrem` (#30), then in `redc_lin` (#32), wearing
different instructions each time. Two independent checks run on every build:
`dudect` measures timing statistically, and `tools/verify/ct_branch_scan.py`
reads the object code for the instruction shapes themselves. They fail
differently on purpose.

**Secret material never reaches a branch, and never reaches the heap.** Keys
live in caller-provided buffers; nothing below `bls.c` allocates.

## How a signature actually flows

`elips_bls_sign(sig, sk, msg, len)`:

1. `bls.c` deserializes `sk` to a scalar and rejects it if out of range.
2. `hash_to_curve.c` maps `msg` to a point in G2 — SHA-256 through
   `expand_message_xmd`, then SSWU with the isogeny, then cofactor clearing.
3. `ec.c` multiplies that point by the scalar with `ep2_mul`, the plain
   constant-time fixed-window ladder — deliberately **not** `ep2_mul_glv`,
   which is about 1.5x faster. The note at the top of `bls.c` records why: when
   it was written the GLV paths read `max|t| = 32.97 -> 69.36` in dudect under
   clang against 0.94 for the plain ladder, and signing is the one place in
   this library where a scalar is a long-term secret. That was #30, and it is
   fixed; see the caveat below.
4. `serialize.c` compresses the result to 96 bytes.

`elips_bls_verify` validates the public key and then makes one multi-pairing
check. `pairing_check` negates the generator into slot 0 and evaluates
`e(-g1, sig) · e(pk, H(m)) == 1` through `elips_pairing_multi`, so a single
final exponentiation covers both sides. It fills slot 0 itself rather than
copying the caller's arrays, because at aggregate width that copy is about
55 KB and WebAssembly's default stack is 64 KB.

> **Open, and the code asks for it.** The `bls.c` note ends "when #30 closes,
> this decision should be revisited with a benchmark, not assumed." #30 and #32
> are both closed, and the GLV paths now read 1.50 to 2.62 under clang across
> three runs. Signing therefore pays about 245 µs of 1510 µs for a constraint
> that no longer holds. Tracked, not silently changed: switching which routine
> handles a long-term secret key is not a documentation edit.

## The three consumers

All three compile the same C. None reimplements anything.

| | how it binds | what it ships |
|---|---|---|
| **C** | `find_package(ELiPS)`, link `ELiPS::bls` or `ELiPS::arith` | headers + static libs |
| **Python** | cffi, API mode, out-of-line — `bindings/python/build_ffi.py` compiles the same sources into one extension | a wheel with no system dependencies |
| **JavaScript** | emscripten — `bindings/js/build.sh`, twelve `elips_bls_*` entry points | one `.mjs` with the wasm embedded |

The Python and wasm builds keep their own source lists, which is a real hazard:
adding `fp_addsub_x86_64.S` to CMake and not to the Python list produced a wheel
that imported and died on an undefined symbol. `tools/verify/asm_sources_listed.py`
now asserts the lists agree. The wasm build is exempt on purpose — it targets
wasm32, where every `.S` is guarded out and `fp.c` takes the portable path.

**The portable C is never a second-class path.** Every assembly routine ships
alongside the C it replaces, both run against the same vectors on every build,
and `ELIPS_NO_ASM=1` switches the whole lot off at process start. That is not
only a fallback for other architectures; it is what makes the assembly testable
at all, because the differential test needs something to differ from.

## Where correctness comes from

Nothing in `test/` compares one C implementation against another.

- `tools/reference/` is an independent Python implementation of the field
  tower, the curve arithmetic and the optimal ate pairing, written from the
  defining equations rather than from `src/`. It generates `test/kat/`.
- RFC 9380's published vectors pin hash-to-curve byte for byte.
- `py_ecc` cross-checks the signature layer on every build, in both directions,
  because an implementation can be self-consistently wrong and pass a one-way
  check.
- The sanitizers, the timing tests and the codegen scan cover what vectors
  cannot.

`CONTRIBUTING.md` has the commands.
