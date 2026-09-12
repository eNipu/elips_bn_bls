#!/usr/bin/env python3
"""
Look in the OBJECT FILES for masked selects the compiler turned into branches.

Why this exists. A masked select is written branch-free in C:

    mask = 0 - condition;            /* 0 or ~0 */
    r    = (a & mask) | (b & ~mask);

and a compiler may notice that `mask` is provably one of two values, decide the
expression means "a or b", and emit a branch on the condition. The source still
reads as constant time. Nothing warns. Only the disassembly says otherwise.

This has happened twice in this repository, in different instructions:

    src/arith/wide.c   clang emitted `neg` then `jae`
    glv_divrem         clang emitted `setb` then `test $1` then `je`, a branch
                       on a borrow derived from the secret scalar (issue #30)

The second hid from a search for the first, because the same defect wears
different instructions depending on how the condition was computed. dudect
caught both, but only statistically and only when it happened to be run on the
affected toolchain -- issue #30 read 66 under clang and 2 under gcc and was
blamed on the platform for weeks, because macOS is where clang runs by default.

So this looks at the instructions instead. It is a structural check, not a
statistical one: it cannot be flaky and it does not need a quiet machine.

WHAT IT IS NOT. It does not prove constant-time execution. It finds one
specific, recurring, compiler-introduced defect. dudect still runs.
"""
import argparse
import glob
import os
import re
import subprocess
import sys

# Routines that are variable time ON PURPOSE. Each one needs a reason, and
# "it kept failing" is not one.
ALLOW = {
    # The dudect negative control: it must branch on its operand, or the
    # timing harness has nothing it is known to catch. See test/dudect_test.c.
    "fp_inv_vartime",
}

# The exact shape, and deliberately not a broad one.
#
#     set<cc>  %r8b            materialise a condition into a byte
#     test $1, %r8b            look at that same byte
#     j<cc>    ...             branch on it
#
# Three instructions, with the REGISTER MATCHING across them. That is what
# "the compiler decided my mask was really a condition" looks like, and it is
# what clang emitted for glv_divrem.
#
# A broader pattern was tried first and is useless: pairing any flag-producing
# instruction with any following conditional jump flags ten sites on both
# compilers, all of them fine -- the deserializers branching on public encoded
# bytes, and loop-unrolling remainder checks on public trip counts. A check
# with ten false positives is not a gate, it is something that gets switched
# off. So this matches the register too.
SETCC = re.compile(r"^set([a-z]+)\s+%([a-z0-9]+)$")
TESTB = re.compile(r"^test\s+\$0x1,%([a-z0-9]+)$")
COND_JUMP = re.compile(r"^j(e|ne|b|ae|be|a|s|ns)\b")

# THE SECOND SHAPE, and the reason there is a second one.
#
#     sub  %rbx,%r10          subtract, on 64-bit registers
#     js   ...                branch on the sign of the result
#
# No setcc in sight, so the pattern above walks straight past it. This is what
# clang emitted for redc_lin's conditional subtraction (issue #32): the mask
# was `0 - ((tmp[hi] >> 63) == 0)` and clang read the select as "take tmp or
# take acc" and branched on the borrow. Four of them, one per unrolled round,
# inside fp_inv. dudect had been reading 3 to 24 on fp_inv for weeks and it
# was very nearly written off as hardware operand latency.
#
# TWO RESTRICTIONS, both of which earn their keep. Without them this shape
# fires on six sites that are all fine:
#
#   64-bit registers only. Secrets here are limb_t, which is 64-bit. Lengths,
#   indices and trip counts are int, and clang does that arithmetic in %eXX --
#   which is exactly what elips_mod_wide's `sub %eax,%esi ; js` is, a check on
#   two public limb counts.
#
#   no immediates. Comparing against a constant is comparing against a public
#   value. fp12_exp_param's `cmp $0xe,%ebp ; js` walks a fixed, published
#   addition chain.
#
# Same rule as above: match narrowly, or it becomes a check somebody turns off.
R64 = r"%r(?:[abcd]x|si|di|bp|sp|8|9|1[0-5])"
SUBCMP = re.compile(rf"^(sub|sbb|cmp)\s+{R64},{R64}$")
SIGN_JUMP = re.compile(r"^j(s|ns)\b")


def scan(obj):
    out = subprocess.run(["objdump", "-d", obj], capture_output=True, text=True)
    if out.returncode != 0:
        return []
    hits, fn = [], None
    window = []                      # (mnemonic, text) for the last two
    for line in out.stdout.splitlines():
        m = re.match(r"^[0-9a-f]+ <(.+)>:", line)
        if m:
            fn, window = m.group(1), []
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ins = re.sub(r"\s+", " ", parts[2].strip())
        ins = re.sub(r"\s*<.*", "", ins)
        window.append(ins)
        if len(window) > 3:
            window.pop(0)
        if len(window) == 3:
            a, b, c = window
            ma, mb = SETCC.match(a), TESTB.match(b)
            if ma and mb and ma.group(2) == mb.group(1) and COND_JUMP.match(c):
                hits.append((fn, f"{a} ; {b} ; {c}"))
        # The second shape spans two instructions, not three, and clang puts
        # unrelated work between them, so look back over the whole window.
        if window and SIGN_JUMP.match(window[-1]):
            for prev in reversed(window[:-1]):
                if SUBCMP.match(prev):
                    hits.append((fn, f"{prev} ; ... ; {window[-1]}"))
                    break
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build", help="a configured build directory")
    args = ap.parse_args()

    pattern = os.path.join(args.build, "CMakeFiles", "elips_arith_*.dir", "src", "*", "*.o")
    objects = sorted(glob.glob(pattern))
    if not objects:
        # A scanner that finds nothing because it looked nowhere is the worst
        # possible outcome: a green check that means nothing.
        print(f"error: no object files under {pattern}", file=sys.stderr)
        return 2

    total, failures = 0, []
    for obj in objects:
        for fn, seq in scan(obj):
            total += 1
            if fn in ALLOW:
                continue
            failures.append((os.path.basename(obj), fn, seq))

    print(f"scanned {len(objects)} objects, {total} matches, "
          f"{len(failures)} outside the allow list")
    for obj, fn, seq in failures:
        print(f"  {obj}: {fn}: {seq}")

    if failures:
        print("\nA masked select may have become a branch. Check the source, and if it "
              "has, launder the mask through ct_mask (src/arith/ct.h).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
