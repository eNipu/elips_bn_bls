#!/usr/bin/env python3
"""
elips_bls against py_ecc, both directions, on randomized inputs.

py_ecc is the Ethereum Foundation's Python implementation of
draft-irtf-cfrg-bls-signature. Nobody here wrote it, so agreement across a
large sample is evidence rather than a restatement of our own code. It is the
same argument that puts GMP in the test suite and PARI in tools/verify.

BOTH DIRECTIONS matter and they catch different things. A signature we produce
and py_ecc accepts proves our signing follows the draft. A signature py_ecc
produces and we accept proves our verification does. An implementation can be
self-consistently wrong and pass a one-directional check.

The rejection cases are not decoration. A verifier that cannot reject is not a
verifier, which is why this repository already carries kat.detects_corruption
and h2c.detects_corruption. Every rejection asserted here fails if the
corresponding check is removed from src/bls/bls.c.

Usage:  crosscheck_bls_pyecc.py <path-to-bls_cli> [iterations]
"""
import os
import subprocess
import sys

# The proof-of-possession scheme, which is the one elips_bls implements.
from py_ecc.bls import G2ProofOfPossession as BLS

CLI = None
failures = 0
checks = 0


def run(*args):
    """One bls_cli invocation. Returns (ok, output)."""
    r = subprocess.run([CLI, *args], capture_output=True, text=True)
    out = r.stdout.strip()
    return r.returncode == 0, out


def check(cond, what):
    global failures, checks
    checks += 1
    if not cond:
        failures += 1
        print(f"  FAIL {what}")


def expect_ok(what, *args):
    ok, out = run(*args)
    check(ok, what)
    return out


def expect_err(what, *args):
    ok, out = run(*args)
    check(not ok, f"{what} (should be rejected, got {out!r})")


def main():
    global CLI
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    CLI = sys.argv[1]
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 8

    print("elips_bls against py_ecc (draft-irtf-cfrg-bls-signature)")

    for it in range(iters):
        ikm = os.urandom(32)
        msg = os.urandom(1 + (it * 7) % 64)

        # --- KeyGen, SkToPk, Sign must agree BYTE FOR BYTE ----------------
        sk_c = expect_ok("keygen", "keygen", ikm.hex())
        sk_p = "%064x" % BLS.KeyGen(ikm)
        check(sk_c == sk_p, f"KeyGen differs\n    c={sk_c}\n    p={sk_p}")

        pk_c = expect_ok("sk_to_pk", "sk_to_pk", sk_c)
        pk_p = BLS.SkToPk(int(sk_c, 16)).hex()
        check(pk_c == pk_p, f"SkToPk differs\n    c={pk_c}\n    p={pk_p}")

        sig_c = expect_ok("sign", "sign", sk_c, msg.hex())
        sig_p = BLS.Sign(int(sk_c, 16), msg).hex()
        check(sig_c == sig_p, f"Sign differs\n    c={sig_c}\n    p={sig_p}")

        # --- cross verification, both directions --------------------------
        expect_ok("we verify our own", "verify", pk_c, msg.hex(), sig_c)
        check(BLS.Verify(bytes.fromhex(pk_c), msg, bytes.fromhex(sig_c)),
              "py_ecc rejects a signature we produced")
        expect_ok("we verify py_ecc's", "verify", pk_c, msg.hex(), sig_p)

        # --- rejection: every one of these must fail -----------------------
        bad_msg = bytes([msg[0] ^ 1]) + msg[1:] if msg else b"x"
        expect_err("tampered message", "verify", pk_c, bad_msg.hex(), sig_c)

        other = "%064x" % BLS.KeyGen(os.urandom(32))
        other_pk = BLS.SkToPk(int(other, 16)).hex()
        expect_err("wrong public key", "verify", other_pk, msg.hex(), sig_c)

        flipped = bytearray(bytes.fromhex(sig_c))
        flipped[-1] ^= 1
        expect_err("mangled signature", "verify", pk_c, msg.hex(), flipped.hex())

        # The all-zero G1 encoding with the compression bit set is the
        # identity: a well-formed encoding of a useless key, and accepting it
        # would let a signer claim membership with no key at all.
        inf_pk = "c0" + "00" * 47
        expect_err("identity public key", "verify", inf_pk, msg.hex(), sig_c)

        expect_err("truncated key", "verify", pk_c[:-2], msg.hex(), sig_c)
        expect_err("truncated signature", "verify", pk_c, msg.hex(), sig_c[:-2])

        # --- aggregation ---------------------------------------------------
        n = 3
        sks = [("%064x" % BLS.KeyGen(os.urandom(32))) for _ in range(n)]
        pks = [BLS.SkToPk(int(s, 16)).hex() for s in sks]
        msgs = [os.urandom(8 + i) for i in range(n)]
        sigs = [expect_ok("sign_i", "sign", sks[i], msgs[i].hex()) for i in range(n)]

        agg_c = expect_ok("aggregate", "aggregate", *sigs)
        agg_p = BLS.Aggregate([bytes.fromhex(s) for s in sigs]).hex()
        check(agg_c == agg_p, "Aggregate differs")

        av = ["aggregate_verify", agg_c]
        for i in range(n):
            av += [pks[i], msgs[i].hex()]
        expect_ok("aggregate_verify", *av)
        check(BLS.AggregateVerify([bytes.fromhex(p) for p in pks], msgs,
                                      bytes.fromhex(agg_c)),
              "py_ecc rejects our aggregate")

        # A DELIBERATE DIVERGENCE, recorded rather than hidden.
        #
        # Under the proof-of-possession scheme, py_ecc's AggregateVerify
        # accepts repeated messages, because verifying the proofs is supposed
        # to have happened already and is what makes it safe. Our
        # aggregate_verify does not take the proofs, so it cannot know that,
        # and a repeated message without them is forgeable. It refuses.
        #
        # The safe shared-message path is fast_aggregate_verify, which demands
        # the proofs in its signature and checks them.
        dup = ["aggregate_verify", agg_c]
        for i in range(n):
            dup += [pks[i], msgs[0].hex()]
        expect_err("repeated message without proofs", *dup)

        # --- proof of possession and the shared-message path ---------------
        pops = [expect_ok("pop_prove", "pop_prove", sks[i]) for i in range(n)]
        for i in range(n):
            check(pops[i] == BLS.PopProve(int(sks[i], 16)).hex(),
                  "PopProve differs")
            expect_ok("pop_verify", "pop_verify", pks[i], pops[i])
        # A proof only proves possession of ITS key.
        expect_err("pop under the wrong key", "pop_verify", pks[1], pops[0])

        shared = os.urandom(16)
        ssigs = [expect_ok("sign shared", "sign", sks[i], shared.hex()) for i in range(n)]
        sagg = expect_ok("aggregate shared", "aggregate", *ssigs)
        fav = ["fast_aggregate_verify", sagg, shared.hex()]
        for i in range(n):
            fav += [pks[i], pops[i]]
        expect_ok("fast_aggregate_verify", *fav)
        check(BLS.FastAggregateVerify(
                  [bytes.fromhex(p) for p in pks], shared, bytes.fromhex(sagg)),
              "py_ecc rejects our shared-message aggregate")

        # A bad proof must stop it, which is the whole point of demanding them.
        bad = bytearray(bytes.fromhex(pops[0]))
        bad[-1] ^= 1
        fav_bad = ["fast_aggregate_verify", sagg, shared.hex(), pks[0], bad.hex()]
        for i in range(1, n):
            fav_bad += [pks[i], pops[i]]
        expect_err("fast_aggregate_verify with a bad proof", *fav_bad)

    print(f"  {checks} checks, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
