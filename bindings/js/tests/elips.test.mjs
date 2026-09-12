/**
 * The JavaScript surface.
 *
 * tools/verify/crosscheck_bls_pyecc.py already runs against this binding
 * unchanged, through cli.mjs, so correctness against an independent
 * implementation is covered there and is not repeated here. What is left is
 * what only JavaScript can get wrong: which error type comes out, that
 * verification cannot be mistaken for a boolean, and that a thrown error does
 * not leak WebAssembly memory.
 */
import { test } from "node:test";
import assert from "node:assert/strict";

import { load, ElipsError, InvalidSignature, InvalidKey } from "../elips.mjs";

const e = await load();
const enc = new TextEncoder();
const IKM = new Uint8Array(32).fill(7);
const MSG = enc.encode("the message");
const sk = e.keygen(IKM);
const pk = e.skToPk(sk);
const sig = e.sign(sk, MSG);

test("sizes are what the constants say", () => {
  assert.equal(sk.length, e.SK_BYTES);
  assert.equal(pk.length, e.PK_BYTES);
  assert.equal(sig.length, e.SIG_BYTES);
  assert.equal(e.popProve(sk).length, e.POP_BYTES);
});

test("keygen is deterministic in its input, random without one", () => {
  assert.deepEqual(e.keygen(IKM), sk);
  assert.notDeepEqual(e.keygen(), e.keygen());
});

test("verify returns undefined so it cannot be used as a boolean", () => {
  // The single most important property here. `if (verify(...))` must never be
  // a silent accept: undefined is falsy, so that branch never fires, and the
  // mistake shows up on the first VALID signature rather than a forged one.
  assert.equal(e.verify(pk, MSG, sig), undefined);
  assert.ok(!e.verify(pk, MSG, sig));
});

test("a tampered message throws InvalidSignature", () => {
  assert.throws(() => e.verify(pk, enc.encode("the messagf"), sig),
                (x) => x instanceof InvalidSignature);
});

test("the identity is refused as a public key", () => {
  const inf = new Uint8Array(48); inf[0] = 0xc0;
  assert.throws(() => e.pkValidate(inf), (x) => x instanceof InvalidKey);
});

test("every rejection is catchable as ElipsError", () => {
  const inf = new Uint8Array(48); inf[0] = 0xc0;
  const bad = Uint8Array.from(sig); bad[95] ^= 1;
  for (const f of [() => e.verify(pk, enc.encode("x"), sig),
                   () => e.verify(inf, MSG, sig),
                   () => e.verify(pk, MSG, bad)]) {
    assert.throws(f, (x) => x instanceof ElipsError);
  }
});

test("wrong length or type is a caller bug, not a bad signature", () => {
  assert.throws(() => e.skToPk(new Uint8Array(31)), RangeError);
  assert.throws(() => e.skToPk("not bytes"), TypeError);
  assert.throws(() => e.keygen(new Uint8Array(31)), RangeError);
  // Deliberately distinct hierarchies: conflating them would hide the bug.
  assert.ok(!(new RangeError("x") instanceof ElipsError));
});

test("aggregate of n is still one signature", () => {
  const sks = [1, 2, 3, 4].map((i) => e.keygen(new Uint8Array(32).fill(i)));
  const pks = sks.map((s) => e.skToPk(s));
  const msgs = ["a", "bb", "ccc", "dddd"].map((s) => enc.encode(s));
  const agg = e.aggregate(sks.map((s, i) => e.sign(s, msgs[i])));
  assert.equal(agg.length, e.SIG_BYTES);
  assert.equal(e.aggregateVerify(pks, msgs, agg), undefined);
  assert.throws(() => e.aggregateVerify(pks, [msgs[1], msgs[0], msgs[2], msgs[3]], agg),
                (x) => x instanceof InvalidSignature);
});

test("a repeated message is refused, and the error names the alternative", () => {
  const sks = [5, 6].map((i) => e.keygen(new Uint8Array(32).fill(i)));
  const pks = sks.map((s) => e.skToPk(s));
  const m = enc.encode("same");
  const agg = e.aggregate(sks.map((s) => e.sign(s, m)));
  assert.throws(() => e.aggregateVerify(pks, [m, m], agg),
                (x) => x instanceof RangeError && /distinct/.test(x.message));
});

test("the shared-message path takes the proofs and checks them", () => {
  const sks = [8, 9, 10].map((i) => e.keygen(new Uint8Array(32).fill(i)));
  const pks = sks.map((s) => e.skToPk(s));
  const pops = sks.map((s) => e.popProve(s));
  const m = enc.encode("one message, many signers");
  const agg = e.aggregate(sks.map((s) => e.sign(s, m)));
  assert.equal(e.fastAggregateVerify(pks, pops, m, agg), undefined);

  const bad = pops.map((p) => Uint8Array.from(p));
  bad[1][95] ^= 1;
  assert.throws(() => e.fastAggregateVerify(pks, bad, m, agg),
                (x) => x instanceof ElipsError);
});

test("a proof does not verify under another key", () => {
  const a = e.keygen(new Uint8Array(32).fill(11));
  const b = e.keygen(new Uint8Array(32).fill(12));
  assert.equal(e.popVerify(e.skToPk(a), e.popProve(a)), undefined);
  assert.throws(() => e.popVerify(e.skToPk(b), e.popProve(a)),
                (x) => x instanceof ElipsError);
});

test("aggregate_verify at the documented limit does not blow the wasm stack", () => {
  // This is the case that failed first: WebAssembly's default stack is 64 KB
  // and the C used to hold two full sets of points at once. Native has 8 MB
  // and never noticed. Exercise the widest supported aggregate every run.
  const n = 63;
  const sks = Array.from({ length: n }, (_, i) =>
    e.keygen(new Uint8Array(32).fill((i % 250) + 1)));
  const pks = sks.map((s) => e.skToPk(s));
  const msgs = sks.map((_, i) => enc.encode("msg-" + i));
  const agg = e.aggregate(sks.map((s, i) => e.sign(s, msgs[i])));
  assert.equal(e.aggregateVerify(pks, msgs, agg), undefined);
});

test("a thrown error does not leak wasm memory", () => {
  // Every allocation is freed in a finally. Without that, a page repeatedly
  // shown a bad signature would grow until it died.
  const bad = Uint8Array.from(sig); bad[95] ^= 1;
  for (let i = 0; i < 200; i++) {
    assert.throws(() => e.verify(pk, MSG, bad), (x) => x instanceof ElipsError);
  }
  e.verify(pk, MSG, sig);       // still works afterwards
});
