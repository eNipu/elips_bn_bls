#!/usr/bin/env node
/**
 * The same hex-in, hex-out interface as tools/verify/bls_cli.c, so that
 * tools/verify/crosscheck_bls_pyecc.py runs against this binding UNCHANGED.
 *
 * That is the point: the JavaScript and the C are held to one bar by one
 * script, rather than by two that could drift apart. Every rejection case the
 * C has to refuse, this has to refuse, in the same words.
 */
import { load } from "./elips.mjs";

const hex = (b) => [...b].map((x) => x.toString(16).padStart(2, "0")).join("");
const unhex = (s) => {
  if (s.length % 2) { process.exit(2); }
  const out = new Uint8Array(s.length / 2);
  for (let i = 0; i < out.length; i++) out[i] = parseInt(s.substr(2 * i, 2), 16);
  return out;
};

const e = await load();
const [op, ...args] = process.argv.slice(2);

function guard(fn) {
  try { return fn(); }
  catch (err) {
    if (err instanceof RangeError || err instanceof TypeError) {
      console.log("ERR -1");
    } else {
      console.error(err.message);
      console.log("ERR -4");
    }
    process.exit(1);
  }
}

switch (op) {
  case "keygen":
    console.log(hex(guard(() => e.keygen(unhex(args[0])))));
    break;
  case "sk_to_pk":
    console.log(hex(guard(() => e.skToPk(unhex(args[0])))));
    break;
  case "sign":
    console.log(hex(guard(() => e.sign(unhex(args[0]), unhex(args[1])))));
    break;
  case "verify":
    guard(() => e.verify(unhex(args[0]), unhex(args[1]), unhex(args[2])));
    console.log("OK");
    break;
  case "aggregate":
    console.log(hex(guard(() => e.aggregate(args.map(unhex)))));
    break;
  case "aggregate_verify": {
    const sig = unhex(args[0]);
    const pks = [], msgs = [];
    for (let i = 1; i < args.length; i += 2) { pks.push(unhex(args[i])); msgs.push(unhex(args[i + 1])); }
    guard(() => e.aggregateVerify(pks, msgs, sig));
    console.log("OK");
    break;
  }
  case "fast_aggregate_verify": {
    const sig = unhex(args[0]), msg = unhex(args[1]);
    const pks = [], pops = [];
    for (let i = 2; i < args.length; i += 2) { pks.push(unhex(args[i])); pops.push(unhex(args[i + 1])); }
    guard(() => e.fastAggregateVerify(pks, pops, msg, sig));
    console.log("OK");
    break;
  }
  case "pop_prove":
    console.log(hex(guard(() => e.popProve(unhex(args[0])))));
    break;
  case "pop_verify":
    guard(() => e.popVerify(unhex(args[0]), unhex(args[1])));
    console.log("OK");
    break;
  default:
    console.error(`unknown op: ${op}`);
    process.exit(2);
}
