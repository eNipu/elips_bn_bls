/*
 * The WebAssembly half of the performance table in README.md.
 *
 * It measures exactly what bench/bench.c measures for the native column --
 * bls_sign and bls_verify, bytes in and bytes out -- with the same interleaved
 * structure and the same reported statistics, so the two columns can honestly
 * be divided by one another. Measuring them with different harnesses is how a
 * ratio ends up describing the harnesses rather than the code.
 *
 *   node bench.mjs [--reps N] [--json]
 *
 * Run it on the SAME MACHINE as the native benchmark. The whole point of the
 * table is a like-for-like comparison, and a number carried over from another
 * machine silently destroys that.
 */
import { load, SK_BYTES } from "./elips.mjs";

const args = process.argv.slice(2);
const asJson = args.includes("--json");
const repsArg = args.indexOf("--reps");
const REPS = repsArg >= 0 ? Math.max(3, parseInt(args[repsArg + 1], 10)) : 9;
const INNER = 8;

const elips = await load();

const sk = elips.keygen();
const pk = elips.skToPk(sk);
const msg = new Uint8Array(32).fill(0x5a);
const sig = elips.sign(sk, msg);
elips.verify(pk, msg, sig);            // fail loudly here, not inside a timer

const cases = [
  ["bls_sign",   () => elips.sign(sk, msg)],
  ["bls_verify", () => elips.verify(pk, msg, sig)],
];

/* Interleaved, for the reason bench/bench.c gives at length: reps taken back
 * to back share a frequency and a cache state, so timing all of A then all of
 * B reports a spread that does not bound the run-to-run difference. One rep
 * times every case in turn. */
const samples = new Map(cases.map(([n]) => [n, []]));
for (let r = 0; r < REPS; r++) {
  for (const [name, fn] of cases) {
    const t0 = process.hrtime.bigint();
    for (let i = 0; i < INNER; i++) fn();
    const t1 = process.hrtime.bigint();
    samples.get(name).push(Number(t1 - t0) / 1000 / INNER);   // us per call
  }
}

const stats = (v) => {
  const s = [...v].sort((a, b) => a - b);
  const median = s[(s.length - 1) >> 1];
  return { median, min: s[0], spread: ((s[s.length - 1] - s[0]) / median) * 100 };
};

if (asJson) {
  const out = {};
  for (const [n, v] of samples) out[n] = stats(v);
  console.log(JSON.stringify({ runtime: `node ${process.version}`, reps: REPS, us: out }, null, 2));
} else {
  console.log(`WebAssembly  ${REPS} reps of ${INNER} calls  (node ${process.version})\n`);
  console.log("operation          median        min   spread");
  for (const [n, v] of samples) {
    const { median, min, spread } = stats(v);
    console.log(`${n.padEnd(16)} ${median.toFixed(2).padStart(9)} ${min.toFixed(2).padStart(10)} ${spread.toFixed(1).padStart(7)}%`);
  }
}
