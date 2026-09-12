/* Split out of index.html so the page carries no inline <script>.
 * See the note at the top of style.css. */
import { load, ElipsError } from "../elips.mjs";

const $ = (id) => document.getElementById(id);
const enc = new TextEncoder();
const hex = (b) => [...b].map((x) => x.toString(16).padStart(2, "0")).join("");

let elips, signers = [], signedText = null, agg = null;

const shortKey = (pk) => hex(pk).slice(0, 8) + "…" + hex(pk).slice(-4);

function render() {
  const n = signers.length;
  $("n").textContent = n;
  $("sep").innerHTML = (n * 96).toLocaleString() + ' <span class="unit">bytes</span>';
  $("agg").innerHTML = n ? '96 <span class="unit">bytes</span>' : "&mdash;";
  $("saved").innerHTML = n > 1
    ? (100 - 100 / n).toFixed(n > 20 ? 1 : 0) + ' <span class="unit">%</span>'
    : "&mdash;";
  $("signers").innerHTML = signers
    .map((s, i) => `<span>${String(i + 1).padStart(2, "0")} ${shortKey(s.pk)}</span>`)
    .join("");
  $("sighex").textContent = agg ? hex(agg) : "—";

  // The input carries a mark once the text in it is the text that was signed,
  // so "I edited this and did not re-sign" is visible before you press verify.
  const stmt = $("stmt");
  if (signedText !== null && stmt.value === signedText) stmt.setAttribute("data-signed", "");
  else stmt.removeAttribute("data-signed");

  for (const id of ["verify", "tamperSig", "tamperMsg"]) $(id).disabled = !n;
}

function show(ok, title, detail) {
  $("result").innerHTML =
    `<span class="badge ${ok ? "ok" : "fail"}">${ok ? "✓" : "✗"} ${title}</span>` +
    `<p class="detail">${detail}</p>`;
}

function clearVerdict() { $("result").innerHTML = ""; }

function signAll(count) {
  const text = $("stmt").value;
  // Everyone must have signed the SAME bytes for the aggregate to mean
  // anything, so adding a signer starts over if the text changed.
  if (signedText !== null && signedText !== text) signers = [];
  signedText = text;
  const msg = enc.encode(text);
  for (let i = 0; i < count; i++) {
    const sk = elips.keygen();
    signers.push({ sk, pk: elips.skToPk(sk), pop: elips.popProve(sk),
                   sig: elips.sign(sk, msg) });
  }
  agg = elips.aggregate(signers.map((s) => s.sig));
  clearVerdict();
  render();
}

function verify() {
  const msg = enc.encode($("stmt").value);
  const t0 = performance.now();
  try {
    // fastAggregateVerify takes the proofs of possession and checks them. It
    // will not run without them: summing public keys over a shared message is
    // forgeable otherwise.
    elips.fastAggregateVerify(signers.map((s) => s.pk),
                              signers.map((s) => s.pop), msg, agg);
    const ms = (performance.now() - t0).toFixed(0);
    show(true, `valid — all ${signers.length} signed this exact text`,
         `Checked in ${ms} ms, from 96 bytes, against ${signers.length} public ` +
         `key${signers.length > 1 ? "s" : ""} at once.`);
  } catch (err) {
    const ms = (performance.now() - t0).toFixed(0);
    const why = err instanceof ElipsError
      ? "The aggregate does not match these keys and this text."
      : err.message;
    show(false, "rejected", `${why} (${ms} ms)`);
  }
  render();
}

/* Any uncaught throw in a handler leaves the page looking like the button did
 * nothing, which on this page is indistinguishable from a check that passed.
 * Every handler reports instead. */
function guarded(fn) {
  return () => {
    try { fn(); }
    catch (err) { show(false, "something went wrong", err.message); }
  };
}

$("add").onclick = guarded(() => signAll(1));
$("add10").onclick = guarded(() => signAll(10));
$("reset").onclick = guarded(() => {
  signers = []; agg = null; signedText = null;
  clearVerdict(); render();
});
$("verify").onclick = guarded(verify);
$("stmt").oninput = render;

$("tamperMsg").onclick = guarded(() => {
  // Change the text WITHOUT re-signing: the signatures now cover something
  // else. This is the case a demo where everything always works would hide.
  const stmt = $("stmt");
  stmt.value = stmt.value.replace(/Friday/, "Monday");
  if (stmt.value === signedText) stmt.value += " (edited)";
  verify();
});

$("tamperSig").onclick = guarded(() => {
  const i = Math.floor(Math.random() * signers.length);
  const bad = Uint8Array.from(signers[i].sig);
  bad[95] ^= 1;                                   // one bit
  try {
    agg = elips.aggregate(signers.map((s, j) => (j === i ? bad : s.sig)));
  } catch {
    // Most single-bit changes stop the 96 bytes being a point on the curve at
    // all, so aggregation refuses them before any pairing runs. That is a real
    // and useful outcome, not an error to swallow: showing nothing here would
    // be a demo whose failure case is invisible.
    render();
    show(false, "rejected before verifying",
         `One bit was flipped in signer ${String(i + 1).padStart(2, "0")}'s ` +
         `signature. Those 96 bytes are no longer a point on the curve, so the ` +
         `signatures cannot even be combined. Nothing reached the pairing check.`);
    return;
  }
  render();
  verify();
});

load().then((api) => {
  elips = api;
  const status = $("status");
  status.className = "status ready";
  status.textContent = "ready — running in this tab, no backend";
  $("app").hidden = false;
  signAll(3);
}).catch((err) => {
  $("status").textContent = "failed to load the WebAssembly module: " + err.message;
});
