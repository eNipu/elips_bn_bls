/**
 * The demo page, driven in a real browser.
 *
 * Every claim this test makes is one the page makes to its reader, and a page
 * that quietly stops keeping them is worse than no page. In particular:
 *
 *   - "There are no network requests after the page loads. Open your
 *     browser's network panel and watch." So this counts them.
 *   - The tamper buttons must visibly fail. A demo where everything always
 *     succeeds teaches that the check does not matter. This was not
 *     hypothetical: "Corrupt one signature" threw out of its handler and the
 *     page silently did nothing, because flipping a bit usually stops the 96
 *     bytes being a point on the curve and aggregation refuses it before any
 *     pairing runs. Driving the page found it; reading it had not.
 *
 * Skipped when playwright is not installed, so `npm test` still works without
 * a browser download.
 */
import { test, before, after } from "node:test";
import assert from "node:assert/strict";
import { createServer } from "node:http";
import { readFileSync, existsSync } from "node:fs";
import { extname, join, normalize, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = join(HERE, "..");
const PORT = 8739;

let chromium;
try { ({ chromium } = await import("playwright")); }
catch { /* not installed */ }

const maybe = chromium ? test : test.skip;

// A wrong content type is not a 404: the browser fetches the file and then
// refuses to use it. Serving style.css as octet-stream left the page unstyled
// and the overflow assertion below is what noticed.
const TYPES = { ".html": "text/html", ".mjs": "text/javascript",
                ".js": "text/javascript", ".css": "text/css",
                ".wasm": "application/wasm" };

let server, browser, page, requestsAfterLoad, pageErrors;

before(async () => {
  if (!chromium) return;
  server = createServer((req, res) => {
    const rel = normalize(decodeURI(req.url.split("?")[0]));
    const p = join(ROOT, rel);
    if (rel.includes("..") || !existsSync(p)) { res.writeHead(404); return res.end(); }
    res.writeHead(200, { "content-type": TYPES[extname(p)] || "application/octet-stream" });
    res.end(readFileSync(p));
  });
  await new Promise((r) => server.listen(PORT, r));

  const launch = {};
  // The container ships browsers at a pinned path; use them rather than
  // downloading a second copy.
  for (const p of ["/opt/pw-browsers/chromium-1194/chrome-linux/chrome"]) {
    if (existsSync(p)) launch.executablePath = p;
  }
  browser = await chromium.launch(launch);
  page = await browser.newPage({ viewport: { width: 390, height: 844 } });

  pageErrors = [];
  page.on("pageerror", (e) => pageErrors.push(e.message));
  page.on("console", (m) => { if (m.type() === "error") pageErrors.push(m.text()); });

  await page.goto(`http://localhost:${PORT}/demo/index.html`);
  await page.waitForSelector("#app:not([hidden])", { timeout: 60000 });
  requestsAfterLoad = [];
  page.on("request", (r) => requestsAfterLoad.push(r.url()));
});

after(async () => {
  await browser?.close();
  server?.close();
});

const text = (s) => page.textContent(s).then((t) => t.trim());

maybe("it loads and signs without being asked twice", async () => {
  assert.equal(await text("#n"), "3");
  assert.match(await text("#status"), /ready/i);
});

maybe("the aggregate stays 96 bytes as signers are added", async () => {
  await page.click("#add10");
  await page.waitForFunction(() => document.getElementById("n").textContent === "13");
  assert.equal(await text("#n"), "13");
  assert.match(await text("#sep"), /1,248/);
  assert.match(await text("#agg"), /^96/);
});

maybe("a genuine aggregate verifies", async () => {
  await page.click("#verify");
  await page.waitForSelector("#result .badge");
  assert.ok(await page.locator("#result .badge.ok").count(),
            "the verdict should be marked ok");
  assert.match(await text("#result .badge"), /valid/i);
});

maybe("changing the statement after signing visibly fails", async () => {
  await page.click("#tamperMsg");
  await page.waitForSelector("#result .badge.fail");
  assert.match(await text("#result .badge"), /rejected/i);
});

maybe("corrupting one signature visibly fails", async () => {
  // The case that silently did nothing before it was driven.
  await page.click("#reset");
  await page.waitForFunction(() => document.getElementById("n").textContent === "0");
  await page.click("#add10");
  await page.waitForFunction(() => document.getElementById("n").textContent === "10",
                             null, { timeout: 60000 });
  await page.click("#tamperSig");
  await page.waitForSelector("#result .badge", { timeout: 60000 });
  assert.ok(await page.locator("#result .badge.fail").count(),
            "the verdict should be marked failed");
  assert.match(await text("#result .badge"), /rejected/i);
});

maybe("no network requests are made after load", async () => {
  assert.deepEqual(requestsAfterLoad, [],
    "the page tells the reader to check this in devtools");
});

maybe("the stylesheet actually applied", async () => {
  // Cheap, and it is the difference between "styled" and "the browser fetched
  // the CSS and declined to use it", which look identical to a 200 check.
  // A rule that only the stylesheet supplies. Unstyled, a <section> has no
  // border at all, so this separates "styled" from "the browser fetched the
  // CSS and declined to use it" -- which look identical to a 200 check.
  const border = await page.evaluate(() =>
    getComputedStyle(document.querySelector("section")).borderTopWidth);
  assert.equal(border, "1px", "the section rule did not take effect");
});

maybe("no horizontal overflow at phone width", async () => {
  const over = await page.evaluate(() =>
    document.documentElement.scrollWidth > window.innerWidth);
  assert.equal(over, false);
});

maybe("the page states what it is not", async () => {
  const body = await page.textContent("footer");
  assert.match(body, /Not audited/i);
  assert.match(body, /slower than the native library/i);
});

maybe("nothing threw", async () => {
  assert.deepEqual(pageErrors, []);
});
