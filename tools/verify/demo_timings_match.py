#!/usr/bin/env python3
"""
The demo page quotes timings. They must be the ones that were measured.

WHY THIS EXISTS. The same two numbers live in three places: bench/baseline.json,
which records them; README.md, which quotes them; and
bindings/js/demo/index.html, which quotes them at a reader who is watching the
library run. Nothing connected the three, so #37 and #38 updated two of them and
left the demo claiming 1.77 ms and 7.41 ms -- figures from before either change,
one of them flattering, on a public page.

A hand-maintained copy of a measured fact goes stale. This makes it fail the
build instead.

Compared against baseline.json rather than against README.md, because
baseline.json is the recorded measurement and the README is itself a copy.
Tolerance is 0.01 ms: both are quoted to two decimals, so anything larger is a
real divergence rather than rounding.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEMO = ROOT / "bindings" / "js" / "demo" / "index.html"
BASE = ROOT / "bench" / "baseline.json"
TOL_MS = 0.01


def main() -> int:
    base = json.loads(BASE.read_text())
    results = base.get("results")
    if not isinstance(results, dict):
        print(f"error: {BASE.name} has no results object", file=sys.stderr)
        return 2

    html = DEMO.read_text()
    rows = re.findall(
        r"<tr><td>(\w+)</td><td>([\d.]+)\s*ms</td><td>([\d.]+)\s*ms</td></tr>", html)
    if not rows:
        # The table was restructured, or removed. Either way this check is no
        # longer looking at anything, which is the failure mode it exists to
        # prevent in the first place.
        print(f"error: no timing rows found in {DEMO.name}. If the table moved, "
              f"update this script; do not delete it.", file=sys.stderr)
        return 2

    bad = []
    for op, native_ms, _wasm_ms in rows:
        key = f"bls_{op}"
        if key not in results:
            bad.append(f"{op}: {BASE.name} has no {key} to check against")
            continue
        want = results[key]["median"] / 1000.0
        got = float(native_ms)
        if abs(got - want) > TOL_MS:
            bad.append(f"{op} native: page says {got:.2f} ms, "
                       f"{BASE.name} measured {want:.2f} ms")

    print(f"{len(rows)} timing rows in the demo, {len(bad)} disagree with "
          f"{BASE.name}")
    for b in bad:
        print(f"  {b}")
    if bad:
        print("\nRe-measure and update both, or the page is telling people "
              "something that is not true.", file=sys.stderr)
        return 1

    # The wasm column has no recorded counterpart -- bench/baseline.json is the
    # native harness. It is checked by eye against bindings/js/bench.mjs when
    # the table is re-measured; stating that here so its absence is a known
    # gap rather than an oversight.
    return 0


if __name__ == "__main__":
    sys.exit(main())
