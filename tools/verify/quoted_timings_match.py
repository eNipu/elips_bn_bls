#!/usr/bin/env python3
"""
Everything that quotes a timing must quote the one that was measured.

WHY THIS EXISTS. The same two numbers live in three places: bench/baseline.json,
which records them; README.md, which quotes them; and
bindings/js/demo/index.html, which quotes them at a reader who is watching the
library run. Nothing connected the three, so #37 and #38 updated two of them and
left the demo claiming 1.77 ms and 7.41 ms -- figures from before either change,
one of them flattering, on a public page.

A hand-maintained copy of a measured fact goes stale. This makes it fail the
build instead.

It checks the README too, and the machine name as well as the numbers, because
both have been wrong here. The README claimed "Intel Xeon at 2.10 GHz" while
baseline.json -- written programmatically from /proc/cpuinfo in the same
sitting -- said 2.80 GHz. The container had been moved to different hardware
between one turn and the next, and the machine was hand-typed from the older
reading. Numbers taken on one machine and attributed to another are worse than
no attribution, because they look reproducible.

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
README = ROOT / "README.md"
BASE = ROOT / "bench" / "baseline.json"
TOL_MS = 0.01


def main() -> int:
    base = json.loads(BASE.read_text())
    results = base.get("results")
    if not isinstance(results, dict):
        print(f"error: {BASE.name} has no results object", file=sys.stderr)
        return 2

    bad = []
    readme = README.read_text()

    # The machine, not only the numbers. Hand-typed once from a stale reading
    # after the container moved hosts mid-session.
    machine = base.get("machine", "")
    if not machine:
        print(f"error: {BASE.name} records no machine", file=sys.stderr)
        return 2
    if machine not in readme:
        bad.append(f"README does not name the machine {BASE.name} recorded "
                   f"({machine})")

    readme_rows = re.findall(r"^\| (sign|verify|pairing) \| ([\d.]+) ms", readme, re.M)
    if not readme_rows:
        print(f"error: no timing table found in README.md. If it moved, update "
              f"this script; do not delete it.", file=sys.stderr)
        return 2
    for op, ms in readme_rows:
        key = "pairing" if op == "pairing" else f"bls_{op}"
        if key not in results:
            bad.append(f"README {op}: {BASE.name} has no {key}")
            continue
        want = results[key]["median"] / 1000.0
        if abs(float(ms) - want) > TOL_MS:
            bad.append(f"README {op}: says {float(ms):.2f} ms, "
                       f"{BASE.name} measured {want:.2f} ms")

    html = DEMO.read_text()
    demo_rows = re.findall(
        r"<tr><td>(\w+)</td><td>([\d.]+)\s*ms</td><td>([\d.]+)\s*ms</td></tr>", html)
    if not demo_rows:
        print(f"error: no timing rows found in {DEMO.name}. If the table moved, "
              f"update this script; do not delete it.", file=sys.stderr)
        return 2
    for op, native_ms, _wasm_ms in demo_rows:
        key = f"bls_{op}"
        if key not in results:
            bad.append(f"demo {op}: {BASE.name} has no {key} to check against")
            continue
        want = results[key]["median"] / 1000.0
        if abs(float(native_ms) - want) > TOL_MS:
            bad.append(f"demo {op} native: page says {float(native_ms):.2f} ms, "
                       f"{BASE.name} measured {want:.2f} ms")

    print(f"checked the README ({len(readme_rows)} rows + machine) and "
          f"{len(demo_rows)} demo rows against {BASE.name}: {len(bad)} disagree")
    for b in bad:
        print(f"  {b}")
    if bad:
        print("\nRe-measure and update both, or something here is telling "
              "people what is not true.", file=sys.stderr)
        return 1

    # The wasm column has no recorded counterpart -- bench/baseline.json is the
    # native harness. It is checked by eye against bindings/js/bench.mjs when
    # the table is re-measured; stating that here so its absence is a known gap
    # rather than an oversight.
    return 0


if __name__ == "__main__":
    sys.exit(main())
