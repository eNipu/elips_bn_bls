#!/usr/bin/env python3
"""Decide whether a change made the library slower, on a machine that is noisy.

Two modes, and the first is the one to use.

  A/B, alternated (recommended)

      python3 bench/compare.py --ab ./bench_old ./bench_new

    Runs the two binaries alternately, one round each, and compares them
    pairwise. This is the only method here that was shown to work: a planted
    12% slowdown was invisible to everything else tried and came out at +9.1%
    under alternation, with every individual pair agreeing on the direction.

  Two stored JSON files

      python3 bench/compare.py before.json after.json

    For baselines recorded at different times. Weaker, and it says so: it
    subtracts the median change across all operations as machine drift, which
    only works if most operations did not change.

WHY ALTERNATION

Absolute timings drift. On the machine this was written on, every one of twelve
operations moved 19-23% between three consecutive runs of the SAME binary.
Worse, the drift is not uniform across operations, so normalising against one
reference operation just imports that operation's noise into every number. Both
were tried and both failed to see a real 12% regression.

Alternation works because the two binaries meet the same machine within seconds
of each other, over and over. Drift moves both members of a pair together, and
the comparison is of pairs.

RUN THIS ON AN IDLE MACHINE. ABBA cancels drift that is smooth across a round;
it does not cancel a load that arrives partway through the sequence. Running a
compile alongside one of these produced +9 to +13% on all twelve operations at
100% agreement, including several the change could not possibly have touched.
100% agreement on everything at once is the signature of that mistake, not of a
real regression: a real one moves one or two operations, not the whole table.

Exit status is 1 if anything regressed, so this can gate a pull request.
"""
import argparse
import json
import statistics
import subprocess
import sys


def run_binary(path, reps, name=None):
    cmd = [path, "--json", "--reps", str(reps)]
    if name:
        cmd.append(name)
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return json.loads(out)["results"]


def ab(old, new, rounds, reps, floor, abs_floor):
    """Alternate the two binaries and compare pairwise."""
    # ABBA ordering. Running old-then-new every round biases the second one by
    # whatever the machine is doing over the round: tried that way, a planted
    # regression in ONE routine reported ten others as slower too, all at +3 to
    # +7% while the real one was +15.8%. Swapping the order on alternate rounds
    # cancels drift that is linear over a round, which is most of it.
    pairs = {}
    for i in range(rounds):
        if i % 2 == 0:
            a = run_binary(old, reps)
            b = run_binary(new, reps)
        else:
            b = run_binary(new, reps)
            a = run_binary(old, reps)
        for k in set(a) & set(b):
            pairs.setdefault(k, []).append((a[k]["median"], b[k]["median"]))
        print(f"  round {i + 1}/{rounds} done", file=sys.stderr)

    if not pairs:
        sys.exit("the two binaries share no operations")

    width = max(len(k) for k in pairs)
    print(f"\nA/B alternated, {rounds} rounds of {reps} reps\n")
    print(f"{'operation':<{width}} {'old':>10} {'new':>10} {'change':>9} {'agree':>7}   verdict")

    regressed = 0
    for k in sorted(pairs):
        ratios = [(b - a) / a * 100.0 for a, b in pairs[k]]
        change = statistics.median(ratios)
        # How many pairs agree with the median's sign. A real change moves
        # every pair the same way; noise splits them.
        same = sum(1 for r in ratios if (r > 0) == (change > 0))
        agree = same / len(ratios)

        oa = statistics.median([a for a, _ in pairs[k]])
        ob = statistics.median([b for _, b in pairs[k]])

        # An operation below the absolute floor cannot be measured here. On BN
        # ep_in_subgroup is a curve equation and runs in 0.8 us; a 0.1 us
        # difference read as "+5.4% SLOWER" is the timer, not the code.
        if min(oa, ob) < abs_floor:
            verdict = "too fast to measure"
        elif abs(change) <= floor or agree < 0.75:
            verdict = "unchanged"
        elif change > 0:
            verdict = "SLOWER"
            regressed += 1
        else:
            verdict = "faster"

        print(f"{k:<{width}} {oa:>10.1f} {ob:>10.1f} {change:>+8.1f}% {agree:>6.0%}   {verdict}")
    return regressed


def files(before_path, after_path, floor, abs_floor):
    """Compare two recorded runs, correcting for drift with the median change."""
    with open(before_path) as f:
        db = json.load(f)
    with open(after_path) as f:
        da = json.load(f)
    if db["curve"] != da["curve"]:
        sys.exit(f"different curves: {db['curve']} vs {da['curve']}")
    before, after = db["results"], da["results"]

    shared = sorted(set(before) & set(after))
    if not shared:
        sys.exit("the two files share no operations")

    change = {k: (after[k]["median"] - before[k]["median"]) / before[k]["median"] * 100.0
              for k in shared}
    drift = statistics.median(change.values())

    width = max(len(k) for k in shared)
    print(f"{db['curve']}   two recorded runs, machine drift {drift:+.1f}% removed")
    print("NOTE: weaker than --ab. Only trust a result larger than the drift.\n")
    print(f"{'operation':<{width}} {'before':>10} {'after':>10} {'change':>9}   verdict")

    regressed = 0
    for k in shared:
        c = change[k] - drift
        noise = max(before[k]["spread_pct"], after[k]["spread_pct"], floor)
        if min(before[k]["median"], after[k]["median"]) < abs_floor:
            verdict = "too fast to measure"
            print(f"{k:<{width}} {before[k]['median']:>10.1f} {after[k]['median']:>10.1f} "
                  f"{c:>+8.1f}%   {verdict}")
            continue
        if abs(c) <= noise:
            verdict = "unchanged"
        elif c > 0:
            verdict = "SLOWER"
            regressed += 1
        else:
            verdict = "faster"
        print(f"{k:<{width}} {before[k]['median']:>10.1f} {after[k]['median']:>10.1f} "
              f"{c:>+8.1f}%   {verdict}")
    return regressed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("first", help="old binary with --ab, else the before JSON")
    ap.add_argument("second", help="new binary with --ab, else the after JSON")
    ap.add_argument("--ab", action="store_true", help="alternate two binaries (recommended)")
    ap.add_argument("--rounds", type=int, default=8,
                    help="A/B rounds, even is better for the ABBA ordering (default 8)")
    ap.add_argument("--reps", type=int, default=9, help="reps inside each run (default 9)")
    ap.add_argument("--floor", type=float, default=3.0,
                    help="smallest change treated as real, percent (default 3)")
    ap.add_argument("--abs-floor", type=float, default=2.0,
                    help="ignore operations faster than this many microseconds "
                         "(default 2); below it the timer noise is the reading")
    args = ap.parse_args()

    regressed = (ab(args.first, args.second, args.rounds, args.reps,
                    args.floor, args.abs_floor)
                 if args.ab else files(args.first, args.second,
                                       args.floor, args.abs_floor))

    print()
    print(f"{regressed} operation(s) regressed" if regressed
          else "nothing regressed")
    return 1 if regressed else 0


if __name__ == "__main__":
    sys.exit(main())
