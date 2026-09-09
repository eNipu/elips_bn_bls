"""Summarize preserved C timing samples; these are local measurements, not CIs."""
from collections import defaultdict
from pathlib import Path
from statistics import median
import argparse
import math


def parse(text):
    groups = defaultdict(lambda: defaultdict(dict))
    for line in text.splitlines():
        if not line.startswith("SAMPLE,"):
            continue
        _, kind, n, rep, method, time = line.split(",")
        bucket = groups[(kind, int(n))][method]
        value = float(time)
        if int(rep) in bucket or not math.isfinite(value) or value <= 0:
            raise ValueError("duplicate or invalid timing sample")
        bucket[int(rep)] = value
    if not groups:
        raise ValueError("no timing samples")
    return groups


def summarize(path):
    groups = parse(path.read_text())
    print(path.name)
    print("workload,n,control,control_us,fft4_us,median_time_reduction_pct,"
          "median_paired_reduction_pct")
    for (kind, n), methods in groups.items():
        candidate = methods["fft4-scaled4"]
        for control in ("repository-square", "cubic-scaled2"):
            baseline = methods[control]
            if set(baseline) != set(candidate):
                raise ValueError("unpaired timing samples")
            old, new = median(baseline.values()), median(candidate.values())
            paired = median(100*(1-candidate[i]/baseline[i]) for i in baseline)
            print(f"{kind},{n},{control},{old:.3f},{new:.3f},"
                  f"{100*(1-new/old):.3f},{paired:.3f}")


def selftest():
    text = "SAMPLE,square,0,0,repository-square,2\nSAMPLE,square,0,0,fft4-scaled4,1.8"
    groups = parse(text)
    if groups[("square", 0)]["fft4-scaled4"][0] != 1.8:
        raise AssertionError("timing parser")
    for invalid in ("", text+"\n"+text, "SAMPLE,square,0,0,fft4-scaled4,nan"):
        try:
            parse(invalid)
        except ValueError:
            continue
        raise AssertionError("invalid timing accepted")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", type=Path, nargs="*")
    args = parser.parse_args()
    selftest()
    for path in args.logs or sorted(Path(__file__).parent.glob("benchmark-*.log")):
        summarize(path)
