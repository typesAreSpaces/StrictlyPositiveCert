"""Compare the unchanged active benchmarking.mpl suite; only toggle factoring."""
import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ["-DENABLE_BINARY_SEARCH_AVKL", "-DENABLE_DEBUGGING",
         "-DEPS_FACTOR=17/20", "-DN_GUESS_AVKL=70", "-DN_GUESS_LS=600"]
DS = ["10", "5", "1", "1/2", "1/3"]
REPEATS = 7


def parse(output, enabled, repeat):
    if "FACTOR_SOS_BENCHMARKING_COMPLETE" not in output:
        raise RuntimeError("Benchmark did not finish")
    blocks = output.split("\n>> Test\n")[1:]
    if len(blocks) != 20:
        raise RuntimeError(f"Expected 20 method calls, got {len(blocks)}")
    rows = []
    for index, block in enumerate(blocks):
        case = index // 2
        elapsed = re.search(r"^>> Time taken\n([0-9.]+)", block, re.M)
        degree = re.search(r"^>> Degree size\n([0-9]+)", block, re.M)
        status = "success"
        if ">> Timeout" in block:
            status = "caught_failure_or_timeout"
        elif ">> RealCertify couldnt find certificates" in block:
            status = "no_certificate"
        if status == "success" and (not elapsed or not degree):
            raise RuntimeError(f"Missing result for call {index}")
        rows.append(dict(repeat=repeat, factor_sos_pre=enabled,
                         method="RealCertify" if index % 2 == 0 else "StrictlyPositiveCert",
                         d=DS[case // 2], k=1 if case % 2 == 0 else 3,
                         status=status,
                         cpu_seconds=float(elapsed[1]) if elapsed else None,
                         degree=int(degree[1]) if degree else None))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-first-pair", action="store_true",
                        help="reuse the initial off/on logs from this configuration")
    args = parser.parse_args()
    rows = []
    for repeat in range(1, REPEATS + 1):
        for enabled in ([False, True] if repeat % 2 else [True, False]):
            name = "on" if enabled else "off"
            log = Path(f"/tmp/factor_sos_benchmarking_{name}"
                       + ("" if repeat == 1 else f"_{repeat}") + ".log")
            cmd = ["maple", "-q", *FLAGS]
            if enabled:
                cmd.append("-DFACTOR_SOS_PRE")
            cmd.append("tests/run_factor_sos_benchmarking.mpl")
            started = time.monotonic()
            if repeat != 1 or not args.reuse_first_pair:
                with log.open("w") as stream:
                    subprocess.run(cmd, cwd=ROOT, stdout=stream,
                                   stderr=subprocess.STDOUT, check=True)
            rows.extend(parse(log.read_text(), enabled, repeat))
            print(f"Run {repeat}/{REPEATS}, flag {name}: complete "
                  f"({time.monotonic()-started:.1f}s)", flush=True)
    dest = ROOT / "tests/results"
    dest.mkdir(exist_ok=True)
    (dest / "benchmarking_factor_sos_pre_runs.json").write_text(
        json.dumps(dict(flags=FLAGS, repetitions=REPEATS, runs=rows), indent=2) + "\n")
    report = """# FACTOR_SOS_PRE comparison on tests/benchmarking.mpl

The benchmark file and its helpers were run unchanged. The only active suite is
`strictly_positive_polynomials_stengle_benchmark`: f=1-x^2+d, basis=[(1-x^2)^k],
d in [10, 5, 1, 1/2, 1/3], k in [1, 3]. Its existing 300-second limits remain.
Both configurations use the current Makefile settings: Averkov binary search and
debugging enabled, EPS_FACTOR=17/20, N_GUESS_AVKL=70, N_GUESS_LS=600.
WEIFENG_OPTIMIZATION and other optional flags remain disabled. The helper calls
the package method “Weifeng” regardless of this optimization flag.

Maple 2022 loaded package source through `tests/run_factor_sos_benchmarking.mpl`
before reading the actual benchmark file from the tests directory. No installed
.mla, benchmark case, or helper was changed. Debug output was redirected in both
configurations. Seven fresh-process runs per configuration were executed
sequentially, alternating order after the initial off/on pair. Tables use median
CPU seconds reported by the original helpers (not wall time); loading and output
after each timed construction are excluded. This differs from the previous
small-test harness, which measured warmed wall times. Degrees are the maximum
coefficient degrees reported by the existing `quickdegree` helper, not degrees
of individual square summands. Neither original helper independently validates
the certificate identity or coefficient nonnegativity during this benchmark.

## StrictlyPositiveCert results

Negative changes mean faster with factoring enabled.

| d | k | Off (s) | On (s) | Difference (s) | Change | Degree off → on | CPU range off / on (s) |
|---|---:|---:|---:|---:|---:|---:|---|
"""
    totals = [0.0, 0.0]
    for d in DS:
        for k in (1, 3):
            groups = [[r for r in rows if r["method"] == "StrictlyPositiveCert"
                       and r["d"] == d and r["k"] == k
                       and r["factor_sos_pre"] == flag] for flag in (False, True)]
            if any(r["status"] != "success" for group in groups for r in group):
                raise RuntimeError(f"Package failed for d={d}, k={k}")
            values = [[r["cpu_seconds"] for r in group] for group in groups]
            med = list(map(statistics.median, values))
            degrees = [sorted(set(r["degree"] for r in group)) for group in groups]
            totals = [a+b for a, b in zip(totals, med)]
            change = f"{100*(med[1]/med[0]-1):+.1f}%" if med[0] else "timer resolution"
            degree_text = " → ".join("/".join(map(str, ds)) for ds in degrees)
            ranges = " / ".join(f"{min(v):.3f}–{max(v):.3f}" for v in values)
            report += (f"| {d} | {k} | {med[0]:.3f} | {med[1]:.3f} | "
                       f"{med[1]-med[0]:+.3f} | {change} | {degree_text} | {ranges} |\n")
    report += (f"\nSum of per-case medians: **{totals[0]:.3f} → {totals[1]:.3f} s**, "
               f"**{100*(totals[1]/totals[0]-1):+.1f}%**. This is not total process wall time.\n")
    report += """
## RealCertify results

RealCertify is a separate package and does not use FACTOR_SOS_PRE. Its timings
provide context, not evidence of an effect of this flag. The benchmark prints
“Succeeds RealCertify” even when the helper returns no certificate; the status
below is taken from the helper's actual output.

| d | k | Off median (s) | On median (s) | Degree off → on | Status off / on |
|---|---:|---:|---:|---:|---|
"""
    for d in DS:
        for k in (1, 3):
            cells, degrees, statuses = [], [], []
            for flag in (False, True):
                group = [r for r in rows if r["method"] == "RealCertify"
                         and r["d"] == d and r["k"] == k and r["factor_sos_pre"] == flag]
                values = [r["cpu_seconds"] for r in group if r["cpu_seconds"] is not None]
                cells.append(f"{statistics.median(values):.3f}" if values else "—")
                degrees.append("/".join(map(str, sorted(set(
                    r["degree"] for r in group if r["degree"] is not None)))) or "—")
                statuses.append("/".join(sorted(set(r["status"] for r in group))))
            report += f"| {d} | {k} | {' | '.join(cells)} | {' → '.join(degrees)} | {' / '.join(statuses)} |\n"
    report += """
Small timing differences are noisy; min/max ranges are included above. Raw
per-run results are in `benchmarking_factor_sos_pre_runs.json`; logs remain in
`/tmp/factor_sos_benchmarking_{off,on}*.log`.

## Reproduce

Run from the repository root, adding only `-DFACTOR_SOS_PRE` for the enabled run:

```sh
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 tests/run_factor_sos_benchmarking.mpl
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 -DFACTOR_SOS_PRE tests/run_factor_sos_benchmarking.mpl
```
"""
    (dest / "benchmarking_factor_sos_pre.md").write_text(report)
    print(report, flush=True)


if __name__ == "__main__":
    main()
