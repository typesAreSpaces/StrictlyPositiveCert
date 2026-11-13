"""Compare the four flag combinations on the active external Archimedean suite."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
from run_stage_division_benchmarking import ROOT, COMBINATIONS, TARGET_FLAGS, source_hashes

REPEATS = 7
METHODS = ["StrictlyPositiveCert", "RealCertify_without_arch", "RealCertify_with_arch"]


def hashes():
    result = source_hashes()
    for path in sorted((ROOT / "tests/benchmarks").glob("*.mpl")):
        result[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def parse(output, factor, division, repeat):
    if "FACTOR_SOS_BENCHMARKING_COMPLETE" not in output:
        raise RuntimeError("Benchmark did not finish")
    rows = []
    for block in output.split("\n>> Test\n")[1:]:
        match = re.match(r"Batch: (\d+) Test - (.+)\n", block)
        if not match:
            raise RuntimeError("Unexpected benchmark test label")
        method = {"Weifeng": "StrictlyPositiveCert",
                  "RealCertify without Archimedean polynomial": "RealCertify_without_arch",
                  "RealCertify with Archimedean polynomial": "RealCertify_with_arch"}[match[2]]
        elapsed = re.search(r"^>> Time taken\n([0-9.]+)", block, re.M)
        degree = re.search(r"^>> Degree size\n([0-9]+)", block, re.M)
        status = "success"
        if ">> Timeout" in block:
            status = "caught_failure_or_timeout"
        elif ">> RealCertify couldnt find certificates" in block:
            status = "no_certificate"
        if status == "success" and (not elapsed or not degree):
            raise RuntimeError(f"Missing result for {match[0]}")
        failed_time = re.search(r"^>> RealCertify couldnt find certificates\n([0-9.]+)", block, re.M)
        if elapsed is None:
            elapsed = failed_time
        poly = re.search(r"^>> input polynomial\n([^\n]+)", block, re.M)
        basis = re.search(r"^>> basis\n([^\n]+)", block, re.M)
        rows.append(dict(case=int(match[1]), method=method, repeat=repeat,
                         factor_sos_pre=factor, enable_stage_division=division,
                         status=status, cpu_seconds=float(elapsed[1]) if elapsed else None,
                         degree=int(degree[1]) if degree else None,
                         polynomial=poly[1] if poly else None,
                         basis=basis[1] if basis else None))
    counts = Counter((r["case"], r["method"]) for r in rows)
    cases = sorted({r["case"] for r in rows})
    if not cases or cases != list(range(1, max(cases)+1)):
        raise RuntimeError("Missing benchmark cases")
    if any(counts[(case, method)] != 1 for case in cases for method in METHODS):
        raise RuntimeError("Expected all three method calls for each case")
    return rows


def write_report(data, dest):
    rows = data["runs"]
    cases = sorted({r["case"] for r in rows})
    def group(case, combo, method="StrictlyPositiveCert"):
        return [r for r in rows if r["case"] == case and r["method"] == method
                and r["factor_sos_pre"] == combo[0] and r["enable_stage_division"] == combo[1]]
    def values(g):
        return [r["cpu_seconds"] for r in g if r["status"] == "success"]
    common = [case for case in cases if all(len(values(group(case, c))) == REPEATS
                                           for c in COMBINATIONS)]
    baseline_degrees = {case: {r["degree"] for r in group(case, COMBINATIONS[0])
                              if r["status"] == "success"} for case in common}
    report = f"""# Four-flag-combination study: updated Archimedean benchmark

Ran the current `tests/benchmarking.mpl` unchanged. Its active call is
`runExternalBenchmark("./benchmarks/arch.mpl")`, with **{len(cases)} cases**.
Each case invokes StrictlyPositiveCert and RealCertify with and without the
additional Archimedean generator. Existing 300-second call limits remain.
Only FACTOR_SOS_PRE (F) and ENABLE_STAGE_DIVISION (D) vary.

Fixed Makefile build command:

```text
{data['make_command']}
```

Seven fresh-process runs per configuration, sequential with rotating order.
Maple source is loaded directly using the existing wrapper before reading the
actual benchmark file. Neither the .mla nor benchmark cases/helpers were
changed. Source and data hashes are saved in the JSON and checked throughout.
Debug output stays enabled and is redirected for all runs. Timings are CPU
seconds reported by the original helpers, excluding loading and subsequent
certificate printing. Degrees are maximum SOS coefficient degrees reported by
the existing quickdegree helper, not degrees of individual square summands.
No extra verification flags were enabled; completion is not an independent
check of the certificate identity or coefficient nonnegativity.

## Aggregate package comparison

All aggregates below use the **{len(common)} cases** that succeeded in all seven
runs in every configuration. They are sums of per-case median CPU times, not
total process wall times. Small differences may reflect measurement noise.

| F | D | CPU time (s) | Change vs both off | Successful calls | Lower / equal / higher coefficient degrees vs both off |
|---|---|---:|---:|---:|---|
"""
    totals = [sum(statistics.median(values(group(case, combo))) for case in common)
              for combo in COMBINATIONS]
    summaries = []
    for combo, total in zip(COMBINATIONS, totals):
        successes = sum(r["status"] == "success" for r in rows
                        if r["method"] == "StrictlyPositiveCert"
                        and r["factor_sos_pre"] == combo[0]
                        and r["enable_stage_division"] == combo[1])
        changes = Counter()
        for case in common:
            degrees = {r["degree"] for r in group(case, combo)}
            if len(degrees) != 1 or len(baseline_degrees[case]) != 1:
                changes["variable"] += 1
            else:
                current, baseline = next(iter(degrees)), next(iter(baseline_degrees[case]))
                changes["lower" if current < baseline else "higher" if current > baseline else "equal"] += 1
        delta = 100*(total/totals[0]-1) if totals[0] else 0
        degrees_text = f"{changes['lower']} / {changes['equal']} / {changes['higher']}"
        if changes["variable"]:
            degrees_text += f" ({changes['variable']} variable)"
        report += (f"| {int(combo[0])} | {int(combo[1])} | {total:.3f} | {delta:+.1f}% | "
                   f"{successes}/{len(cases)*REPEATS} | {degrees_text} |\n")
        summaries.append(dict(factor=combo[0], division=combo[1], cpu_seconds=total,
                              change_percent=delta, degree_changes=dict(changes)))
    report += "\n## Per-case median CPU seconds and degrees\n\n"
    report += "Each cell is `seconds; coefficient degree`. Multiple degrees indicate variation between runs.\n\n"
    report += "| Case | F0 D0 | F1 D0 | F0 D1 | F1 D1 |\n|---:|---|---|---|---|\n"
    for case in cases:
        cells = []
        for combo in COMBINATIONS:
            g = group(case, combo)
            v = values(g)
            degrees = "/".join(map(str, sorted({r["degree"] for r in g if r["degree"] is not None})))
            cell = f"{statistics.median(v):.3f}; {degrees}" if v else "—"
            if len(v) != REPEATS:
                cell += f" ({len(v)}/{REPEATS} successful)"
            cells.append(cell)
        report += f"| {case} | " + " | ".join(cells) + " |\n"
    report += "\n## Package CPU ranges (minimum–maximum seconds)\n\n"
    report += "| Case | F0 D0 | F1 D0 | F0 D1 | F1 D1 |\n|---:|---|---|---|---|\n"
    for case in cases:
        cells = []
        for combo in COMBINATIONS:
            v = values(group(case, combo))
            cells.append(f"{min(v):.3f}–{max(v):.3f}" if v else "—")
        report += f"| {case} | " + " | ".join(cells) + " |\n"
    report += "\n## Method status counts\n\n"
    report += "RealCertify is independent of both flags. Its unconditional “Succeeds” print is ignored when the helper reports no certificate. Caught failures are labelled by the original benchmark as timeouts and may also represent errors.\n\n"
    report += "| F | D | Method | Status | Calls |\n|---|---|---|---|---:|\n"
    for key, count in sorted(Counter((r["factor_sos_pre"], r["enable_stage_division"],
                                      r["method"], r["status"]) for r in rows).items()):
        report += f"| {int(key[0])} | {int(key[1])} | {key[2]} | {key[3]} | {count} |\n"
    report += f"""
## Reproduce

```sh
python3 tests/run_arch_flag_benchmarking.py
```

The runner reads fixed flags from `make -Bn` and independently toggles only the
two requested flags. Raw per-run data and input polynomials are in
`benchmarking_arch_factor_stage_runs.json`. Raw logs: `{data['log_directory']}`.
"""
    (dest / "benchmarking_arch_factor_stage.md").write_text(report)
    print(json.dumps(dict(cases=len(cases), common_cases=len(common),
                          configurations=summaries), indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-off-log", type=Path)
    args = parser.parse_args()
    build = subprocess.check_output(["make", "-Bn"], cwd=ROOT, text=True).strip()
    command = shlex.split(build)
    if len(build.splitlines()) != 1 or command[0] != "maple":
        raise RuntimeError("Unexpected Makefile command")
    fixed = [flag for flag in command[1:-1] if flag not in TARGET_FLAGS]
    before = hashes()
    logdir = Path(tempfile.mkdtemp(prefix="arch_flag_benchmarking_"))
    dest = ROOT / "tests/results"
    dest.mkdir(exist_ok=True)
    data = dict(make_command=build, fixed_flags=fixed, repetitions=REPEATS,
                source_sha256=before, log_directory=str(logdir), runs=[])
    expected = None
    for repeat in range(1, REPEATS+1):
        shift = (repeat-1) % len(COMBINATIONS)
        for factor, division in COMBINATIONS[shift:] + COMBINATIONS[:shift]:
            label = f"F{int(factor)} D{int(division)} run {repeat}/{REPEATS}"
            log = logdir / f"f{int(factor)}_d{int(division)}_{repeat}.log"
            cmd = ["maple", "-q", *fixed]
            if factor: cmd.append("-DFACTOR_SOS_PRE")
            if division: cmd.append("-DENABLE_STAGE_DIVISION")
            cmd.append("tests/run_factor_sos_benchmarking.mpl")
            started = time.monotonic()
            if repeat == 1 and not factor and not division and args.first_off_log:
                log.write_text(args.first_off_log.read_text())
            else:
                with log.open("w") as stream:
                    process = subprocess.Popen(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
                    last_update = started
                    while process.poll() is None:
                        time.sleep(1)
                        if time.monotonic()-last_update >= 20:
                            labels = re.findall(r"Batch: (\d+) Test - ([^\n]+)", log.read_text())
                            print(f"{label}: latest {labels[-1] if labels else 'starting'}, "
                                  f"{time.monotonic()-started:.0f}s elapsed", flush=True)
                            last_update = time.monotonic()
                    if process.returncode:
                        raise RuntimeError(f"Maple exited with {process.returncode}; see {log}")
            parsed = parse(log.read_text(), factor, division, repeat)
            signature = [(r["case"], r["method"], r["polynomial"], r["basis"]) for r in parsed]
            if expected is None: expected = signature
            if signature != expected: raise RuntimeError("Benchmark inputs differ between runs")
            data["runs"].extend(parsed)
            (dest / "benchmarking_arch_factor_stage_runs.json").write_text(json.dumps(data, indent=2)+"\n")
            print(f"{label}: completed {len(parsed)//3} cases in "
                  f"{time.monotonic()-started:.1f}s", flush=True)
            if hashes() != before: raise RuntimeError("Sources/data changed during the benchmark")
    write_report(data, dest)


if __name__ == "__main__":
    main()
