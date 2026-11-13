"""Run benchmarking.mpl for all factoring/stage-division flag combinations."""
import hashlib
import json
from pathlib import Path
import shlex
import statistics
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
from run_factor_sos_benchmarking import ROOT, DS, REPEATS, parse

COMBINATIONS = [(False, False), (True, False), (False, True), (True, True)]
TARGET_FLAGS = {"-DFACTOR_SOS_PRE", "-DENABLE_STAGE_DIVISION"}


def source_hashes():
    paths = [ROOT / "Makefile", ROOT / "tests/benchmarking.mpl",
             ROOT / "tests/run_factor_sos_benchmarking.mpl"]
    paths += sorted((ROOT / "src").glob("*.mpl"))
    paths += sorted((ROOT / "tests/util").glob("*.mpl"))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def main():
    build = subprocess.check_output(["make", "-Bn"], cwd=ROOT, text=True).strip()
    command = shlex.split(build)
    if len(build.splitlines()) != 1 or command[0] != "maple":
        raise RuntimeError("Unexpected Makefile command; inspect before benchmarking")
    fixed = [flag for flag in command[1:-1] if flag not in TARGET_FLAGS]
    hashes = source_hashes()
    logdir = Path(tempfile.mkdtemp(prefix="stage_division_benchmarking_"))
    rows = []
    for repeat in range(1, REPEATS + 1):
        # Rotate sequential execution order to distribute startup/load effects.
        shift = (repeat - 1) % len(COMBINATIONS)
        order = COMBINATIONS[shift:] + COMBINATIONS[:shift]
        for factor, division in order:
            name = f"factor_{int(factor)}_division_{int(division)}"
            cmd = ["maple", "-q", *fixed]
            if factor:
                cmd.append("-DFACTOR_SOS_PRE")
            if division:
                cmd.append("-DENABLE_STAGE_DIVISION")
            cmd.append("tests/run_factor_sos_benchmarking.mpl")
            log = logdir / f"{name}_run_{repeat}.log"
            started = time.monotonic()
            with log.open("w") as stream:
                subprocess.run(cmd, cwd=ROOT, stdout=stream,
                               stderr=subprocess.STDOUT, check=True)
            parsed = parse(log.read_text(), factor, repeat)
            for row in parsed:
                row["enable_stage_division"] = division
            rows.extend(parsed)
            print(f"Run {repeat}/{REPEATS}: {name}, "
                  f"{time.monotonic()-started:.1f}s", flush=True)
            if source_hashes() != hashes:
                raise RuntimeError("Sources changed during measurement; do not mix results")
    dest = ROOT / "tests/results"
    dest.mkdir(exist_ok=True)
    data = dict(make_command=build, fixed_flags=fixed, repetitions=REPEATS,
                source_sha256=hashes, log_directory=str(logdir), runs=rows)
    (dest / "benchmarking_factor_stage_runs.json").write_text(
        json.dumps(data, indent=2) + "\n")

    def group(method, d, k, combo):
        return [r for r in rows if r["method"] == method and r["d"] == d
                and r["k"] == k and r["factor_sos_pre"] == combo[0]
                and r["enable_stage_division"] == combo[1]]

    def median_text(group):
        values = [r["cpu_seconds"] for r in group if r["status"] == "success"]
        return f"{statistics.median(values):.3f}" if values else "—"

    def degree_text(group):
        degrees = sorted({r["degree"] for r in group if r["degree"] is not None})
        return "/".join(map(str, degrees)) or "—"

    report = """# Factoring and stage division: four-configuration benchmark

Ran the unchanged active suite in `tests/benchmarking.mpl`: ten Stengle cases
f=1-x^2+d, basis=[(1-x^2)^k], d in [10, 5, 1, 1/2, 1/3], k in [1, 3].
The file's existing 300-second limits and RealCertify calls remain unchanged.
Only FACTOR_SOS_PRE and ENABLE_STAGE_DIVISION vary. All other compiler flags
and values come from the current Makefile's dry-run command:

```text
""" + build + "\n```\n\n"
    report += """Maple 2022; seven fresh processes per combination, run sequentially with
rotating configuration order. Source is loaded through the existing wrapper,
then the actual benchmark file is read from the tests directory. Installed
.mla files and benchmark source/helpers are untouched. Source hashes were
checked throughout; they and per-run results are saved in the JSON artifact.
Debug output remains enabled and is redirected to logs for every run.

Times are median CPU seconds reported by the existing checkWeifeng helper,
excluding source loading and subsequent certificate printing. Degrees are
the maximum certificate coefficient degrees reported by the existing
quickdegree helper, not degrees of individual square summands. No extra
certificate verification flags were enabled. The original benchmark does
not independently verify identities/nonnegativity. Small timing differences
are noisy; per-case ranges appear below.

Columns identify the two flags as F=FACTOR_SOS_PRE, D=ENABLE_STAGE_DIVISION.

## StrictlyPositiveCert median CPU times (seconds)

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---:|---:|---:|---:|
"""
    totals = [0.0] * 4
    failures = []
    for d in DS:
        for k in (1, 3):
            groups = [group("StrictlyPositiveCert", d, k, c) for c in COMBINATIONS]
            report += f"| {d} | {k} | " + " | ".join(map(median_text, groups)) + " |\n"
            for i, g in enumerate(groups):
                values = [r["cpu_seconds"] for r in g if r["status"] == "success"]
                if len(values) == REPEATS:
                    totals[i] += statistics.median(values)
                else:
                    failures.append((d, k, COMBINATIONS[i], len(values)))
    report += "\n## Certificate coefficient degrees\n\n"
    report += "| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |\n|---|---:|---:|---:|---:|---:|\n"
    for d in DS:
        for k in (1, 3):
            report += f"| {d} | {k} | " + " | ".join(
                degree_text(group("StrictlyPositiveCert", d, k, c))
                for c in COMBINATIONS) + " |\n"
    report += "\n## Aggregate\n\n"
    report += "| Configuration | Sum of per-case medians (s) | Change vs F0 D0 |\n|---|---:|---:|\n"
    for combo, total in zip(COMBINATIONS, totals):
        report += (f"| F{int(combo[0])} D{int(combo[1])} | {total:.3f} | "
                   f"{100*(total/totals[0]-1):+.1f}% |\n")
    report += "\nThese aggregates are sums of per-case medians, not total process wall times.\n"
    report += f"\nPackage failures: {failures if failures else 'none; all 280 package calls completed successfully'}.\n"
    report += "\n## Package CPU ranges (minimum–maximum seconds)\n\n"
    report += "| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |\n|---|---:|---:|---:|---:|---:|\n"
    for d in DS:
        for k in (1, 3):
            cells = []
            for combo in COMBINATIONS:
                values = [r["cpu_seconds"] for r in group("StrictlyPositiveCert", d, k, combo)
                          if r["status"] == "success"]
                cells.append(f"{min(values):.3f}–{max(values):.3f}" if values else "—")
            report += f"| {d} | {k} | " + " | ".join(cells) + " |\n"
    report += """
## RealCertify

This separate package does not use either flag. Each cell below gives median
CPU seconds and certificate degree for its successful calls, or failure status.
The script's unconditional “Succeeds RealCertify” messages are not treated as
proof of success; the parser checks for the helper's no-certificate output.

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---|---|---|---|
"""
    for d in DS:
        for k in (1, 3):
            cells = []
            for combo in COMBINATIONS:
                g = group("RealCertify", d, k, combo)
                statuses = "/".join(sorted({r["status"] for r in g}))
                cells.append(f"{median_text(g)} s, degree {degree_text(g)} ({statuses})")
            report += f"| {d} | {k} | " + " | ".join(cells) + " |\n"
    report += "\n## Reproduce\n\n```sh\npython3 tests/run_stage_division_benchmarking.py\n```\n"
    report += "\nThe runner reads the fixed flags from `make -Bn` and independently toggles only the two requested flags.\n"
    report += f"\nRaw logs: `{logdir}`. Per-run data: `benchmarking_factor_stage_runs.json`.\n"
    (dest / "benchmarking_factor_stage.md").write_text(report)
    print(report, flush=True)


if __name__ == "__main__":
    main()
