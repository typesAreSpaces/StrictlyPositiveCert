"""Compare main, the working tree, and the working tree with FACTOR_SOS_PRE.

All variants run the current tests/benchmarking.mpl unchanged (main's src/ is
extracted into a scratch tree next to a copy of the working tree's tests/), with
the fixed flags from `make -Bn`. After timing, every certificate is recomputed
and checked exactly by tests/verify_external_certificates.mpl, which also gives
the exact expanded degrees (quickdegree does not see cancellation).
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
from run_arch_flag_benchmarking import hashes
from run_stage_division_benchmarking import ROOT

REPEATS = 7
METHODS = {"Weifeng": "StrictlyPositiveCert",
           "RealCertify without Archimedean polynomial": "RealCertify_without_arch",
           "RealCertify with Archimedean polynomial": "RealCertify_with_arch"}


def main_tree(ref):
    tree = Path(tempfile.mkdtemp(prefix="lbp_main_tree_"))
    archive = subprocess.run(["git", "archive", ref, "src"], cwd=ROOT,
                             check=True, capture_output=True).stdout
    subprocess.run(["tar", "-x", "-C", str(tree)], input=archive, check=True)
    shutil.copytree(ROOT / "tests", tree / "tests", ignore=shutil.ignore_patterns("results"))
    return tree


def active_suites():
    text = (ROOT / "tests/benchmarking.mpl").read_text()
    return re.findall(r'^runExternalBenchmark\("\./benchmarks/([\w]+)\.mpl"\)', text, re.M)


def parse(output, suites, variant, repeat):
    if "FACTOR_SOS_BENCHMARKING_COMPLETE" not in output:
        raise RuntimeError("Benchmark did not finish")
    segments = output.split("\n>> Start benchmark\n")[1:]
    if len(segments) != len(suites):
        raise RuntimeError(f"Expected {len(suites)} suites, found {len(segments)}")
    rows = []
    for suite, segment in zip(suites, segments):
        for block in segment.split("\n>> Test\n")[1:]:
            match = re.match(r"Batch: (\d+) Test - (.+)\n", block)
            if not match:
                raise RuntimeError("Unexpected benchmark test label")
            elapsed = re.search(r"^>> Time taken\n([0-9.]+)", block, re.M)
            if elapsed is None:
                elapsed = re.search(r"^>> RealCertify couldnt find certificates\n([0-9.]+)", block, re.M)
            degree = re.search(r"^>> Degree size\n([0-9]+)", block, re.M)
            status = "success"
            if ">> Timeout" in block:
                status = "caught_failure_or_timeout"
            elif ">> RealCertify couldnt find certificates" in block:
                status = "no_certificate"
            poly = re.search(r"^>> input polynomial\n([^\n]+)", block, re.M)
            rows.append(dict(suite=suite, case=int(match[1]), method=METHODS[match[2]],
                             variant=variant, repeat=repeat, status=status,
                             cpu_seconds=float(elapsed[1]) if elapsed else None,
                             quickdegree=int(degree[1]) if degree else None,
                             polynomial=poly[1] if poly else None))
    counts = Counter((r["suite"], r["case"], r["method"]) for r in rows)
    if any(n != 1 for n in counts.values()) or \
       len(counts) != len({(r["suite"], r["case"]) for r in rows})*len(METHODS):
        raise RuntimeError("Expected all three method calls for each case")
    return rows


def parse_degrees(text):
    degrees = [None if d == "-infinity" else int(d) for d in text.strip("[]").split(", ")]
    return degrees, max(d for d in degrees if d is not None)


def verify(root, flags, suites, logdir, variant, timeout=300):
    """valid is True/False, or None when the certificate was not computed or not checked."""
    log = logdir / f"{variant}_verify.log"
    with log.open("w") as stream:
        subprocess.run(["maple", "-q", *flags, "tests/verify_external_certificates.mpl"],
                       cwd=root, stdout=stream, stderr=subprocess.STDOUT, check=True,
                       env=dict(os.environ, BENCH_SUITES=",".join(suites),
                                VERIFY_TIMEOUT=str(timeout)))
    result = {}
    for line in log.read_text().splitlines():
        if not line.startswith("VERIFY "):
            continue
        parts = line.split(" ", 4)
        key = (parts[1], int(parts[2]))
        if parts[3] in ("error", "timeout"):
            result[key] = dict(valid=None if parts[3] == "timeout" else False,
                               status=parts[3], degrees=None)
        elif parts[3] == "check_timeout":
            degrees, top = parse_degrees(parts[4])
            result[key] = dict(valid=None, status="check_timeout", degrees=degrees, max_degree=top)
        else:
            flag, degree_text = parts[4].split(" ", 1)
            degrees, top = parse_degrees(degree_text)
            result[key] = dict(valid=parts[3] == "true" and flag == "true", status="checked",
                               degrees=degrees, max_degree=top)
    return result


def write_report(data, dest):
    variants = list(data["variants"])
    rows = [r for r in data["runs"] if r["method"] == "StrictlyPositiveCert"]
    report = f"""# lower_bound_poly and FACTOR_SOS_PRE vs main on tests/benchmarking.mpl

All variants ran the working tree's `tests/benchmarking.mpl` and helpers unchanged
(active suites: {", ".join(f"`{s}`" for s in data["suites"])}). `main` is `{data['main_ref']}`; its
`src/` was extracted into a scratch tree next to a copy of the working tree's `tests/`.
Fixed Makefile flags for every variant:

```text
{data['make_command']}
```

Variants: {"; ".join(f"`{v}` adds `{' '.join(f) or 'nothing'}`" for v, f in data["variants"].items())}.

{REPEATS} fresh-process runs per variant, rotating order. Timings are the CPU seconds
printed by `checkWeifeng` (StrictlyPositiveCert only); aggregates are sums of per-case
medians over cases that succeeded in every run of every variant. Degrees are **exact**:
each certificate was recomputed and expanded by `tests/verify_external_certificates.mpl`,
which also checks `f = sigma_0 + sum(sigma_i*g_i)` and the nonnegativity of every
sigma (columns "valid"). The benchmark's own quickdegree misses cancellation.
"""
    for suite in data["suites"]:
        srows = [r for r in rows if r["suite"] == suite]
        cases = sorted({r["case"] for r in srows})
        def values(case, v):
            return [r["cpu_seconds"] for r in srows
                    if r["case"] == case and r["variant"] == v and r["status"] == "success"]
        common = [c for c in cases if all(len(values(c, v)) == REPEATS for v in variants)]
        totals = {v: sum(statistics.median(values(c, v)) for c in common) for v in variants}
        ver = {v: {c: data["verification"][v].get(f"{suite}:{c}") for c in cases} for v in variants}
        report += f"\n## Suite `{suite}` ({len(cases)} cases)\n\n"
        report += "| Variant | CPU time (s) | Change vs main | Successful calls | Valid certificates | Max degree vs main (lower / equal / higher) | Max-degree histogram |\n"
        report += "|---|---:|---:|---:|---:|---|---|\n"
        for v in variants:
            ok = sum(r["status"] == "success" for r in srows if r["variant"] == v)
            valid = sum(1 for c in cases if ver[v][c] and ver[v][c]["valid"])
            cmp = Counter()
            for c in cases:
                a, b = ver["main"][c], ver[v][c]
                if not (a and b and a["degrees"] and b["degrees"]):
                    cmp["n/a"] += 1
                else:
                    cmp["lower" if b["max_degree"] < a["max_degree"] else
                        "higher" if b["max_degree"] > a["max_degree"] else "equal"] += 1
            hist = Counter(ver[v][c]["max_degree"] for c in cases if ver[v][c] and ver[v][c]["degrees"])
            delta = 100*(totals[v]/totals["main"]-1) if totals["main"] else 0
            report += (f"| {v} | {totals[v]:.3f} | {delta:+.1f}% | {ok}/{len(cases)*REPEATS} | "
                       f"{valid}/{len(cases)} | "
                       + ("—" if v == "main" else f"{cmp['lower']} / {cmp['equal']} / {cmp['higher']}"
                          + (f" ({cmp['n/a']} n/a)" if cmp["n/a"] else ""))
                       + " | " + ", ".join(f"{d}: {n}" for d, n in sorted(hist.items())) + " |\n")
        report += f"\n<details><summary>Per-case median CPU seconds; exact degrees</summary>\n\n"
        report += "| Case | " + " | ".join(variants) + " |\n|---:|" + "---|"*len(variants) + "\n"
        for c in cases:
            cells = []
            for v in variants:
                vals, info = values(c, v), ver[v][c]
                deg = info["degrees"] if info and info["degrees"] else "error"
                mark = "" if info and info["valid"] else " INVALID"
                cells.append((f"{statistics.median(vals):.3f}" if vals else "—") + f"; {deg}{mark}")
            report += f"| {c} | " + " | ".join(cells) + " |\n"
        report += "\n</details>\n"
    report += "\n## Status counts (all methods)\n\n| Suite | Variant | Method | Status | Calls |\n|---|---|---|---|---:|\n"
    for key, n in sorted(Counter((r["suite"], r["variant"], r["method"], r["status"])
                                 for r in data["runs"]).items()):
        report += f"| {' | '.join(key)} | {n} |\n"
    report += f"\nRaw data: `lower_bound_comparison_runs.json`. Logs: `{data['log_directory']}`.\n"
    (dest / "lower_bound_comparison.md").write_text(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-ref", default="main")
    args = parser.parse_args()
    build = subprocess.check_output(["make", "-Bn"], cwd=ROOT, text=True).strip()
    command = shlex.split(build)
    if len(build.splitlines()) != 1 or command[0] != "maple":
        raise RuntimeError("Unexpected Makefile command")
    fixed = command[1:-1]
    suites = active_suites()
    ref = subprocess.check_output(["git", "rev-parse", "--short", args.main_ref],
                                  cwd=ROOT, text=True).strip()
    variant_flags = {"main": [], "working_tree": [], "working_tree+FACTOR_SOS_PRE": ["-DFACTOR_SOS_PRE"]}
    roots = {"main": main_tree(ref), "working_tree": ROOT, "working_tree+FACTOR_SOS_PRE": ROOT}
    variants = list(variant_flags)
    before = hashes()
    logdir = Path(tempfile.mkdtemp(prefix="lbp_comparison_"))
    dest = ROOT / "tests/results"
    dest.mkdir(exist_ok=True)
    runs, verification = [], {}
    data = dict(make_command=build, fixed_flags=fixed, main_ref=ref, suites=suites,
                variants=variant_flags, repetitions=REPEATS, source_sha256=before,
                log_directory=str(logdir), runs=runs, verification=verification)
    expected = None
    for repeat in range(1, REPEATS+1):
        shift = (repeat-1) % len(variants)
        for variant in variants[shift:] + variants[:shift]:
            log = logdir / f"{variant}_{repeat}.log"
            started = time.monotonic()
            with log.open("w") as stream:
                subprocess.run(["maple", "-q", *fixed, *variant_flags[variant],
                                "tests/run_factor_sos_benchmarking.mpl"],
                               cwd=roots[variant], stdout=stream, stderr=subprocess.STDOUT, check=True)
            parsed = parse(log.read_text(), suites, variant, repeat)
            signature = [(r["suite"], r["case"], r["method"], r["polynomial"]) for r in parsed]
            if expected is None:
                expected = signature
            if signature != expected:
                raise RuntimeError("Benchmark inputs differ between runs")
            runs.extend(parsed)
            (dest / "lower_bound_comparison_runs.json").write_text(json.dumps(data, indent=2)+"\n")
            print(f"{variant} run {repeat}/{REPEATS}: {len(parsed)//3} cases in "
                  f"{time.monotonic()-started:.1f}s", flush=True)
            if hashes() != before:
                raise RuntimeError("Sources/data changed during the benchmark")
    for variant in variants:
        result = verify(roots[variant], fixed + variant_flags[variant], suites, logdir, variant)
        verification[variant] = {f"{s}:{c}": v for (s, c), v in result.items()}
        print(f"{variant} verification: {sum(v['valid'] for v in result.values())}/{len(result)} valid",
              flush=True)
    (dest / "lower_bound_comparison_runs.json").write_text(json.dumps(data, indent=2)+"\n")
    write_report(data, dest)


if __name__ == "__main__":
    main()
