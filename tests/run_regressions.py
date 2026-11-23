"""Run source-based Maple regressions and detect errors even if Maple exits zero."""
from pathlib import Path
import os
import re
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parent.parent
scripts = (
    "regression_monogenic.mpl",
    "regression_epsilon.mpl",
    "regression_positivity.mpl",
    "regression_signs.mpl",
    "regression_scheduler.mpl",
    "regression_quick_certificates.mpl",
    "regression_logging.mpl",
    "regression_univariate_heuristics.mpl",
)
for name in scripts:
    flags = list(sys.argv[1:])
    if name == "regression_quick_certificates.mpl":
        flags.append("-DENABLE_QUICK_CERT")
    if name == "regression_logging.mpl":
        flags.extend(("-DLOG_TIME", "-DENABLE_QUICK_CERT"))
    with tempfile.TemporaryDirectory(prefix="spcert-regression-") as log_dir:
        result = subprocess.run(
            ["maple", "-q", *flags, str(root / "tests" / name)],
            cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=300,
            env={**os.environ, "SPCERT_LOG_DIR": log_dir},
        )
        print(result.stdout, end="")
        if result.returncode or "PASS:" not in result.stdout or re.search(
            r"(?m)^Error(?:,|:|\s)", result.stdout
        ):
            sys.exit(f"FAILED: {name}")
        if name == "regression_logging.mpl":
            stack = []
            expect_time = False
            entries = 0
            for line in (Path(log_dir) / "log_time.txt").read_text().splitlines():
                if expect_time:
                    assert line.startswith("Time: ") and float(line[6:]) >= 0, line
                    expect_time = False
                    continue
                event = re.fullmatch(r"(Start|End): (\S+) (\d+)", line)
                assert event, line
                kind, routine, depth = event.groups()
                if kind == "Start":
                    assert int(depth) == len(stack), line
                    stack.append(routine)
                    entries += 1
                else:
                    assert stack and stack.pop() == routine, line
                    assert int(depth) == len(stack), line
                    expect_time = True
            assert entries and not stack and not expect_time, "incomplete timing log"
