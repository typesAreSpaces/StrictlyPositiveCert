"""Run source-based Maple regressions and detect errors even if Maple exits zero."""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parent.parent
scripts = (
    "regression_monogenic.mpl",
    "regression_lower_bound.mpl",
    "regression_quick_certificates.mpl",
)
for name in scripts:
    flags = list(sys.argv[1:])
    if name == "regression_quick_certificates.mpl":
        flags.append("-DENABLE_QUICK_CERT")
    result = subprocess.run(
        ["maple", "-q", *flags, str(root / "tests" / name)],
        cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=300,
    )
    print(result.stdout, end="")
    if result.returncode or "PASS:" not in result.stdout or re.search(
        r"(?m)^Error(?:,|:|\s)", result.stdout
    ):
        sys.exit(f"FAILED: {name}")
