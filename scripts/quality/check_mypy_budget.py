#!/usr/bin/env python3
"""Run mypy and fail when the documented error budget increases."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=Path, default=Path("config/quality/mypy_error_budget.txt"))
    parser.add_argument("paths", nargs="*", default=["src/metainformant"])
    args = parser.parse_args()

    budget = int(args.budget.read_text(encoding="utf-8").strip())
    result = subprocess.run(
        ["mypy", "--show-error-codes", *args.paths],
        text=True,
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    errors = len(re.findall(r": error(?:\[[^]]+\])?:", output))
    print(output, end="")
    print(f"mypy error budget: observed={errors} allowed={budget}")
    if result.returncode not in (0, 1):
        return result.returncode
    if errors > budget:
        print("mypy error budget increased; fix or explicitly lower the baseline with a review.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
