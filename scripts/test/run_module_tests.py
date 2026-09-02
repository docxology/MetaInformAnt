#!/usr/bin/env python3
"""Per-module test runner for METAINFORMANT.

Enforces the one-pytest-directory-per-invocation rule: each module's tests
run in their own subprocess pytest invocation, so per-directory ``conftest``
plugins can never collide.

Exit codes (CI-consumable):
    0 - all requested modules passed
    1 - one or more modules failed
    2 - usage/configuration error (unknown module, invalid arguments)

Examples:
    uv run python scripts/test/run_module_tests.py --module rna
    uv run python scripts/test/run_module_tests.py --module rna --coverage
    uv run python scripts/test/run_module_tests.py --module rna --module dna --json results.json
    uv run python scripts/test/run_module_tests.py --list
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DIR = REPO_ROOT / "tests"

# Directories under tests/ that are not per-source-module suites.
NON_MODULE_DIRS = {
    "__pycache__",
    "_support",
    "data",
    "fixtures",
    "infrastructure",
    "integration",
    "other",
    "scripts",
    "cstmm_fixture",
}

DEFAULT_TIMEOUT_SECONDS = 1800


@dataclass
class ModuleResult:
    """Outcome of one module's pytest run."""

    module: str
    test_dir: str
    returncode: int
    duration_seconds: float
    timed_out: bool = False
    collected: Optional[int] = None
    passed: Optional[int] = None
    failed: Optional[int] = None
    skipped: Optional[int] = None
    errors: Optional[int] = None
    command: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.timed_out

    def to_dict(self) -> dict:
        return {
            "module": self.module,
            "test_dir": self.test_dir,
            "ok": self.ok,
            "returncode": self.returncode,
            "timed_out": self.timed_out,
            "duration_seconds": round(self.duration_seconds, 2),
            "collected": self.collected,
            "passed": self.passed,
            "failed": self.failed,
            "skipped": self.skipped,
            "errors": self.errors,
            "command": self.command,
        }


def discover_modules() -> List[str]:
    """List per-module test directories (directories mirroring src modules)."""
    modules = []
    for entry in sorted(TESTS_DIR.iterdir()):
        if entry.is_dir() and entry.name not in NON_MODULE_DIRS:
            if any(entry.glob("test_*.py")):
                modules.append(entry.name)
    return modules


def parse_summary(output: str) -> dict:
    """Extract the pytest summary line counts from captured output."""
    import re

    counts = {"collected": None, "passed": None, "failed": None, "skipped": None, "errors": None}
    for line in output.splitlines():
        if "collected" in line and " items" in line:
            match = re.search(r"collected (\d+) items", line)
            if match:
                counts["collected"] = int(match.group(1))
        if line.startswith("=") and ("passed" in line or "failed" in line or "skipped" in line):
            for key in ("passed", "failed", "skipped", "errors"):
                match = re.search(rf"(\d+) {key}", line)
                if match:
                    counts[key] = int(match.group(1))
    return counts


def run_module_tests(
    module: str,
    coverage: bool = False,
    timeout: int = DEFAULT_TIMEOUT_SECONDS,
    extra_args: Optional[List[str]] = None,
    verbose: bool = False,
) -> ModuleResult:
    """Run one module's test directory in a dedicated pytest subprocess."""
    test_dir = TESTS_DIR / module
    if not test_dir.is_dir():
        raise SystemExit(f"Unknown module '{module}': no tests/{module}/ directory")

    cmd = [sys.executable, "-m", "pytest", str(test_dir), "-p", "no:cacheprovider"]
    if verbose:
        cmd.append("-v")
    if coverage:
        cmd += ["--cov", f"src/metainformant/{module}", "--cov-report=term-missing", "--cov-append", "--no-cov-on-fail"]
    if extra_args:
        cmd += extra_args

    start = time.monotonic()
    timed_out = False
    try:
        proc = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        returncode, stdout = proc.returncode, proc.stdout + proc.stderr[-4000:]
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        returncode = 124
        out = exc.stdout or b""
        err = exc.stderr or b""
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        if isinstance(err, bytes):
            err = err.decode(errors="replace")
        stdout = str(out) + str(err)

    duration = time.monotonic() - start
    summary = parse_summary(stdout)
    return ModuleResult(
        module=module,
        test_dir=str(test_dir.relative_to(REPO_ROOT)),
        returncode=returncode,
        duration_seconds=duration,
        timed_out=timed_out,
        collected=summary["collected"],
        passed=summary["passed"],
        failed=summary["failed"],
        skipped=summary["skipped"],
        errors=summary["errors"],
        command=cmd,
    )


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Per-module METAINFORMANT test runner")
    parser.add_argument("--module", action="append", default=[], help="module name under tests/ (repeatable)")
    parser.add_argument("--all", action="store_true", help="run every discovered module, one pytest per module")
    parser.add_argument("--coverage", action="store_true", help="per-module coverage of src/metainformant/<module>")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS, help="per-module timeout in seconds")
    parser.add_argument("--json", dest="json_path", default=None, help="write JSON results to this path")
    parser.add_argument("--list", action="store_true", help="list discovered modules and exit")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose pytest output")
    args = parser.parse_args(argv)

    modules = discover_modules()
    if args.list:
        print("\n".join(modules))
        return 0

    if args.all:
        selected = modules
    elif args.module:
        known = set(modules)
        unknown = [m for m in args.module if m not in known]
        if unknown:
            print(f"Unknown module(s): {', '.join(unknown)}", file=sys.stderr)
            print(f"Known modules: {', '.join(modules)}", file=sys.stderr)
            return 2
        selected = args.module
    else:
        parser.print_usage()
        print("\nNothing to do: pass --module <name>, --all, or --list.", file=sys.stderr)
        return 2

    results = []
    for module in selected:
        print(f"[run_module_tests] pytest tests/{module}/ ...", flush=True)
        result = run_module_tests(
            module,
            coverage=args.coverage,
            timeout=args.timeout,
            verbose=args.verbose,
        )
        results.append(result)
        status = "OK" if result.ok else ("TIMEOUT" if result.timed_out else f"FAIL (rc={result.returncode})")
        print(
            f"[run_module_tests] {module}: {status} in {result.duration_seconds:.1f}s "
            f"(passed={result.passed} failed={result.failed} skipped={result.skipped} errors={result.errors})",
            flush=True,
        )

    if args.json_path:
        payload = {
            "modules": [r.to_dict() for r in results],
            "all_ok": all(r.ok for r in results),
        }
        Path(args.json_path).write_text(json.dumps(payload, indent=2) + "\n")

    return 0 if all(r.ok for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
