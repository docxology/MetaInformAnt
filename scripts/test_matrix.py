#!/usr/bin/env python3
"""Run independent test groups in parallel.

The groups are deliberately process-isolated: tests that change CWD, import
optional dependencies, or use module-level state cannot contaminate another
group.  This is a fast developer/CI orchestration layer; the canonical full
suite remains ``uv run pytest -q``.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class TestGroup:
    """A non-overlapping pytest selection."""

    name: str
    paths: tuple[str, ...]
    description: str


GROUPS: tuple[TestGroup, ...] = (
    TestGroup("core", ("tests/core",), "core I/O, validation, execution, and compatibility contracts"),
    TestGroup("dna", ("tests/dna",), "DNA sequence, alignment, and external-access contracts"),
    TestGroup("rna", ("tests/rna",), "RNA orchestration, provenance, retry, and readiness contracts"),
    TestGroup("mcp", ("tests/mcp",), "MCP-adjacent monitor and capability contracts"),
    TestGroup(
        "domains",
        ("tests/ecology", "tests/gwas", "tests/information", "tests/visualization"),
        "domain analysis contracts",
    ),
    TestGroup(
        "package",
        ("tests/infrastructure", "tests/quality", "tests/other"),
        "packaging, import, and implementation-policy contracts",
    ),
)


def _group_map() -> dict[str, TestGroup]:
    return {group.name: group for group in GROUPS}


def _run_group(group: TestGroup, *, pytest_args: tuple[str, ...]) -> tuple[str, int, str]:
    command = [sys.executable, "-m", "pytest", *pytest_args, *group.paths]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + environment.get("PYTHONPATH", "")
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    output = completed.stdout + ("\n" + completed.stderr if completed.stderr else "")
    return group.name, completed.returncode, output.strip()


def main(argv: list[str] | None = None) -> int:
    """Parse options, run selected groups, and return a shell-friendly status."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List available groups and exit")
    parser.add_argument("--group", action="append", dest="groups", help="Group to run; repeatable (default: all)")
    parser.add_argument("--workers", type=int, default=0, help="Maximum concurrent groups (default: one per group)")
    parser.add_argument("--slow", action="store_true", help="Show pytest durations")
    parser.add_argument(
        "--include-network",
        action="store_true",
        help="Include tests marked network or external_tool (default: skip optional lanes)",
    )
    args = parser.parse_args(argv)

    if args.list:
        for group in GROUPS:
            print(f"{group.name}: {group.description} [{', '.join(group.paths)}]")
        return 0

    available = _group_map()
    names = args.groups or [group.name for group in GROUPS]
    unknown = sorted(set(names) - available.keys())
    if unknown:
        parser.error(f"unknown group(s): {', '.join(unknown)}")

    selected = [available[name] for name in names]
    worker_count = args.workers or len(selected)
    if worker_count < 1:
        parser.error("--workers must be positive")
    pytest_args = ("-q", "--durations=10") if args.slow else ("-q",)
    if not args.include_network:
        pytest_args += ("-m", "not network and not external_tool")

    failures = 0
    with ThreadPoolExecutor(max_workers=min(worker_count, len(selected))) as executor:
        futures = [executor.submit(_run_group, group, pytest_args=pytest_args) for group in selected]
        for future in as_completed(futures):
            name, returncode, output = future.result()
            print(f"\n=== {name} ({'PASS' if returncode == 0 else 'FAIL'}) ===")
            if output:
                print(output)
            failures += returncode != 0

    print(f"\n{len(selected) - failures}/{len(selected)} test groups passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
