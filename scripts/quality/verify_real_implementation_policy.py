#!/usr/bin/env python3
"""Verify real-implementation policy references and test-double guardrails."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]

SKIP_DIRS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".uv-cache",
    ".venv",
    "__pycache__",
    "build",
    "output",
}

SKIP_PATHS = {
    Path("docs/REAL_IMPLEMENTATION_POLICY.md"),
    Path("tests/REAL_IMPLEMENTATION_TESTING_POLICY.md"),
    Path("scripts/quality/verify_real_implementation_policy.py"),
    Path("tests/quality/test_real_implementation_policy.py"),
}

TEXT_SUFFIXES = {
    ".md",
    ".py",
    ".toml",
    ".yaml",
    ".yml",
    ".json",
    ".sh",
    ".txt",
    ".cursorrules",
}

OLD_POLICY_PATTERNS = [
    re.compile(r"NO_MOCKING_POLICY"),
    re.compile(r"\bNO_MOCKING\b"),
    re.compile(r"\bno_mock\b"),
    re.compile(r"\bZero[- ]Mock(?:ing)?\b", re.IGNORECASE),
    re.compile(r"\bno[- ]mock(?:ing)?\b", re.IGNORECASE),
]

TEST_DOUBLE_PATTERNS = [
    re.compile(r"\bunittest\.mock\b"),
    re.compile(r"\bfrom\s+unittest\s+import\s+mock\b"),
    re.compile(r"(?<![A-Za-z0-9_])import\s+mock\b"),
    re.compile(r"\bpytest-mock\b"),
    re.compile(r"\bMagicMock\b"),
    re.compile(r"(?<![A-Za-z0-9_])Mock\("),
    re.compile(r"@\s*patch\b"),
    re.compile(r"(?<![A-Za-z0-9_])patch\("),
    re.compile(r"\bmocker\.patch\b"),
]


@dataclass(frozen=True)
class PolicyViolation:
    """One real-implementation policy violation."""

    path: Path
    line: int
    rule: str
    text: str


def iter_policy_files(repo_root: Path = REPO_ROOT) -> Iterable[Path]:
    """Yield text files covered by the policy scanner."""
    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(repo_root)
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        if rel in SKIP_PATHS:
            continue
        if path.suffix in TEXT_SUFFIXES or path.name in {"AGENTS.md", "SPEC.md", "PAI.md"}:
            yield path


def _scan_patterns(
    path: Path,
    patterns: Sequence[re.Pattern[str]],
    rule: str,
    repo_root: Path,
) -> list[PolicyViolation]:
    violations: list[PolicyViolation] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        return violations

    for line_number, line in enumerate(lines, start=1):
        for pattern in patterns:
            if pattern.search(line):
                violations.append(
                    PolicyViolation(
                        path=path.relative_to(repo_root),
                        line=line_number,
                        rule=rule,
                        text=line.strip()[:180],
                    )
                )
                break
    return violations


def scan_repo(repo_root: Path = REPO_ROOT) -> list[PolicyViolation]:
    """Return all current policy violations."""
    violations: list[PolicyViolation] = []
    repo_root = repo_root.resolve()
    for path in iter_policy_files(repo_root):
        violations.extend(_scan_patterns(path, OLD_POLICY_PATTERNS, "old-policy-reference", repo_root))
        violations.extend(_scan_patterns(path, TEST_DOUBLE_PATTERNS, "test-double-api", repo_root))
    return violations


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description="Verify real-implementation policy compliance")
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args(argv)

    violations = scan_repo(args.repo_root.resolve())
    if not violations:
        print("Real-implementation policy scan passed.")
        return 0

    print(f"Real-implementation policy scan found {len(violations)} violation(s):", file=sys.stderr)
    for violation in violations[:200]:
        print(
            f"{violation.path}:{violation.line}: {violation.rule}: {violation.text}",
            file=sys.stderr,
        )
    if len(violations) > 200:
        print(f"... {len(violations) - 200} more violation(s) omitted", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
