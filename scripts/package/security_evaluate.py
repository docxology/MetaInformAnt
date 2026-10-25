#!/usr/bin/env python3
"""Evaluate security-scan reports against the accepted-findings baseline.

``security_checks.sh`` runs the scanners and calls this module with the names
of the scanners it ran. Each named scanner must leave a complete report: a
missing, unreadable, empty or partial scan is a blocking failure, never a clean
result. Exit status 0 means no blocking findings, 1 means blocking findings or
missing evidence, and 2 means the invocation itself was invalid.

Baseline format (``security_baseline.txt``): ``SEVERITY|OWNER|ID|LOCATION``.
Bandit and pattern LOCATIONs match as path substrings; pip-audit LOCATIONs
must equal the package name.
"""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

SCANNERS = ("bandit", "deps", "patterns")


@dataclass(frozen=True)
class BaselineEntry:
    """One accepted finding."""

    severity: str
    owner: str
    finding_id: str
    location: str


@dataclass
class Evaluation:
    """Blocking findings and informational counts from one invocation."""

    blocking: list[str] = field(default_factory=list)
    review: list[str] = field(default_factory=list)
    bandit_low: int = 0
    audited: int = 0


def load_baseline(path: Path) -> list[BaselineEntry]:
    """Parse the baseline; malformed lines raise ``ValueError``."""
    entries = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = [part.strip() for part in line.split("|")]
        if len(fields) != 4:
            raise ValueError(f"Malformed baseline line (need SEVERITY|OWNER|ID|LOCATION): {line!r}")
        entries.append(BaselineEntry(*fields))
    return entries


def is_baselined(
    baseline: Sequence[BaselineEntry], severity: str, finding_id: str, location: str, *, exact: bool = False
) -> bool:
    """Return True when an entry accepts this finding."""
    for entry in baseline:
        if entry.severity not in ("*", severity) or entry.finding_id not in ("-", "*", finding_id):
            continue
        if entry.location in ("-", "*"):
            return True
        if (entry.location == location) if exact else (entry.location in location):
            return True
    return False


def _json_object(path: Path, scanner: str, result: Evaluation) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        result.blocking.append(f"{scanner} produced no readable report ({exc.__class__.__name__}); the scan failed")
        return None
    if not isinstance(data, dict):
        result.blocking.append(f"{scanner} report is not a JSON object; the scan failed")
        return None
    return data


def _lines(path: Path, scanner: str, result: Evaluation) -> list[str]:
    try:
        return [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    except OSError as exc:
        result.blocking.append(f"{scanner} produced no readable report ({exc.__class__.__name__}); the scan failed")
        return []


def evaluate_bandit(path: Path, baseline: Sequence[BaselineEntry], result: Evaluation) -> None:
    """Block on medium/high findings and on any sign the scan did not cover the tree."""
    report = _json_object(path, "bandit", result)
    if report is None:
        return
    results, errors = report.get("results"), report.get("errors")
    if not isinstance(results, list) or not isinstance(errors, list):
        result.blocking.append("bandit report lacks its results/errors lists; the scan failed")
        return
    # Bandit exits 0 when it cannot read the scan root or parse a file; those
    # failures appear only here and in a zero line count.
    for error in errors:
        detail = error if not isinstance(error, dict) else f"{error.get('filename', '?')}: {error.get('reason', '?')}"
        result.blocking.append(f"bandit could not scan {detail}")
    totals = report.get("metrics", {}).get("_totals", {}) if isinstance(report.get("metrics"), dict) else {}
    loc = totals.get("loc") if isinstance(totals, dict) else None
    if not isinstance(loc, int) or loc <= 0:
        result.blocking.append("bandit scanned no source lines; the scan root is missing or empty")
    for finding in results:
        severity = finding["issue_severity"].lower()
        if severity not in ("medium", "high"):
            result.bandit_low += 1
            continue
        if not is_baselined(baseline, severity, finding["test_id"], finding["filename"]):
            result.blocking.append(
                f"bandit {severity.upper()} {finding['test_id']} {finding['test_name']} "
                f"at {finding['filename']}:{finding['line_number']}"
            )


def evaluate_pip_audit(path: Path, pylock: Path, baseline: Sequence[BaselineEntry], result: Evaluation) -> None:
    """Block on vulnerabilities, unaudited packages and incomplete lock coverage."""
    report = _json_object(path, "pip-audit", result)
    if report is None:
        return
    dependencies = report.get("dependencies")
    if not isinstance(dependencies, list) or not dependencies:
        result.blocking.append("pip-audit audited no dependencies; the lock export was empty")
        return
    try:
        locked = len(tomllib.loads(pylock.read_text(encoding="utf-8")).get("packages", []))
    except (OSError, ValueError) as exc:
        result.blocking.append(f"exported pylock.toml is unreadable ({exc.__class__.__name__}); coverage unknown")
        locked = None
    if locked is not None and len(dependencies) != locked:
        result.blocking.append(
            f"pip-audit reported {len(dependencies)} of {locked} locked packages; coverage incomplete"
        )
    for dep in dependencies:
        name = dep.get("name", "?")
        if dep.get("skip_reason"):
            if not is_baselined(baseline, "vuln", "unaudited", name, exact=True):
                result.blocking.append(f"pip-audit could not audit {name}: {dep['skip_reason']}")
            continue
        result.audited += 1
        for vuln in dep.get("vulns", []):
            vuln_id = vuln.get("id", "unknown")
            if not is_baselined(baseline, "vuln", vuln_id, name, exact=True):
                fixes = ", ".join(vuln.get("fix_versions") or ["none"])
                result.blocking.append(f"pip-audit vuln {vuln_id} on {name}=={dep.get('version', '?')} (fix: {fixes})")


def evaluate_patterns(
    pattern_path: Path, traversal_path: Path, baseline: Sequence[BaselineEntry], result: Evaluation
) -> None:
    """Block on shell-execution hits; path-traversal candidates are review-only."""
    for hit in _lines(pattern_path, "pattern audit", result):
        if not is_baselined(baseline, "pattern", "PATTERN", hit.split(":", 1)[0]):
            result.blocking.append(f"pattern audit hit: {hit}")
    result.review.extend(_lines(traversal_path, "path-traversal audit", result))


def main(argv: Sequence[str] | None = None) -> int:
    """Evaluate the named scanners' reports; see the module docstring for exit codes."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ran", required=True, help="comma-separated scanners run by this invocation")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--bandit", type=Path, required=True)
    parser.add_argument("--pip-audit", type=Path, required=True)
    parser.add_argument("--pylock", type=Path, required=True)
    parser.add_argument("--patterns", type=Path, required=True)
    parser.add_argument("--traversal", type=Path, required=True)
    args = parser.parse_args(argv)

    ran = [name for name in args.ran.split(",") if name]
    if not ran or set(ran) - set(SCANNERS):
        print(f"ERROR: scanners must be a non-empty subset of {', '.join(SCANNERS)}; got {args.ran!r}")
        return 2
    try:
        baseline = load_baseline(args.baseline)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1

    result = Evaluation()
    if "bandit" in ran:
        evaluate_bandit(args.bandit, baseline, result)
    if "deps" in ran:
        evaluate_pip_audit(args.pip_audit, args.pylock, baseline, result)
    if "patterns" in ran:
        evaluate_patterns(args.patterns, args.traversal, baseline, result)

    for hit in result.review:
        print(f"REVIEW (non-blocking): path-traversal candidate: {hit}")
    print(
        f"Scanners evaluated: {', '.join(sorted(ran))}; baseline entries: {len(baseline)}; "
        f"locked packages audited: {result.audited}; "
        f"bandit low-severity (non-blocking): {result.bandit_low}; "
        f"path-traversal review candidates (non-blocking): {len(result.review)}"
    )
    if result.blocking:
        print(f"BLOCKING: {len(result.blocking)} security finding(s) not covered by the baseline:")
        for finding in result.blocking:
            print(f"  - {finding}")
        return 1
    print("No blocking security findings.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
