"""Table-driven checks for the blocking security-report evaluator.

Zero-mocks: every case writes real report files and runs the evaluator's
command-line entry point exactly as ``security_checks.sh`` invokes it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from scripts.package import security_evaluate as gate

MISSING = object()
CLEAN_BANDIT: dict[str, Any] = {"errors": [], "results": [], "metrics": {"_totals": {"loc": 120}}}
PYLOCK = (
    'lock-version = "1.0"\ncreated-by = "fixture"\n'
    '[[packages]]\nname = "alpha"\nversion = "1.0"\n'
    '[[packages]]\nname = "beta"\nversion = "2.0"\n'
)
CLEAN_AUDIT: dict[str, Any] = {
    "dependencies": [
        {"name": "alpha", "version": "1.0", "vulns": []},
        {"name": "beta", "version": "2.0", "vulns": []},
    ],
    "fixes": [],
}


def _finding(severity: str) -> dict[str, Any]:
    return {
        "issue_severity": severity,
        "test_id": "B602",
        "test_name": "subprocess_popen_with_shell_equals_true",
        "filename": "src/metainformant/core/run.py",
        "line_number": 7,
    }


def _vulnerable(vuln_id: str = "GHSA-aaaa-bbbb-cccc") -> dict[str, Any]:
    alpha = {"name": "alpha", "version": "1.0", "vulns": [{"id": vuln_id, "fix_versions": ["1.1"]}]}
    return {"dependencies": [alpha, CLEAN_AUDIT["dependencies"][1]], "fixes": []}


def _write(path: Path, value: object) -> None:
    if value is MISSING:
        return
    if isinstance(value, bytes):
        path.write_bytes(value)
    elif isinstance(value, str):
        path.write_text(value, encoding="utf-8")
    else:
        path.write_text(json.dumps(value), encoding="utf-8")


def run_gate(
    tmp_path: Path,
    *,
    ran: str = "bandit,deps,patterns",
    bandit: object = CLEAN_BANDIT,
    audit: object = CLEAN_AUDIT,
    pylock: object = PYLOCK,
    patterns: object = "",
    traversal: object = "",
    baseline: str = "",
) -> int:
    paths = {name: tmp_path / name for name in ("bandit.json", "audit.json", "pylock.toml", "patterns", "traversal")}
    for name, value in zip(paths, (bandit, audit, pylock, patterns, traversal)):
        _write(paths[name], value)
    (tmp_path / "baseline.txt").write_text(baseline, encoding="utf-8")
    return gate.main(
        [
            "--ran",
            ran,
            "--baseline",
            str(tmp_path / "baseline.txt"),
            "--bandit",
            str(paths["bandit.json"]),
            "--pip-audit",
            str(paths["audit.json"]),
            "--pylock",
            str(paths["pylock.toml"]),
            "--patterns",
            str(paths["patterns"]),
            "--traversal",
            str(paths["traversal"]),
        ]
    )


CASES = {
    "clean": ({}, 0),
    # Bandit exits 0 on these; only the report shows the scan did not happen.
    "bandit-missing-root": (
        {"bandit": {"errors": [{"filename": "src", "reason": "No such file or directory"}], "results": []}},
        1,
    ),
    "bandit-parse-error": (
        {
            "bandit": {
                "errors": [{"filename": "src/x.py", "reason": "syntax error while parsing AST from file"}],
                "results": [],
                "metrics": {"_totals": {"loc": 120}},
            }
        },
        1,
    ),
    "bandit-zero-lines": ({"bandit": {"errors": [], "results": [], "metrics": {"_totals": {"loc": 0}}}}, 1),
    "bandit-empty-object": ({"bandit": {}}, 1),
    "bandit-null": ({"bandit": None}, 1),
    "bandit-zero-bytes": ({"bandit": b""}, 1),
    "bandit-truncated": ({"bandit": '{"errors": [], "res'}, 1),
    "bandit-missing-file": ({"bandit": MISSING}, 1),
    "bandit-medium": ({"bandit": {**CLEAN_BANDIT, "results": [_finding("MEDIUM")]}}, 1),
    "bandit-medium-baselined": (
        {
            "bandit": {**CLEAN_BANDIT, "results": [_finding("MEDIUM")]},
            "baseline": "medium|owners|B602|src/metainformant/core/run.py\n",
        },
        0,
    ),
    "bandit-low-is-review-only": ({"bandit": {**CLEAN_BANDIT, "results": [_finding("LOW")]}}, 0),
    "audit-null": ({"audit": None}, 1),
    "audit-empty-object": ({"audit": {}}, 1),
    "audit-no-dependencies": ({"audit": {"dependencies": [], "fixes": []}}, 1),
    "audit-zero-bytes": ({"audit": b""}, 1),
    "audit-missing-file": ({"audit": MISSING}, 1),
    "audit-partial-coverage": ({"audit": {"dependencies": CLEAN_AUDIT["dependencies"][:1], "fixes": []}}, 1),
    "audit-unreadable-lock": ({"pylock": MISSING}, 1),
    "audit-vulnerable": ({"audit": _vulnerable()}, 1),
    "audit-vulnerable-baselined": (
        {"audit": _vulnerable(), "baseline": "vuln|owners|GHSA-aaaa-bbbb-cccc|alpha\n"},
        0,
    ),
    # pip-audit locations name a package exactly; a prefix must not match.
    "audit-baseline-prefix-does-not-match": (
        {"audit": _vulnerable(), "baseline": "vuln|owners|GHSA-aaaa-bbbb-cccc|alp\n"},
        1,
    ),
    "audit-skipped": (
        {
            "audit": {
                "dependencies": [
                    {"name": "alpha", "skip_reason": "Dependency not found on PyPI"},
                    CLEAN_AUDIT["dependencies"][1],
                ],
                "fixes": [],
            }
        },
        1,
    ),
    "audit-skipped-baselined": (
        {
            "audit": {
                "dependencies": [
                    {"name": "alpha", "skip_reason": "Dependency not found on PyPI"},
                    CLEAN_AUDIT["dependencies"][1],
                ],
                "fixes": [],
            },
            "baseline": "vuln|owners|unaudited|alpha\n",
        },
        0,
    ),
    "pattern-hit": ({"patterns": "src/metainformant/core/run.py:7:os.system(cmd)\n"}, 1),
    "pattern-hit-baselined": (
        {
            "patterns": "src/metainformant/core/run.py:7:os.system(cmd)\n",
            "baseline": "pattern|owners|PATTERN|src/metainformant/core/run.py\n",
        },
        0,
    ),
    "pattern-report-missing": ({"patterns": MISSING}, 1),
    "traversal-is-review-only": ({"traversal": 'src/metainformant/core/paths.py:5:    "../"\n'}, 0),
    # Only the scanners named in --ran are evaluated.
    "unrun-scanner-reports-are-ignored": ({"ran": "bandit", "audit": MISSING, "patterns": MISSING}, 0),
    "empty-scanner-list": ({"ran": ""}, 2),
    "unknown-scanner": ({"ran": "bandit,typo"}, 2),
    "malformed-baseline": ({"baseline": "medium|owners|B602\n"}, 1),
}


@pytest.mark.parametrize(("overrides", "expected"), list(CASES.values()), ids=list(CASES))
def test_evaluator_blocks_on_findings_and_missing_evidence(
    tmp_path: Path, overrides: dict[str, Any], expected: int
) -> None:
    assert run_gate(tmp_path, **overrides) == expected


def test_clean_run_reports_every_audited_package(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert run_gate(tmp_path) == 0
    assert "locked packages audited: 2" in capsys.readouterr().out
