#!/usr/bin/env python3
"""Verify the source-derived core API reference and topic-guide declarations.

``docs/core/API_REFERENCE.md`` is the canonical inventory of public core
symbols.  It is generated from the AST by
``scripts/package/generate_core_api_reference.py``.  This checker compares the
generated inventory with the current source tree and separately scans topic
guide headings for stale API declarations.  Prose and fenced examples are not
treated as declarations, which avoids false positives from ordinary examples.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
SRC_DIR = WORKSPACE / "src" / "metainformant" / "core"
DOCS_DIR = WORKSPACE / "docs" / "core"
DEFAULT_API = DOCS_DIR / "API_REFERENCE.md"
DEFAULT_REPORT = DOCS_DIR / "CROSS_CHECK_REPORT.md"

# Make the source-derived generator importable when this file is executed as a
# script (`uv run python scripts/core_docs_cross_check.py`).
sys.path.insert(0, str(WORKSPACE))
from scripts.package.generate_core_api_reference import (  # noqa: E402
    Symbol,
    collect_symbols,
)


@dataclass(frozen=True)
class ApiEntry:
    """One parsed row from the generated API reference."""

    module: str
    qualified_name: str
    kind: str
    signature: str
    line: int


@dataclass(frozen=True)
class Issue:
    """A check failure written to the human-readable report."""

    check: str
    subject: str
    location: str
    detail: str


MODULE_HEADING = re.compile(r"^##\s+(metainformant\.core\.[^\s{]+)(?:\s+\{#[^}]+\})?\s*$")
API_ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|\s*([^|]+?)\s*\|\s*`((?:\\\||[^`])*)`\s*\|\s*(.*?)\s*\|\s*$")
TOPIC_HEADING = re.compile(r"^#{2,6}\s+`?([A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*)\([^`]*\)`?\s*$")


def _unescape_cell(value: str) -> str:
    """Undo the minimal escaping performed by the API generator."""

    return value.replace(r"\|", "|").strip()


def parse_api_reference(path: Path) -> tuple[list[ApiEntry], list[Issue]]:
    """Parse generated API rows and report malformed or duplicate rows."""

    if not path.exists():
        return [], [Issue("api_reference_present", "API_REFERENCE.md", str(path), "file is missing")]

    entries: list[ApiEntry] = []
    issues: list[Issue] = []
    current_module: str | None = None
    seen: set[tuple[str, str]] = set()

    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        heading = MODULE_HEADING.match(line)
        if heading:
            current_module = heading.group(1)
            continue
        row = API_ROW.match(line)
        if not row or current_module is None:
            continue
        qualified_name = _unescape_cell(row.group(1))
        kind = _unescape_cell(row.group(2))
        signature = _unescape_cell(row.group(3))
        key = (current_module, qualified_name)
        if key in seen:
            issues.append(
                Issue(
                    "api_reference_unique",
                    f"{current_module}.{qualified_name}",
                    f"{path}:{line_number}",
                    "duplicate generated API row",
                )
            )
            continue
        seen.add(key)
        entries.append(ApiEntry(current_module, qualified_name, kind, signature, line_number))

    if not entries:
        issues.append(
            Issue(
                "api_reference_parseable",
                path.name,
                str(path),
                "no module table rows were parsed",
            )
        )
    return entries, issues


def topic_declarations(path: Path) -> list[tuple[str, int]]:
    """Return API-like headings from a hand-written topic guide."""

    declarations: list[tuple[str, int]] = []
    in_fence = False
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        match = TOPIC_HEADING.match(line)
        if match:
            declarations.append((match.group(1), line_number))
    return declarations


def _source_key(symbol: Symbol) -> tuple[str, str]:
    return symbol.module, symbol.qualified_name


def _display_signature(signature: str) -> str:
    """Keep report rows readable without introducing Markdown table breaks."""

    return signature.replace("|", r"\|").replace("`", "'")


def _write_report(
    path: Path,
    *,
    symbols: list[Symbol],
    entries: list[ApiEntry],
    issues: list[Issue],
    topic_count: int,
) -> None:
    """Write a deterministic cross-check report."""

    source_counts = Counter(symbol.module for symbol in symbols)
    api_counts = Counter(entry.module for entry in entries)
    generated = datetime.now(UTC).isoformat(timespec="seconds")
    status = "PASS" if not issues else "REVIEW REQUIRED"

    lines = [
        "# Core Documentation Cross-Check Report",
        "",
        "This report is generated by `scripts/core_docs_cross_check.py`.",
        "The source-derived inventory is generated by",
        "`scripts/package/generate_core_api_reference.py`; topic guides remain",
        "the place for behavioral guidance, examples, and caveats.",
        "",
        f"- Generated (UTC): `{generated}`",
        f"- Status: **{status}**",
        f"- Public source symbols: **{len(symbols)}**",
        f"- API-reference rows: **{len(entries)}**",
        f"- Topic-guide declarations inspected: **{topic_count}**",
        f"- Issues: **{len(issues)}**",
        "",
        "## Module coverage",
        "",
        "| Module | Source symbols | API rows |",
        "| --- | ---: | ---: |",
    ]
    for module in sorted(set(source_counts) | set(api_counts)):
        lines.append(f"| `{module}` | {source_counts[module]} | {api_counts[module]} |")

    lines.extend(["", "## Issues", ""])
    if not issues:
        lines.append("No source/API-reference or stale topic-guide issues were detected.")
    else:
        for index, issue in enumerate(issues, 1):
            lines.extend(
                [
                    f"### {index}. `{issue.check}` — `{issue.subject}`",
                    "",
                    f"- Location: `{issue.location}`",
                    f"- Detail: {issue.detail}",
                    "",
                ]
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_check(api_path: Path) -> tuple[list[Symbol], list[ApiEntry], list[Issue], int]:
    """Run source/API and topic-guide checks."""

    symbols = collect_symbols(SRC_DIR)
    entries, issues = parse_api_reference(api_path)
    source_by_key = {_source_key(symbol): symbol for symbol in symbols}
    api_by_key = {(entry.module, entry.qualified_name): entry for entry in entries}

    if len(source_by_key) != len(symbols):
        issues.append(
            Issue(
                "source_inventory_unique",
                "source symbols",
                str(SRC_DIR),
                "duplicate module-qualified public symbol",
            )
        )

    missing = sorted(set(source_by_key) - set(api_by_key))
    extra = sorted(set(api_by_key) - set(source_by_key))
    for module, qualified_name in missing:
        symbol = source_by_key[(module, qualified_name)]
        issues.append(
            Issue(
                "api_reference_complete",
                f"{module}.{qualified_name}",
                f"{symbol.module}:{qualified_name}",
                "public source symbol is absent from API_REFERENCE.md",
            )
        )
    for module, qualified_name in extra:
        entry = api_by_key[(module, qualified_name)]
        issues.append(
            Issue(
                "api_reference_current",
                f"{module}.{qualified_name}",
                f"{api_path}:{entry.line}",
                "API row has no corresponding public source symbol",
            )
        )

    for key in sorted(set(source_by_key) & set(api_by_key)):
        source = source_by_key[key]
        entry = api_by_key[key]
        if source.signature != entry.signature:
            issues.append(
                Issue(
                    "api_reference_signature",
                    f"{source.module}.{source.qualified_name}",
                    f"{api_path}:{entry.line}",
                    "source="
                    f"`{_display_signature(source.signature)}`; "
                    "api="
                    f"`{_display_signature(entry.signature)}`",
                )
            )

    known_names = {qualified_name for _, qualified_name in source_by_key}
    known_names.update(qualified_name.rsplit(".", 1)[-1] for _, qualified_name in source_by_key)
    topic_count = 0
    excluded = {api_path.name, DEFAULT_REPORT.name, "README.md", "index.md", "SPEC.md"}
    for topic_path in sorted(DOCS_DIR.glob("*.md")):
        if topic_path.name in excluded:
            continue
        for declaration, line_number in topic_declarations(topic_path):
            topic_count += 1
            if declaration not in known_names:
                issues.append(
                    Issue(
                        "topic_declaration_current",
                        declaration,
                        f"{topic_path}:{line_number}",
                        "heading looks like a public API declaration but is not in the source inventory",
                    )
                )

    return symbols, entries, issues, topic_count


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--api", type=Path, default=DEFAULT_API, help="Generated API reference to check")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT, help="Markdown report path")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="return non-zero when any source/API/topic issue is found",
    )
    args = parser.parse_args(argv)
    api_path = args.api.expanduser().resolve()
    report_path = args.report.expanduser().resolve()

    symbols, entries, issues, topic_count = run_check(api_path)
    _write_report(
        report_path,
        symbols=symbols,
        entries=entries,
        issues=issues,
        topic_count=topic_count,
    )

    source_modules = len({symbol.module for symbol in symbols})
    print(f"Source inventory: {len(symbols)} public symbols across {source_modules} modules")
    print(f"API reference: {len(entries)} rows")
    print(f"Topic declarations inspected: {topic_count}")
    print(f"Issues: {len(issues)}")
    print(f"Report: {report_path}")
    if issues:
        for issue in issues[:10]:
            print(f"- {issue.check}: {issue.subject} ({issue.location})")
    return 1 if args.strict and issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
