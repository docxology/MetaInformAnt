#!/usr/bin/env python3
"""
Project completeness audit — check documentation and test infrastructure across all modules.

Usage:
  python scripts/validate/project_completeness_audit.py [--root PATH]

Outputs:
  Completeness scorecard (Markdown)

This is the canonical implementation of the "project-completeness-audit" pattern.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def discover_modules(src_root: Path) -> list[str]:
    return sorted(
        d.name
        for d in src_root.iterdir()
        if d.is_dir() and not d.name.startswith("_")
    )


def check_doc_files(module_dir: Path) -> dict[str, bool]:
    return {
        "AGENTS.md": (module_dir / "AGENTS.md").exists(),
        "README.md": (module_dir / "README.md").exists(),
        "SPEC.md": (module_dir / "SPEC.md").exists(),
    }


def sample_docstring_coverage(py_files: list[Path], sample_size: int = 100) -> float:
    if not py_files:
        return 0.0
    sample = py_files[:sample_size] if len(py_files) > sample_size else py_files
    with_docs = 0
    for f in sample:
        try:
            text = f.read_text(errors="ignore")
            if '"""' in text or "'''" in text:
                with_docs += 1
        except Exception:
            pass
    return (with_docs / len(sample)) * 100


def find_test_files(tests_root: Path, module: str) -> list[Path]:
    test_dir = tests_root / module
    if not test_dir.exists():
        return []
    return list(test_dir.glob("test_*.py"))


def check_placeholder_content(md_file: Path) -> list[str]:
    try:
        text = md_file.read_text(errors="ignore").lower()
    except Exception:
        return []
    markers = []
    for keyword in ["todo", "fixme", "xxx", "to be written", "coming soon", "tbd", "placeholder"]:
        if keyword in text:
            markers.append(keyword)
    return markers


def run_validation_scripts(repo_root: Path) -> dict[str, str]:
    results = {}
    scripts = [
        "scripts/maintenance/check_docs_spec.py",
        "scripts/quality/check_exports.py",
    ]
    for script_rel in scripts:
        script = repo_root / script_rel
        if script.exists():
            try:
                result = subprocess.run(
                    [sys.executable, str(script)],
                    capture_output=True,
                    text=True,
                    timeout=30,
                    cwd=repo_root,
                )
                results[script_rel] = (
                    f"exit={result.returncode}\\n"
                    f"stdout:\\n{result.stdout[:500]}\\n"
                    f"stderr:\\n{result.stderr[:200]}"
                )
            except Exception as e:
                results[script_rel] = f"error: {e}"
    return results


def can_import_pytest() -> bool:
    try:
        import pytest  # noqa: F401
        return True
    except ImportError:
        return False


def generate_report(
    repo_root: Path,
    modules: list[str],
    doc_matrix: dict[str, dict[str, bool]],
    test_counts: dict[str, int],
    docstring_coverage: float,
    placeholder_issues: dict[str, list[str]],
    validation_outputs: dict[str, str],
    pytest_available: bool,
) -> str:
    lines = []
    lines.append("# Project Completeness Audit")
    lines.append("")
    lines.append(f"**Repository**: `{repo_root}`")
    lines.append(f"**Modules**: {len(modules)}")
    lines.append(f"**Docstring coverage**: {docstring_coverage:.1f}%")
    lines.append(f"**pytest available**: {'yes' if pytest_available else 'no'}")
    lines.append("")

    required_checks = [
        all(doc_matrix[m]["AGENTS.md"] for m in modules),
        all(doc_matrix[m]["README.md"] for m in modules),
        docstring_coverage >= 95.0,
        all(test_counts[m] > 0 for m in modules),
        pytest_available,
        not any(placeholder_issues.values()),
    ]
    score_pct = int(sum(required_checks) / len(required_checks) * 100)
    lines.append(f"## Summary\\n\\n**Overall completeness: {score_pct}%**")
    if score_pct == 100:
        lines.append("✅ All required checks passed.")
    elif score_pct >= 80:
        lines.append("⚠️  Mostly complete with minor gaps.")
    else:
        lines.append("❌ Significant completeness gaps detected.")
    lines.append("")

    lines.append("## Documentation Coverage\\n")
    lines.append("| Module | AGENTS.md | README.md | SPEC.md | Placeholders |")
    lines.append("|--------|-----------|-----------|---------|--------------|")
    for m in modules:
        agents = "✓" if doc_matrix[m]["AGENTS.md"] else "✗"
        readme = "✓" if doc_matrix[m]["README.md"] else "✗"
        spec = "✓" if doc_matrix[m]["SPEC.md"] else "✗"
        ph = ", ".join(placeholder_issues.get(m, [])) or "—"
        lines.append(f"| {m} | {agents} | {readme} | {spec} | {ph} |")
    lines.append("")

    lines.append("## Test Coverage\\n")
    lines.append("| Module | Test files |")
    lines.append("|--------|------------|")
    for m in modules:
        count = test_counts[m]
        lines.append(f"| {m} | {count} |")
    lines.append("")

    gaps = []
    if not all(doc_matrix[m]["AGENTS.md"] for m in modules):
        gaps.append("AGENTS.md missing from one or more modules")
    if not all(doc_matrix[m]["README.md"] for m in modules):
        gaps.append("README.md missing from one or more modules")
    if not all(doc_matrix[m]["SPEC.md"] for m in modules):
        gaps.append("SPEC.md missing from some modules (optional but recommended)")
    if docstring_coverage < 95:
        gaps.append(f"Docstring coverage {docstring_coverage:.1f}% below 95% threshold")
    if any(test_counts[m] == 0 for m in modules):
        no_tests = [m for m in modules if test_counts[m] == 0]
        gaps.append(f"No test files for: {', '.join(no_tests)}")
    if not pytest_available:
        gaps.append("pytest not importable — cannot execute test suite")
    if any(placeholder_issues.values()):
        gaps.append("Placeholder markers found in some AGENTS.md files")

    if gaps:
        lines.append("## Gaps\\n")
        for g in gaps:
            lines.append(f"- ❌ {g}")
        lines.append("")

    if validation_outputs:
        lines.append("## Validation Scripts\\n")
        for script, output in validation_outputs.items():
            lines.append(f"### {script}\\n")
            lines.append("```")
            lines.append(output[:1000])
            lines.append("```")
            lines.append("")

    lines.append("## Recommended Actions\\n")
    if not pytest_available:
        lines.append("1. Install pytest: `uv pip install pytest` or `pip install pytest`")
    if not all(doc_matrix[m]["SPEC.md"] for m in modules):
        missing_spec = [m for m in modules if not doc_matrix[m]["SPEC.md"]]
        lines.append(f"2. Create SPEC.md for: {', '.join(missing_spec)}")
    if any(test_counts[m] == 0 for m in modules):
        no_tests = [m for m in modules if test_counts[m] == 0]
        lines.append(f"3. Add tests for: {', '.join(no_tests)}")
    if docstring_coverage < 95:
        lines.append("4. Improve docstring coverage (target: ≥95%)")
    if any(placeholder_issues.values()):
        lines.append("5. Replace placeholder text in AGENTS.md files with actual content")

    lines.append("")
    lines.append("---")
    lines.append("*Audit completed*")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Project completeness audit")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository root")
    args = parser.parse_args()

    repo_root = args.root.resolve()
    src_root = repo_root / "src" / "metainformant"
    tests_root = repo_root / "tests"
    if not src_root.exists():
        print(f"Error: {src_root} does not exist", file=sys.stderr)
        return 1

    modules = discover_modules(src_root)
    print(f"Discovered {len(modules)} modules: {', '.join(modules)}")

    doc_matrix = {}
    docstring_coverage = 0.0
    placeholder_issues = {}

    all_py_files = list(src_root.rglob("*.py"))
    docstring_coverage = sample_docstring_coverage(all_py_files)

    for module in modules:
        module_dir = src_root / module
        doc_matrix[module] = check_doc_files(module_dir)
        agents_file = module_dir / "AGENTS.md"
        if agents_file.exists():
            ph = check_placeholder_content(agents_file)
            if ph:
                placeholder_issues[module] = ph

    test_counts = {}
    for module in modules:
        test_counts[module] = len(find_test_files(tests_root, module))

    validation_outputs = run_validation_scripts(repo_root)
    pytest_available = can_import_pytest()

    report = generate_report(
        repo_root,
        modules,
        doc_matrix,
        test_counts,
        docstring_coverage,
        placeholder_issues,
        validation_outputs,
        pytest_available,
    )
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
