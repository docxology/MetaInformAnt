#!/usr/bin/env python3
"""Verify the current RNA source modules, producer scripts, and documentation."""

from __future__ import annotations

import ast
import importlib
import inspect
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

CURRENT_MODULES = (
    "metainformant.rna.amalgkit",
    "metainformant.rna.engine.workflow",
    "metainformant.rna.engine.progress_db",
    "metainformant.rna.engine.provenance",
    "metainformant.rna.engine.streaming_orchestrator",
    "metainformant.rna.core.environment",
)
CURRENT_SCRIPTS = (
    REPO_ROOT / "scripts" / "rna" / "run_all_species.py",
    REPO_ROOT / "scripts" / "rna" / "process_species.py",
    REPO_ROOT / "scripts" / "rna" / "check_pipeline_status.py",
    REPO_ROOT / "scripts" / "rna" / "validate_all_species_workflow.py",
)
CURRENT_DOCS = (
    REPO_ROOT / "docs" / "rna" / "README.md",
    REPO_ROOT / "docs" / "rna" / "workflow.md",
    REPO_ROOT / "docs" / "rna" / "CONFIGURATION.md",
    REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "README.md",
)


def verify_modules() -> list[str]:
    """Import current public modules and report undocumented public members."""

    issues: list[str] = []
    for module_name in CURRENT_MODULES:
        try:
            module = importlib.import_module(module_name)
        except Exception as exc:
            issues.append(f"{module_name}: import failed: {exc}")
            continue
        for name, value in inspect.getmembers(module):
            if name.startswith("_") or not inspect.isfunction(value) and not inspect.isclass(value):
                continue
            if getattr(value, "__module__", module_name) != module_name:
                continue
            if not inspect.getdoc(value):
                issues.append(f"{module_name}.{name}: missing docstring")
    return issues


def verify_scripts() -> list[str]:
    """Parse each current production script without executing it."""

    issues: list[str] = []
    for script in CURRENT_SCRIPTS:
        if not script.is_file():
            issues.append(f"missing script: {script}")
            continue
        try:
            ast.parse(script.read_text(encoding="utf-8"), filename=str(script))
        except SyntaxError as exc:
            issues.append(f"{script}: syntax error: {exc}")
    return issues


def verify_docs() -> list[str]:
    """Check that the current documentation anchors and commands exist."""

    issues: list[str] = []
    for document in CURRENT_DOCS:
        if not document.is_file():
            issues.append(f"missing documentation: {document}")
    script_readme = REPO_ROOT / "scripts" / "rna" / "README.md"
    content = script_readme.read_text(encoding="utf-8") if script_readme.is_file() else ""
    for required in ("run_all_species.py", "process_species.py", "check_pipeline_status.py", "AMALGKIT_DATA_ROOT"):
        if required not in content:
            issues.append(f"scripts/rna/README.md: missing {required}")
    return issues


def main() -> int:
    """Run all read-only RNA verification checks."""

    issues = verify_modules() + verify_scripts() + verify_docs()
    if issues:
        print("RNA verification failed:")
        for issue in issues:
            print(f"  - {issue}")
        return 1
    print("RNA verification passed: current modules, scripts, and documentation are complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
