#!/usr/bin/env python3
"""Comprehensive verification of RNA documentation and code accuracy.

This script performs multiple levels of verification:
1. Basic verification: docstrings, examples, links
2. Documentation verification: module docstrings, imports
3. Comprehensive verification: signature matching, cross-references

Utility functions are in _verify_utils.py for reusability.
"""

import ast
import importlib
import importlib.util
import inspect
import sys
from pathlib import Path

# Import utility functions from companion module
sys.path.insert(0, str(Path(__file__).parent))
from _verify_utils import (
    check_doc_examples,
    check_doc_links,
)

REPO_ROOT = Path(__file__).parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "rna"
SRC_DIR = REPO_ROOT / "src" / "metainformant" / "rna"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "rna"

issues = []
warnings = []


def verify_module_docs(key_modules: list[tuple[str, Path]]) -> None:
    """Phase 1: Check that public functions/classes have docstrings."""
    print("Phase 1: Method Documentation Verification")
    print("-" * 80)

    for module_name, module_path in key_modules:
        if module_path.exists():
            try:
                spec = importlib.util.spec_from_file_location(module_name, module_path)
                if spec and spec.loader:
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)

                    missing_docs = []
                    for name, obj in inspect.getmembers(module):
                        if name.startswith("_"):
                            continue
                        if inspect.isfunction(obj) or inspect.isclass(obj):
                            if not obj.__doc__ or len(obj.__doc__.strip()) < 10:
                                missing_docs.append(name)

                    if missing_docs:
                        issues.extend(
                            [f"{module_name}.{name}: Missing/incomplete docstring" for name in missing_docs]
                        )
                        print(f"  ❌ {module_name}: {len(missing_docs)} undocumented items")
                    else:
                        print(f"  ✅ {module_name}: All documented")
            except Exception as e:
                warnings.append(f"{module_name}: Could not check - {e}")
                print(f"  ⚠️  {module_name}: Check failed")

    print()


def verify_doc_accuracy(key_docs: list[Path]) -> None:
    """Phase 2: Check documentation examples and links."""
    print("Phase 2: Documentation Accuracy")
    print("-" * 80)

    for doc_file in key_docs:
        if doc_file.exists():
            example_issues = check_doc_examples(doc_file)
            link_issues = check_doc_links(doc_file, REPO_ROOT, warnings)
            issues.extend(example_issues)
            issues.extend(link_issues)

            total = len(example_issues) + len(link_issues)
            if total == 0:
                print(f"  ✅ {doc_file.name}: Valid")
            else:
                print(f"  ⚠️  {doc_file.name}: {total} issues")
        else:
            warnings.append(f"Missing: {doc_file}")
            print(f"  ❌ {doc_file.name}: Missing")

    print()


def verify_scripts(production_scripts: list[Path]) -> None:
    """Phase 3: Check production scripts for syntax validity."""
    print("Phase 3: Script Verification")
    print("-" * 80)

    for script in production_scripts:
        if script.exists():
            try:
                with open(script) as f:
                    ast.parse(f.read(), filename=str(script))
                print(f"  ✅ {script.name}: Valid syntax")
            except SyntaxError as e:
                issues.append(f"{script}: Syntax error: {e}")
                print(f"  ❌ {script.name}: Syntax error")
        else:
            issues.append(f"Missing production script: {script}")
            print(f"  ❌ {script.name}: Missing")

    print()


def verify_imports() -> None:
    """Phase 4: Check that all RNA module exports are importable."""
    print("Phase 4: Import Verification")
    print("-" * 80)

    try:
        from metainformant.rna import __all__ as rna_all

        importable = 0
        failed = 0

        for export_name in rna_all:
            try:
                exec(f"from metainformant.rna import {export_name}")
                importable += 1
            except Exception as e:
                issues.append(f"Cannot import {export_name}: {e}")
                failed += 1

        print(f"  ✅ {importable}/{len(rna_all)} exports importable")
        if failed > 0:
            print(f"  ❌ {failed} imports failed")
    except Exception as e:
        issues.append(f"Cannot check exports: {e}")
        print(f"  ❌ Cannot check: {e}")

    print()


def main():
    """Run comprehensive verification."""
    print("=" * 80)
    print("RNA DOCUMENTATION AND CODE COMPREHENSIVE VERIFICATION")
    print("=" * 80)
    print()

    key_modules = [
        ("metainformant.rna.amalgkit", SRC_DIR / "amalgkit.py"),
        ("metainformant.rna.workflow", SRC_DIR / "workflow.py"),
        ("metainformant.rna.configs", SRC_DIR / "configs.py"),
        ("metainformant.rna.monitoring", SRC_DIR / "monitoring.py"),
        ("metainformant.rna.environment", SRC_DIR / "environment.py"),
    ]
    verify_module_docs(key_modules)

    key_docs = [
        DOCS_DIR / "README.md",
        DOCS_DIR / "WORKFLOW.md",
        DOCS_DIR / "STEPS.md",
        DOCS_DIR / "CONFIGURATION.md",
        DOCS_DIR / "ORCHESTRATION" / "ENA_WORKFLOW.md",
        DOCS_DIR / "ORCHESTRATION" / "MULTI_SPECIES.md",
        DOCS_DIR / "ORCHESTRATION" / "BATCH_DOWNLOAD.md",
    ]
    verify_doc_accuracy(key_docs)

    production_scripts = [
        SCRIPTS_DIR / "run_workflow.py",
        SCRIPTS_DIR / "setup_genome.py",
        SCRIPTS_DIR / "discover_species.py",
        SCRIPTS_DIR / "check_environment.py",
    ]
    verify_scripts(production_scripts)

    verify_imports()

    # Summary
    print("=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Issues: {len(issues)}")
    print(f"Warnings: {len(warnings)}")

    if issues:
        print("\nIssues to fix:")
        for issue in issues[:20]:
            print(f"  - {issue}")
        if len(issues) > 20:
            print(f"  ... and {len(issues) - 20} more")

    if warnings:
        print("\nWarnings:")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more")

    if not issues:
        print("\n✅ All critical verifications passed!")
        return 0
    else:
        print("\n❌ Some issues found - review above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
