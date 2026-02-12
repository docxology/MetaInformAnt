#!/usr/bin/env python3
"""Build import mapping from __init__.py re-exports to canonical submodule paths.

This script:
1. Parses each module's __init__.py to find function/class re-exports
2. Traces each export to its canonical submodule path
3. Scans test files and source files for shortcut imports
4. Generates canonical replacement imports

Usage:
    python scripts/reorganize/build_import_map.py [--module MODULE] [--fix] [--dry-run]
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent.parent / "src" / "metainformant"
TESTS = Path(__file__).resolve().parent.parent.parent / "tests"
ROOT = Path(__file__).resolve().parent.parent.parent

MODULES = [
    "menu", "quality", "singlecell", "phenotype", "multiomics",
    "ontology", "protein", "ml", "structural_variants", "metagenomics",
    "networks", "dna", "epigenome", "pharmacogenomics", "longread", "spatial",
    "ecology", "simulation", "life_events", "math", "information",
    "rna", "gwas", "visualization", "core",
]


def find_canonical_source(module_name: str, symbol: str) -> str | None:
    """Find canonical source path for a symbol exported from a module's __init__.py.

    Searches all .py files in the module for the symbol definition.
    Returns the canonical import path like 'metainformant.ecology.analysis.community'.
    """
    mod_dir = SRC / module_name

    for py_file in sorted(mod_dir.rglob("*.py")):
        if py_file.name == "__init__.py":
            continue
        if "__pycache__" in str(py_file):
            continue

        try:
            source = py_file.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(py_file))
        except (SyntaxError, UnicodeDecodeError):
            continue

        for node in ast.walk(tree):
            # Check function/class definitions
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                if node.name == symbol:
                    # Build canonical path
                    rel = py_file.relative_to(SRC)
                    parts = list(rel.with_suffix("").parts)
                    return "metainformant." + ".".join(parts)

            # Check top-level assignments (constants, dataclasses, etc.)
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == symbol:
                        rel = py_file.relative_to(SRC)
                        parts = list(rel.with_suffix("").parts)
                        return "metainformant." + ".".join(parts)

    return None


def parse_init_exports(module_name: str) -> dict[str, str]:
    """Parse a module's __init__.py to extract all re-exported symbols.

    Returns {symbol_name: 'from_path'} for imports like:
        from .analysis.community import shannon_diversity
    -> {'shannon_diversity': '.analysis.community'}
    """
    init_file = SRC / module_name / "__init__.py"
    if not init_file.exists():
        return {}

    try:
        source = init_file.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(init_file))
    except (SyntaxError, UnicodeDecodeError):
        return {}

    exports: dict[str, str] = {}

    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.ImportFrom) and node.level >= 1:
            mod_path = node.module or ""
            for alias in node.names:
                name = alias.asname if alias.asname else alias.name
                exports[name] = mod_path

        # Also check inside try/except blocks
        if isinstance(node, ast.Try):
            for child in node.body:
                if isinstance(child, ast.ImportFrom) and child.level >= 1:
                    mod_path = child.module or ""
                    for alias in child.names:
                        name = alias.asname if alias.asname else alias.name
                        exports[name] = mod_path

    return exports


def find_shortcut_imports_in_file(filepath: Path, module_name: str) -> list[dict]:
    """Find imports in a file that use shortcut paths through __init__.py.

    Returns list of {line_num, original_line, symbol, from_module}.
    """
    results = []
    try:
        source = filepath.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(filepath))
    except (SyntaxError, UnicodeDecodeError):
        return results

    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            # Match: from metainformant.<module> import <symbol>
            if node.module == f"metainformant.{module_name}":
                for alias in node.names:
                    results.append({
                        "line_num": node.lineno,
                        "module": node.module,
                        "symbol": alias.name,
                        "asname": alias.asname,
                    })

    return results


def build_full_map() -> dict[str, dict[str, str]]:
    """Build complete import mapping for all modules.

    Returns {module_name: {symbol: canonical_path}}.
    """
    full_map = {}

    for mod in MODULES:
        exports = parse_init_exports(mod)
        symbol_map = {}

        # Filter out submodule imports (we only care about function/class re-exports)
        for symbol, from_path in exports.items():
            # Skip if it's just a submodule import (e.g., from . import analysis)
            # These are fine to keep
            if from_path == ".":
                continue

            # Find canonical source
            canonical = find_canonical_source(mod, symbol)
            if canonical:
                symbol_map[symbol] = canonical

        if symbol_map:
            full_map[mod] = symbol_map

    return full_map


def scan_files_for_shortcuts(full_map: dict[str, dict[str, str]]) -> list[dict]:
    """Scan all test and source files for shortcut imports that need updating."""
    issues = []

    # Scan test files
    for test_file in sorted(TESTS.glob("*.py")):
        for mod, symbol_map in full_map.items():
            shortcuts = find_shortcut_imports_in_file(test_file, mod)
            for s in shortcuts:
                if s["symbol"] in symbol_map:
                    issues.append({
                        "file": str(test_file),
                        "line": s["line_num"],
                        "module": mod,
                        "symbol": s["symbol"],
                        "asname": s["asname"],
                        "canonical_path": symbol_map[s["symbol"]],
                    })

    # Scan source files (cross-module imports)
    for py_file in sorted(SRC.rglob("*.py")):
        if "__pycache__" in str(py_file):
            continue
        if py_file.name == "__init__.py":
            continue

        for mod, symbol_map in full_map.items():
            # Don't flag intra-module imports
            if str(py_file).startswith(str(SRC / mod)):
                continue

            shortcuts = find_shortcut_imports_in_file(py_file, mod)
            for s in shortcuts:
                if s["symbol"] in symbol_map:
                    issues.append({
                        "file": str(py_file),
                        "line": s["line_num"],
                        "module": mod,
                        "symbol": s["symbol"],
                        "asname": s["asname"],
                        "canonical_path": symbol_map[s["symbol"]],
                    })

    return issues


def generate_fix(issue: dict) -> tuple[str, str]:
    """Generate old and new import line for a fix."""
    symbol = issue["symbol"]
    canonical = issue["canonical_path"]
    asname = issue["asname"]

    old_import = f"from metainformant.{issue['module']} import {symbol}"
    if asname:
        new_import = f"from {canonical} import {symbol} as {asname}"
    else:
        new_import = f"from {canonical} import {symbol}"

    return old_import, new_import


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build import mapping and find shortcut imports")
    parser.add_argument("--module", "-m", help="Only process specific module")
    parser.add_argument("--fix", action="store_true", help="Apply fixes to files")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be changed")
    parser.add_argument("--map-only", action="store_true", help="Only show the import map")
    args = parser.parse_args()

    print("Building import map...")
    full_map = build_full_map()

    if args.module:
        if args.module in full_map:
            full_map = {args.module: full_map[args.module]}
        else:
            print(f"Module {args.module} has no re-exports to map")
            return

    # Print the map
    total_symbols = 0
    for mod, symbols in sorted(full_map.items()):
        total_symbols += len(symbols)
        if args.map_only or not args.fix:
            print(f"\n=== {mod} ({len(symbols)} re-exports) ===")
            for sym, path in sorted(symbols.items()):
                print(f"  {sym} -> {path}")

    print(f"\nTotal: {total_symbols} re-exported symbols across {len(full_map)} modules")

    if args.map_only:
        return

    # Scan for shortcut imports
    print("\nScanning for shortcut imports...")
    issues = scan_files_for_shortcuts(full_map)

    # Group by file
    by_file: dict[str, list[dict]] = {}
    for issue in issues:
        by_file.setdefault(issue["file"], []).append(issue)

    print(f"Found {len(issues)} shortcut imports in {len(by_file)} files")

    for filepath, file_issues in sorted(by_file.items()):
        rel = Path(filepath).relative_to(ROOT)
        print(f"\n  {rel}:")
        for issue in file_issues:
            old, new = generate_fix(issue)
            print(f"    L{issue['line']}: {old}")
            print(f"         -> {new}")

    if args.fix and not args.dry_run:
        print("\n\nApplying fixes...")
        fix_files(by_file)
        print("Done!")
    elif args.dry_run:
        print("\n(dry run - no changes made)")


def fix_files(by_file: dict[str, list[dict]]):
    """Apply import fixes to files."""
    for filepath, issues in by_file.items():
        path = Path(filepath)
        content = path.read_text(encoding="utf-8")
        lines = content.split("\n")

        # Group imports by line to handle multi-import lines
        # Sort by line number in reverse to avoid offset issues
        issues_sorted = sorted(issues, key=lambda x: x["line"], reverse=True)

        for issue in issues_sorted:
            line_idx = issue["line"] - 1
            if line_idx < len(lines):
                old_line = lines[line_idx]
                old, new = generate_fix(issue)
                # Simple replacement - may need refinement for multi-symbol imports
                if old in old_line:
                    lines[line_idx] = old_line.replace(old, new)

        path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
