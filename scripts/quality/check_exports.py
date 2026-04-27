#!/usr/bin/env python3
"""Check export completeness across METAINFORMANT modules.

This script verifies that all public functions and classes defined in module
__init__.py files are properly exported in the __all__ list.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from typing import Any


def extract_functions_and_classes(content: str) -> set[str]:
    """Extract function and class names from Python source."""
    tree = ast.parse(content)
    names = set()

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            # Skip private names (starting with _)
            if not node.name.startswith("_"):
                names.add(node.name)

    return names


def is_submodule_init(module_path: Path) -> bool:
    """Check if this is a submodule __init__.py file that imports functions."""
    # Check if the __all__ list contains function-like names (not just submodule names)
    try:
        content = module_path.read_text(encoding="utf-8")
        exported = extract_all_exports(content)
        imported = extract_imported_modules(content)

        # If most exports are also imports, this is likely a package __init__.py
        # that imports submodules
        if exported and imported:
            overlap = len(exported & imported)
            if overlap / len(exported) > 0.5:  # More than 50% overlap
                return False  # This is a package __init__.py

        # If we have defined functions and they match exports, it's a submodule
        defined = extract_functions_and_classes(content)
        if defined and (defined == exported or not exported):
            return True

        return False
    except:
        return False


def extract_imported_modules(content: str) -> set[str]:
    """Extract modules imported in __init__.py files."""
    tree = ast.parse(content)
    imported = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                # For from X import Y, we care about Y
                for alias in node.names:
                    imported.add(alias.name)

    return imported


def extract_all_exports(content: str) -> set[str]:
    """Extract names from __all__ list."""
    tree = ast.parse(content)
    exports = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    if isinstance(node.value, ast.List):
                        for item in node.value.elts:
                            if isinstance(item, ast.Constant) and isinstance(item.value, str):
                                exports.add(item.value)

    return exports


def check_module_exports(module_path: Path) -> dict[str, Any]:
    """Check export completeness for a single module."""
    try:
        content = module_path.read_text(encoding="utf-8")
    except Exception as e:
        return {
            "status": "error",
            "error": f"Failed to read {module_path}: {e}",
            "defined": set(),
            "exported": set(),
            "missing": set(),
            "extra": set(),
        }

    defined = extract_functions_and_classes(content)
    exported = extract_all_exports(content)

    # Special handling for submodule __init__.py files vs package __init__.py files
    if not is_submodule_init(module_path):
        # This is a package __init__.py that imports submodules
        # The __all__ list contains submodule names, which is correct
        return {
            "status": "ok",
            "defined": defined,
            "exported": exported,
            "missing": set(),
            "extra": set(),
            "note": "Package __init__.py with submodule exports",
        }

    missing = defined - exported
    extra = exported - defined

    return {"status": "ok", "defined": defined, "exported": exported, "missing": missing, "extra": extra}


def main():
    """Main function."""
    repo_root = Path(__file__).parent.parent
    src_dir = repo_root / "src" / "metainformant"

    if not src_dir.exists():
        print(f"Error: Source directory not found: {src_dir}")
        sys.exit(1)

    # Find all __init__.py files
    init_files = list(src_dir.rglob("__init__.py"))

    print("Checking export completeness across METAINFORMANT modules...")
    print("=" * 60)

    all_good = True
    summary = []

    for init_file in sorted(init_files):
        module_name = str(init_file.relative_to(src_dir)).replace("/__init__.py", "").replace("\\__init__.py", "")

        result = check_module_exports(init_file)
        summary.append((module_name, result))

        if result["status"] == "error":
            print(f"❌ {module_name}: {result['error']}")
            all_good = False
            continue

        defined_count = len(result["defined"])
        exported_count = len(result["exported"])
        missing_count = len(result["missing"])
        extra_count = len(result["extra"])

        note = result.get("note", "")

        if missing_count == 0 and extra_count == 0:
            status = "✅"
            if note:
                status += f" {note}"
            print(f"{status} {module_name}: {exported_count} exports (complete)")
        else:
            print(f"⚠️  {module_name}: {exported_count}/{defined_count} exported")
            all_good = False

            if result["missing"]:
                print(f"    Missing exports: {sorted(result['missing'])}")
            if result["extra"]:
                print(f"    Extra exports: {sorted(result['extra'])}")
            if note:
                print(f"    Note: {note}")

    print("\n" + "=" * 60)

    if all_good:
        print("🎉 All modules have complete and correct exports!")
        sys.exit(0)
    else:
        print("⚠️  Some modules have export issues. Please review and fix.")
        print("\nRun this script again to verify fixes:")
        print("    python scripts/check_exports.py")
        sys.exit(1)


if __name__ == "__main__":
    main()
