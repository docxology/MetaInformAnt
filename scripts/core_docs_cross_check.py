#!/usr/bin/env python3
"""
Cross-check consistency between core documentation and source code.
Performs AST-based extraction of functions from Python files and
parses markdown documentation for API signatures, then compares.
"""

import ast
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

# Base paths
WORKSPACE = Path("/home/trim/Documents/Git/MetaInformAnt")
SRC_DIR = WORKSPACE / "src" / "metainformant" / "core"
DOCS_DIR = WORKSPACE / "docs" / "core"

# Core documentation modules to check
DOC_MODULES = ["cache", "config", "db", "download", "hash", "io", "logging", "parallel", "paths", "text", "workflow"]

# Map each documented module to its source file paths (relative to src/metainformant/core/)
MODULE_SOURCE_FILES = {
    "cache": [
        "io/cache.py",
    ],
    "config": [
        "utils/config.py",
    ],
    "db": [
        "data/db.py",
        "data/validation.py",  # likely part of db module
    ],
    "download": [
        "io/download.py",
        "io/download_manager.py",
        "io/download_robust.py",
    ],
    "hash": [
        "utils/hash.py",
    ],
    "io": [
        "io/io.py",
        "io/atomic.py",
        "io/disk.py",
        "io/checksums.py",
        "io/errors.py",
        # Note: paths, cache, download are separate modules
    ],
    "logging": [
        "utils/logging.py",
    ],
    "parallel": [
        "execution/parallel.py",
    ],
    "paths": [
        "io/paths.py",
    ],
    "text": [
        "utils/text.py",
    ],
    "workflow": [
        "execution/workflow.py",
        "execution/discovery.py",
        "engine/workflow_manager.py",
    ],
}

# Extra files not in main modules (like ui, utils non-documented) - skip from main report
OTHER_SOURCE_FILES = [
    "ui/tui.py",
    "utils/errors.py",
    "utils/optional_deps.py",
    "utils/progress.py",
    "utils/symbols.py",
    "utils/timing.py",
    "utils/watchdog.py",
    "data/__init__.py",
    "engine/__init__.py",
    "execution/__init__.py",
    "utils/__init__.py",
    "io/__init__.py",
]


def extract_python_symbols(filepath: Path) -> Tuple[Dict[str, Dict[str, Any]], Set[str]]:
    """
    Extract all function definitions and class names from a Python file using AST.
    Returns (functions_dict, class_names_set).
    """
    try:
        with open(filepath, "r") as f:
            content = f.read()
        tree = ast.parse(content, filename=str(filepath))
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return {}, set()

    functions = {}
    class_names = set()

    class SymbolVisitor(ast.NodeVisitor):
        def __init__(self):
            self.parent_stack = []

        def generic_visit(self, node):
            for child in ast.iter_child_nodes(node):
                self.parent_stack.append(node)
                self.visit(child)
                self.parent_stack.pop()

        def visit_FunctionDef(self, node):
            # Skip private/internal functions starting with _
            if node.name.startswith("_") and not node.name.startswith("__"):
                return

            # Check if immediate parent is a ClassDef (i.e., it's a method)
            is_method = any(isinstance(p, ast.ClassDef) for p in self.parent_stack)

            # Extract parameters
            args = node.args
            params = []
            defaults_start = len(args.args) - len(args.defaults)

            for i, arg in enumerate(args.args):
                param_name = arg.arg
                # Skip 'self' and 'cls' for methods
                if is_method and param_name in ("self", "cls"):
                    continue
                default = None
                if i >= defaults_start:
                    default_idx = i - defaults_start
                    default_node = args.defaults[default_idx]
                    default = ast.unparse(default_node) if hasattr(ast, "unparse") else str(default_node)
                params.append((param_name, default))

            # Handle *args and **kwargs
            if args.vararg:
                params.append((f"*{args.vararg.arg}", None))
            if args.kwarg:
                params.append((f"**{args.kwarg.arg}", None))

            # Build signature string
            sig_parts = []
            for pname, pdefault in params:
                if pdefault is not None:
                    sig_parts.append(f"{pname}={pdefault}")
                else:
                    sig_parts.append(pname)
            signature = f"({', '.join(sig_parts)})"

            functions[node.name] = {
                "signature": signature,
                "line_no": node.lineno,
                "is_async": isinstance(node, ast.AsyncFunctionDef),
                "is_method": is_method,
                "file": filepath.name,
                "rel_path": str(filepath.relative_to(SRC_DIR.parent.parent)),
            }

        def visit_ClassDef(self, node):
            class_names.add(node.name)
            # Continue visiting class body to capture methods
            self.generic_visit(node)

    visitor = SymbolVisitor()
    visitor.visit(tree)
    return functions, class_names


def extract_documented_functions(filepath: Path) -> Dict[str, Dict[str, Any]]:
    """
    Extract documented function signatures from a markdown file.
    Searches the entire file for function definitions in code blocks, inline code, and headers/tables.
    """
    try:
        with open(filepath, "r") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return {}

    functions = {}

    # Pattern 1: ```python ... def function(...):```
    # Only capture from code blocks that are signature-only (single-line or just def line)
    code_blocks = re.findall(r"```python\s*\n(.*?)\n```", content, re.DOTALL | re.IGNORECASE)
    for block in code_blocks:
        # Skip multi-line code blocks (likely full examples)
        lines = [ln.strip() for ln in block.splitlines() if ln.strip()]
        if len(lines) > 1:
            continue
        # Single-line block, look for function definition
        line = lines[0] if lines else ""
        m = re.search(r"def\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*(\(.*\))", line)
        if m:
            fname = m.group(1)
            if fname.startswith("_"):
                continue
            params_str = m.group(2).strip()
            signature = params_str
            functions[fname] = {
                "signature": signature,
                "line_no": content[: content.find(m.group(0))].count("\n") + 1,
                "source": "code_block",
            }
    # Pattern 2: Inline code with full signature (type-annotated): `function_name(param: Type = default, ...) -> Return`
    # Only capture if parameters contain type hints (':') or return arrow ('->') to avoid capturing function calls
    inline_pattern = r"`([a-zA-Z_][a-zA-Z0-9_]*)\s*\(([^`]*?)\)`"
    for m in re.finditer(inline_pattern, content):
        fname = m.group(1)
        params_str = m.group(2).strip()
        if fname in functions:
            continue
        if fname.startswith("_"):
            continue
        # Require type hint presence to distinguish signatures from calls
        if not (":" in params_str or "->" in params_str):
            continue
        signature = f"({params_str})" if params_str else "()"
        functions[fname] = {
            "signature": signature,
            "line_no": content[: content.find(m.group(0))].count("\n") + 1,
            "source": "inline_code",
        }

    # Pattern 3: Function headers like `#### function_name(params)`
    header_pattern = r"^(#{1,4})\s+(.+?)\s*$"
    for m in re.finditer(header_pattern, content, re.MULTILINE):
        line = m.group(2).strip()
        if line.startswith("`") and line.endswith("`"):
            line = line[1:-1].strip()
        # Skip markdown table rows (lines starting with '|')
        if line.startswith("|"):
            continue
        paren_match = re.search(r"^([a-zA-Z_][a-zA-Z0-9_]*)\s*(\(.*\))", line)
        if not paren_match:
            continue
        fname = paren_match.group(1)
        signature = paren_match.group(2)
        if not (signature.startswith("(") and signature.endswith(")")):
            continue
        if fname not in functions and not fname.startswith("_"):
            functions[fname] = {
                "signature": signature,
                "line_no": content[: content.find(m.group(0))].count("\n") + 1,
                "source": "header",
            }

    # Pattern 4: Markdown table rows
    table_pattern = r"^\|\s*`?([a-zA-Z_][a-zA-Z0-9_]*?)`?\s*\|\s*`?([^`|]+?)\s*(?:`\s*\|?|$)\s*\|"
    for m in re.finditer(table_pattern, content, re.MULTILINE):
        fname = m.group(1).strip()
        sig_raw = m.group(2).strip()
        has_signature_markers = "(" in sig_raw or "->" in sig_raw
        if not has_signature_markers:
            continue
        sig_match = re.search(r"\([^)]*\)", sig_raw)
        if sig_match:
            signature = sig_match.group(0)
        else:
            signature = "()"
        if fname not in functions and not fname.startswith("_"):
            functions[fname] = {
                "signature": signature,
                "line_no": content[: content.find(m.group(0))].count("\n") + 1,
                "source": "table",
            }

    return functions


def normalize_signature(sig: str) -> str:
    """
    Normalize signature for comparison: strip spaces around = and commas,
    remove type hints, convert defaults to canonical form.
    """
    # Remove type hints in params: `x: int` -> `x`
    sig = re.sub(r":\s*[^,)=]+", "", sig)
    # Normalize spaces around = and commas
    sig = re.sub(r"\s*=\s*", "=", sig)
    sig = re.sub(r"\s*,\s*", ", ", sig)
    # Strip outer parentheses
    sig = sig.strip()
    if sig.startswith("(") and sig.endswith(")"):
        sig = sig[1:-1]
    return sig.strip()


def compare_signatures(code_sig: str, doc_sig: str) -> Tuple[bool, List[str]]:
    """
    Compare two signatures, allowing for some flexibility.
    Returns (match, differences)
    """
    norm_code = normalize_signature(code_sig)
    norm_doc = normalize_signature(doc_sig)

    if norm_code == norm_doc:
        return True, []

    differences = []
    # Parse params
    code_params = [p.strip() for p in norm_code.split(",") if p.strip()]
    doc_params = [p.strip() for p in norm_doc.split(",") if p.strip()]

    code_names = set()
    for p in code_params:
        if "=" in p:
            name = p.split("=")[0].strip()
            code_names.add(name)
        else:
            code_names.add(p)

    doc_names = set()
    for p in doc_params:
        if "=" in p:
            name = p.split("=")[0].strip()
            doc_names.add(name)
        else:
            doc_names.add(p)

    missing_in_doc = code_names - doc_names
    extra_in_doc = doc_names - code_names

    if missing_in_doc:
        differences.append(f"Missing in docs: {sorted(missing_in_doc)}")
    if extra_in_doc:
        differences.append(f"Extra in docs: {sorted(extra_in_doc)}")

    return False, differences


def main():
    print("=" * 80)
    print("CORE DOCUMENTATION CROSS-CHECK")
    print("=" * 80)

    # Step 1: Build source code inventory grouped by documented module name
    print("\n[1/3] Parsing source code...")
    code_inventory = {mod: {} for mod in DOC_MODULES}  # module -> {func_name: info}
    all_class_names = set()  # Set of all class names to exclude from doc comparison
    all_code_files = []

    # Process each documented module's source files
    for mod_name, src_files in MODULE_SOURCE_FILES.items():
        for src_file in src_files:
            src_path = SRC_DIR / src_file
            if src_path.exists():
                funcs, classes = extract_python_symbols(src_path)
                code_inventory[mod_name].update(funcs)
                all_class_names.update(classes)
                all_code_files.append(src_path)
            else:
                print(f"  Warning: {src_path} not found")

    # Also process other source files (not part of documented modules)
    for extra in OTHER_SOURCE_FILES:
        src_path = SRC_DIR / extra
        if src_path.exists():
            # These are not part of the core module comparison
            pass

    total_code_funcs = sum(len(v) for v in code_inventory.values())
    print(f"  Found {total_code_funcs} functions across {len(DOC_MODULES)} documented modules")

    # Step 2: Build documentation inventory
    print("\n[2/3] Parsing documentation...")
    doc_inventory = defaultdict(dict)
    all_doc_files = []

    for mod_name in DOC_MODULES:
        doc_path = DOCS_DIR / f"{mod_name}.md"
        if doc_path.exists():
            funcs = extract_documented_functions(doc_path)
            doc_inventory[mod_name].update(funcs)
            all_doc_files.append(doc_path)

    print(
        f"  Found {sum(len(v) for v in doc_inventory.values())} documented functions across {len(doc_inventory)} modules"
    )

    # Filter out class names (documented classes are not functions)
    for mod, funcs in doc_inventory.items():
        for fname in list(funcs.keys()):
            if fname in all_class_names:
                del funcs[fname]

    # Step 3: Compare and generate report
    print("\n[3/3] Comparing inventories...")
    results = {}
    total_mismatches = 0
    all_issues = []

    all_modules = DOC_MODULES  # Only check documented modules

    for module in all_modules:
        code_funcs = code_inventory.get(module, {})
        doc_funcs = doc_inventory.get(module, {})

        code_names = set(code_funcs.keys())
        doc_names = set(doc_funcs.keys())

        # Metrics
        missing_in_docs = code_names - doc_names
        extra_in_docs = doc_names - code_names
        matches = code_names & doc_names

        # Check signature mismatches
        mismatches = []
        for fname in matches:
            code_sig = code_funcs[fname]["signature"]
            doc_sig = doc_funcs[fname]["signature"]
            match, diffs = compare_signatures(code_sig, doc_sig)
            if not match:
                mismatches.append((fname, code_sig, doc_sig, diffs))

        total_mod_issues = len(missing_in_docs) + len(extra_in_docs) + len(mismatches)
        total_mismatches += total_mod_issues

        # Calculate accuracy (functions that match exactly)
        perfect_matches = len(matches) - len(mismatches)
        total_functions = len(code_names)  # Base on code
        accuracy = (perfect_matches / total_functions * 100) if total_functions > 0 else 100.0

        results[module] = {
            "accuracy": accuracy,
            "code_count": len(code_funcs),
            "doc_count": len(doc_funcs),
            "missing_in_docs": sorted(missing_in_docs),
            "extra_in_docs": sorted(extra_in_docs),
            "signature_mismatches": mismatches,
            "total_issues": total_mod_issues,
        }

        # Collect individual issues for priority ranking
        for fname in missing_in_docs:
            all_issues.append(
                (
                    f"{module}.{fname}",
                    "MISSING_IN_DOCS",
                    code_funcs[fname]["rel_path"],
                    code_funcs[fname]["line_no"],
                    code_funcs[fname]["signature"],
                    None,
                    [],
                )
            )

        for fname in extra_in_docs:
            all_issues.append(
                (
                    f"{module}.{fname}",
                    "EXTRA_IN_DOCS",
                    f"docs/core/{module}.md",
                    doc_funcs[fname]["line_no"],
                    None,
                    doc_funcs[fname]["signature"],
                    [],
                )
            )

        for fname, code_sig, doc_sig, diffs in mismatches:
            all_issues.append(
                (
                    f"{module}.{fname}",
                    "SIGNATURE_MISMATCH",
                    code_funcs[fname]["rel_path"],
                    code_funcs[fname]["line_no"],
                    code_sig,
                    doc_sig,
                    diffs,
                )
            )

    # Sort issues by priority: missing_in_docs > signature_mismatch > extra_in_docs
    # Then by line number (older/earlier = higher priority to fix)
    priority_order = {"MISSING_IN_DOCS": 0, "SIGNATURE_MISMATCH": 1, "EXTRA_IN_DOCS": 2}
    all_issues.sort(key=lambda x: (priority_order.get(x[1], 3), -x[3] if x[3] else 0, x[0]))

    # Generate console output
    print("\n" + "=" * 80)
    print("SUMMARY REPORT")
    print("=" * 80)
    print(f"\nTotal functions in code: {sum(len(v) for v in code_inventory.values())}")
    print(f"Total functions in docs: {sum(len(v) for v in doc_inventory.values())}")
    print(f"Total mismatches: {total_mismatches}")
    print("\nPer-module accuracy:")
    print(f"  {'Module':<15} {'Accuracy':>10} {'Code':>6} {'Docs':>6} {'Issues':>7}")
    print(f"  {'-'*15} {'-'*10} {'-'*6} {'-'*6} {'-'*7}")
    for module in sorted(results.keys()):
        r = results[module]
        print(f"  {module:<15} {r['accuracy']:>9.1f}% {r['code_count']:>6} {r['doc_count']:>6} {r['total_issues']:>7}")

    print("\n" + "=" * 80)
    print("TOP 10 PRIORITY ISSUES")
    print("=" * 80)
    print(f"\n{'Type':<22} {'File:Line':<35} {'Function':<30}")
    print(f"{'-'*22} {'-'*35} {'-'*30}")

    for i, (fqname, issue_type, filepath, lineno, code_sig, doc_sig, diffs) in enumerate(all_issues[:10]):
        file_line = f"{filepath}:{lineno}" if lineno else filepath
        print(f"{issue_type:<22} {file_line:<35} {fqname:<30}")
        if code_sig and doc_sig:
            print(f"  Code:    {code_sig}")
            print(f"  Doc:     {doc_sig}")
            if diffs:
                print(f"  Diff:    {'; '.join(diffs)}")
        elif code_sig:
            print(f"  Code:    {code_sig}")
        elif doc_sig:
            print(f"  Doc:     {doc_sig}")

    # Save detailed report
    report_path = DOCS_DIR / "CROSS_CHECK_REPORT.md"
    print(f"\n\nDetailed report saved to: {report_path}")

    report_lines = []
    report_lines.append("# Core Documentation Cross-Check Report\n")
    report_lines.append(f"Generated: {os.popen('date').read().strip()}\n")
    report_lines.append(f"**Total mismatches:** {total_mismatches}\n")
    report_lines.append(f"**Total functions in code:** {sum(len(v) for v in code_inventory.values())}\n")
    report_lines.append(f"**Total functions in docs:** {sum(len(v) for v in doc_inventory.values())}\n\n")

    report_lines.append("## Per-Module Accuracy\n\n")
    report_lines.append("| Module | Accuracy | Code Funcs | Doc Funcs | Issues |\n")
    report_lines.append("|--------|----------|------------|-----------|--------|\n")
    for module in sorted(results.keys()):
        r = results[module]
        report_lines.append(
            f"| {module} | {r['accuracy']:.1f}% | {r['code_count']} | {r['doc_count']} | {r['total_issues']} |\n"
        )

    report_lines.append("\n## Detailed Issues\n\n")

    for idx, (fqname, issue_type, filepath, lineno, code_sig, doc_sig, diffs) in enumerate(all_issues, 1):
        report_lines.append(f"### {idx}. {issue_type}: `{fqname}`\n\n")
        report_lines.append(f"- **Location:** `{filepath}:{lineno}`\n")
        if code_sig:
            report_lines.append(f"- **Code signature:** `{code_sig}`\n")
        if doc_sig:
            report_lines.append(f"- **Doc signature:** `{doc_sig}`\n")
        if diffs:
            report_lines.append(f"- **Differences:** {'; '.join(diffs)}\n")
        report_lines.append("\n")

    with open(report_path, "w") as f:
        f.writelines(report_lines)

    print("\nCross-check complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
