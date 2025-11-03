#!/usr/bin/env python3
"""Script to audit docstring coverage in source files."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any


def audit_file(file_path: Path) -> dict[str, Any]:
    """Audit a Python file for missing docstrings.
    
    Args:
        file_path: Path to Python file to audit
        
    Returns:
        Dictionary with audit results
    """
    results = {
        "file": str(file_path),
        "functions": [],
        "classes": [],
        "missing_docstrings": [],
    }
    
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
            tree = ast.parse(content, filename=str(file_path))
    except Exception as e:
        results["error"] = str(e)
        return results
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            # Skip private methods (starting with _) unless they're __init__ or __str__ etc
            if node.name.startswith("_") and not node.name.startswith("__"):
                continue
            
            has_docstring = (
                ast.get_docstring(node) is not None
                or (node.body and isinstance(node.body[0], ast.Expr) 
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str))
            )
            
            results["functions"].append({
                "name": node.name,
                "line": node.lineno,
                "has_docstring": has_docstring,
            })
            
            if not has_docstring:
                results["missing_docstrings"].append({
                    "type": "function",
                    "name": node.name,
                    "line": node.lineno,
                })
        
        elif isinstance(node, ast.ClassDef):
            has_docstring = ast.get_docstring(node) is not None
            
            results["classes"].append({
                "name": node.name,
                "line": node.lineno,
                "has_docstring": has_docstring,
            })
            
            if not has_docstring:
                results["missing_docstrings"].append({
                    "type": "class",
                    "name": node.name,
                    "line": node.lineno,
                })
    
    return results


def main():
    """Main audit function."""
    src_dir = Path("src/metainformant")
    
    all_results = []
    total_missing = 0
    
    for py_file in src_dir.rglob("*.py"):
        if py_file.name.startswith("__") and py_file.name != "__init__.py":
            continue
        
        results = audit_file(py_file)
        if "error" not in results:
            missing_count = len(results["missing_docstrings"])
            if missing_count > 0:
                all_results.append(results)
                total_missing += missing_count
    
    # Print summary
    print(f"Found {total_missing} missing docstrings across {len(all_results)} files\n")
    
    for result in sorted(all_results, key=lambda x: len(x["missing_docstrings"]), reverse=True):
        if result["missing_docstrings"]:
            print(f"\n{result['file']}:")
            for item in result["missing_docstrings"]:
                print(f"  {item['type']}: {item['name']} (line {item['line']})")
    
    # Write detailed report
    report_path = Path("output/audit_docstring_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    
    import json
    with open(report_path, "w") as f:
        json.dump({
            "summary": {
                "total_files_with_missing": len(all_results),
                "total_missing_docstrings": total_missing,
            },
            "details": all_results,
        }, f, indent=2)
    
    print(f"\n\nDetailed report written to {report_path}")


if __name__ == "__main__":
    main()

