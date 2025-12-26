#!/usr/bin/env python3
"""Fast validation script for METAINFORMANT examples.

This script performs quick validation checks on examples without running them,
useful for pre-commit hooks and CI/CD validation.

Usage:
    python scripts/validate_examples.py [--fast] [--verbose] [--fix]

Arguments:
    --fast: Only perform syntax and import checks (default)
    --verbose: Enable verbose output
    --fix: Attempt to auto-fix simple issues
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import sys
from pathlib import Path
from typing import Any, Dict, List, Set


class ExampleValidator:
    """Fast validator for METAINFORMANT examples."""

    def __init__(self, verbose: bool = False, fix: bool = False):
        self.verbose = verbose
        self.fix = fix
        self.examples_dir = Path("examples")
        self.issues: List[Dict[str, Any]] = []

    def validate_all_examples(self) -> Dict[str, Any]:
        """Validate all examples and return results."""
        if not self.examples_dir.exists():
            return {"error": f"Examples directory not found: {self.examples_dir}"}

        all_examples = list(self.examples_dir.rglob("example_*.py"))

        for example_path in all_examples:
            if self.verbose:
                print(f"Validating {example_path.relative_to(self.examples_dir)}...")
            self.validate_example(example_path)

        # Categorize issues
        errors = [i for i in self.issues if i["severity"] == "error"]
        warnings = [i for i in self.issues if i["severity"] == "warning"]

        return {
            "total_examples": len(all_examples),
            "errors": len(errors),
            "warnings": len(warnings),
            "issues": self.issues
        }

    def validate_example(self, example_path: Path) -> None:
        """Validate a single example file."""
        try:
            # 1. Syntax check
            self._check_syntax(example_path)

            # 2. Import validation
            self._check_imports(example_path)

            # 3. Basic structure validation
            self._check_structure(example_path)

        except Exception as e:
            self._add_issue(example_path, "error", "validation_failed", f"Validation failed: {e}")

    def _check_syntax(self, example_path: Path) -> None:
        """Check Python syntax."""
        try:
            with open(example_path, 'r', encoding='utf-8') as f:
                ast.parse(f.read())
        except SyntaxError as e:
            self._add_issue(example_path, "error", "syntax_error",
                          f"Syntax error at line {e.lineno}: {e.msg}")
        except Exception as e:
            self._add_issue(example_path, "error", "syntax_check_failed",
                          f"Syntax check failed: {e}")

    def _check_imports(self, example_path: Path) -> None:
        """Check import statements for common issues."""
        try:
            with open(example_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Check for problematic imports
            if "from metainformant." in content:
                # This is expected - check if imports are valid
                self._validate_metainformant_imports(content, example_path)

        except Exception as e:
            self._add_issue(example_path, "warning", "import_check_failed",
                          f"Import check failed: {e}")

    def _validate_metainformant_imports(self, content: str, example_path: Path) -> None:
        """Validate METAINFORMANT imports exist."""
        import re

        # Extract import statements
        imports = re.findall(r'from metainformant\.([^\'"\s]+)', content)

        for module_path in imports:
            # Convert module path to file path
            module_parts = module_path.split('.')
            file_path = Path("src/metainformant") / Path(*module_parts)

            # Check if it's a Python file or package
            if not (file_path.with_suffix('.py').exists() or (file_path / '__init__.py').exists()):
                self._add_issue(example_path, "warning", "invalid_import",
                              f"Import 'metainformant.{module_path}' may not exist")

    def _check_structure(self, example_path: Path) -> None:
        """Check basic example structure."""
        try:
            with open(example_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Check for required elements
            if 'if __name__ == "__main__":' not in content:
                self._add_issue(example_path, "warning", "missing_main_guard",
                              "Missing 'if __name__ == \"__main__\":' guard")

            if 'def main():' not in content:
                self._add_issue(example_path, "info", "missing_main_function",
                              "Consider adding a main() function")

            # Check for output directory usage
            if 'output/examples/' not in content:
                self._add_issue(example_path, "info", "no_output_directory",
                              "Example doesn't write to output/examples/")

        except Exception as e:
            self._add_issue(example_path, "warning", "structure_check_failed",
                          f"Structure check failed: {e}")

    def _add_issue(self, example_path: Path, severity: str, issue_type: str, message: str) -> None:
        """Add an issue to the issues list."""
        self.issues.append({
            "example": str(example_path.relative_to(self.examples_dir)),
            "severity": severity,
            "type": issue_type,
            "message": message,
            "file": str(example_path)
        })


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Validate METAINFORMANT examples")
    parser.add_argument("--fast", action="store_true", help="Fast validation only (default)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    parser.add_argument("--fix", action="store_true", help="Attempt to auto-fix simple issues")

    args = parser.parse_args()

    # Create validator
    validator = ExampleValidator(verbose=args.verbose, fix=args.fix)

    try:
        # Validate all examples
        results = validator.validate_all_examples()

        # Print results
        if results.get("error"):
            print(f"❌ {results['error']}")
            sys.exit(1)

        print(f"Validated {results['total_examples']} examples")
        print(f"Errors: {results['errors']}")
        print(f"Warnings: {results['warnings']}")

        # Print issues
        for issue in results["issues"]:
            if issue["severity"] == "error":
                print(f"❌ {issue['example']}: {issue['message']}")
            elif issue["severity"] == "warning":
                print(f"⚠️  {issue['example']}: {issue['message']}")
            elif args.verbose and issue["severity"] == "info":
                print(f"ℹ️  {issue['example']}: {issue['message']}")

        # Exit with appropriate code
        if results["errors"] > 0:
            print("❌ Validation failed due to errors")
            sys.exit(1)
        elif results["warnings"] > 0:
            print("⚠️  Validation passed with warnings")
            sys.exit(0)
        else:
            print("✅ All examples validated successfully")
            sys.exit(0)

    except Exception as e:
        print(f"❌ Validation failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
