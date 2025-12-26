#!/usr/bin/env python3
"""Analyze coverage of METAINFORMANT functions demonstrated in examples.

This script analyzes which METAINFORMANT functions are demonstrated by the examples,
providing insights into test coverage and identifying gaps in example coverage.

Usage:
    python scripts/analyze_example_coverage.py [--output DIR] [--detailed] [--report]

Arguments:
    --output: Output directory for coverage reports (default: output/examples/coverage)
    --detailed: Generate detailed per-function coverage analysis
    --report: Generate human-readable coverage report
"""

from __future__ import annotations

import argparse
import ast
import importlib
import inspect
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Set


class CoverageAnalyzer:
    """Analyze coverage of METAINFORMANT functions in examples."""

    def __init__(self, output_dir: Path | None = None):
        self.output_dir = output_dir or Path("output/examples/coverage")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Core modules to analyze
        self.core_modules = [
            "metainformant.core",
            "metainformant.dna",
            "metainformant.rna",
            "metainformant.gwas",
            "metainformant.protein",
            "metainformant.epigenome",
            "metainformant.ontology",
            "metainformant.phenotype",
            "metainformant.ecology",
            "metainformant.math",
            "metainformant.information",
            "metainformant.life_events",
            "metainformant.multiomics",
            "metainformant.singlecell",
            "metainformant.quality",
            "metainformant.networks",
            "metainformant.ml",
            "metainformant.simulation",
            "metainformant.visualization"
        ]

    def analyze_coverage(self) -> Dict[str, Any]:
        """Analyze coverage of METAINFORMANT functions in examples."""
        print("Analyzing METAINFORMANT function coverage in examples...")

        # Discover all functions in METAINFORMANT
        all_functions = self._discover_functions()

        # Analyze which functions are used in examples
        examples_dir = Path("examples")
        covered_functions = self._analyze_examples(examples_dir, all_functions)

        # Calculate coverage statistics
        coverage_stats = self._calculate_coverage_stats(all_functions, covered_functions)

        # Generate reports
        coverage_data = {
            "timestamp": __import__('time').time(),
            "total_functions": len(all_functions),
            "covered_functions": len(covered_functions),
            "coverage_percentage": coverage_stats["overall_coverage"],
            "module_coverage": coverage_stats["module_coverage"],
            "uncovered_functions": coverage_stats["uncovered_functions"],
            "covered_functions_detail": covered_functions,
            "examples_analyzed": len(list(examples_dir.rglob("example_*.py")))
        }

        return coverage_data

    def _discover_functions(self) -> Dict[str, Dict[str, Any]]:
        """Discover all public functions in METAINFORMANT modules."""
        functions = {}

        for module_name in self.core_modules:
            try:
                module = importlib.import_module(module_name)
                module_functions = self._get_module_functions(module, module_name)
                functions.update(module_functions)
            except ImportError:
                print(f"Warning: Could not import {module_name}")
                continue

        return functions

    def _get_module_functions(self, module, module_name: str) -> Dict[str, Dict[str, Any]]:
        """Get all public functions from a module."""
        functions = {}

        for name, obj in inspect.getmembers(module):
            if (inspect.isfunction(obj) or inspect.ismethod(obj)) and not name.startswith('_'):
                # Skip private functions and methods
                if hasattr(obj, '__module__') and obj.__module__ == module_name:
                    functions[f"{module_name}.{name}"] = {
                        "name": name,
                        "module": module_name,
                        "signature": str(inspect.signature(obj)),
                        "docstring": obj.__doc__ or "",
                        "file": getattr(obj, '__code__', None) and obj.__code__.co_filename
                    }

        return functions

    def _analyze_examples(self, examples_dir: Path, all_functions: Dict[str, Dict[str, Any]]) -> Dict[str, List[str]]:
        """Analyze which functions are used in examples."""
        covered_functions = {}

        # Get function names for pattern matching
        function_names = {name.split('.')[-1]: full_name for full_name, _ in all_functions.items()}

        for example_file in examples_dir.rglob("example_*.py"):
            try:
                with open(example_file, 'r', encoding='utf-8') as f:
                    content = f.read()

                # Extract function calls
                called_functions = self._extract_function_calls(content)

                # Map to full function names
                for func_name in called_functions:
                    if func_name in function_names:
                        full_name = function_names[func_name]
                        if full_name not in covered_functions:
                            covered_functions[full_name] = []
                        covered_functions[full_name].append(str(example_file.relative_to(examples_dir)))

            except Exception as e:
                print(f"Warning: Could not analyze {example_file}: {e}")
                continue

        return covered_functions

    def _extract_function_calls(self, content: str) -> Set[str]:
        """Extract function calls from Python code."""
        called_functions = set()

        try:
            tree = ast.parse(content)

            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    if isinstance(node.func, ast.Name):
                        # Direct function call
                        called_functions.add(node.func.id)
                    elif isinstance(node.func, ast.Attribute):
                        # Method call
                        called_functions.add(node.func.attr)

        except SyntaxError:
            # If parsing fails, use regex as fallback
            func_pattern = r'\b([a-zA-Z_][a-zA-Z0-9_]*)\s*\('
            matches = re.findall(func_pattern, content)
            called_functions.update(matches)

        return called_functions

    def _calculate_coverage_stats(self, all_functions: Dict[str, Dict[str, Any]], covered_functions: Dict[str, List[str]]) -> Dict[str, Any]:
        """Calculate coverage statistics."""
        total_functions = len(all_functions)
        covered_count = len(covered_functions)

        # Module-level coverage
        module_coverage = {}
        uncovered_functions = {}

        for full_name, func_info in all_functions.items():
            module = func_info["module"]
            if module not in module_coverage:
                module_coverage[module] = {"total": 0, "covered": 0, "functions": []}

            module_coverage[module]["total"] += 1
            module_coverage[module]["functions"].append(full_name)

            if full_name in covered_functions:
                module_coverage[module]["covered"] += 1
            else:
                uncovered_functions[full_name] = func_info

        # Calculate percentages
        for module_stats in module_coverage.values():
            module_stats["coverage_percentage"] = (
                module_stats["covered"] / module_stats["total"] * 100
                if module_stats["total"] > 0 else 0
            )

        overall_coverage = covered_count / total_functions * 100 if total_functions > 0 else 0

        return {
            "overall_coverage": overall_coverage,
            "module_coverage": module_coverage,
            "uncovered_functions": uncovered_functions
        }

    def save_coverage_report(self, coverage_data: Dict[str, Any], detailed: bool = False, report: bool = False) -> None:
        """Save coverage analysis results."""
        # Save JSON data
        coverage_file = self.output_dir / "coverage.json"
        with open(coverage_file, 'w') as f:
            json.dump(coverage_data, f, indent=2)

        print(f"Coverage analysis saved to: {coverage_file}")

        # Generate detailed report if requested
        if detailed:
            detailed_file = self.output_dir / "coverage_detailed.json"
            with open(detailed_file, 'w') as f:
                json.dump({
                    "covered_functions": coverage_data["covered_functions_detail"],
                    "uncovered_functions": coverage_data["uncovered_functions"]
                }, f, indent=2)
            print(f"Detailed coverage saved to: {detailed_file}")

        # Generate human-readable report
        if report:
            report_file = self.output_dir / "coverage_report.md"
            self._generate_readable_report(coverage_data, report_file)
            print(f"Coverage report saved to: {report_file}")

    def _generate_readable_report(self, coverage_data: Dict[str, Any], report_file: Path) -> None:
        """Generate a human-readable coverage report."""
        with open(report_file, 'w') as f:
            f.write("# METAINFORMANT Examples Coverage Report\n\n")

            f.write(f"## Overview\n\n")
            f.write(f"- **Total Functions:** {coverage_data['total_functions']}\n")
            f.write(f"- **Covered Functions:** {coverage_data['covered_functions']}\n")
            f.write(".1f")
            f.write(f"- **Examples Analyzed:** {coverage_data['examples_analyzed']}\n\n")

            f.write("## Module Coverage\n\n")
            f.write("| Module | Functions | Covered | Coverage |\n")
            f.write("|--------|-----------|---------|----------|\n")

            for module, stats in sorted(coverage_data['module_coverage'].items()):
                module_short = module.replace('metainformant.', '')
                f.write(".1f")

            f.write("\n## Top Uncovered Functions\n\n")
            uncovered = coverage_data['uncovered_functions']
            if uncovered:
                f.write("| Function | Module |\n")
                f.write("|----------|--------|\n")

                # Show top 20 uncovered functions
                for i, (func_name, func_info) in enumerate(list(uncovered.items())[:20]):
                    module_short = func_info['module'].replace('metainformant.', '')
                    func_short = func_name.split('.')[-1]
                    f.write(f"| `{func_short}` | {module_short} |\n")

                if len(uncovered) > 20:
                    f.write(f"| ... and {len(uncovered) - 20} more | |\n")
            else:
                f.write("All functions are covered! 🎉\n")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Analyze METAINFORMANT examples coverage")
    parser.add_argument("--output", type=Path, help="Output directory for coverage reports")
    parser.add_argument("--detailed", action="store_true", help="Generate detailed per-function coverage analysis")
    parser.add_argument("--report", action="store_true", help="Generate human-readable coverage report")

    args = parser.parse_args()

    # Create analyzer
    analyzer = CoverageAnalyzer(output_dir=args.output)

    try:
        # Analyze coverage
        coverage_data = analyzer.analyze_coverage()

        # Save reports
        analyzer.save_coverage_report(coverage_data, detailed=args.detailed, report=args.report)

        # Print summary
        print("\nCoverage Analysis Summary")
        print("=" * 40)
        print(f"Total functions: {coverage_data['total_functions']}")
        print(f"Covered functions: {coverage_data['covered_functions']}")
        print(".1f")
        print(f"Examples analyzed: {coverage_data['examples_analyzed']}")

        if coverage_data['uncovered_functions']:
            print(f"\n⚠️  {len(coverage_data['uncovered_functions'])} functions not covered by examples")
        else:
            print("\n✅ All functions are covered by examples!")

    except Exception as e:
        print(f"❌ Coverage analysis failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
