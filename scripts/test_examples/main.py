#!/usr/bin/env python3
"""Main orchestration module for METAINFORMANT test runner."""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

from .discover import discover_examples
from .reporting import generate_html_report, generate_junit_xml, print_summary, save_json_report
from .runner import run_examples_parallel, run_examples_sequential


class ExampleTester:
    """Test runner for METAINFORMANT examples."""

    def __init__(self, verbose: bool = False, domain_filter: str | None = None, continue_on_error: bool = False, parallel: bool = False, max_workers: int | None = None):
        self.verbose = verbose
        self.domain_filter = domain_filter
        self.continue_on_error = continue_on_error
        self.parallel = parallel
        self.max_workers = max_workers or min(4, os.cpu_count() or 1)  # Default to 4 or CPU count
        self.results: List[Dict[str, Any]] = []
        self.examples_dir = Path("examples")
        self.repo_root = Path(__file__).parent.parent.parent

    def run_all_examples(self) -> Dict[str, Any]:
        """Run all discovered examples and return results."""
        examples = discover_examples(self.examples_dir, self.domain_filter)

        if not examples:
            return {
                "summary": {"total": 0, "passed": 0, "failed": 0, "error": 0, "timeout": 0},
                "results": [],
                "error": f"No examples found{' for domain ' + self.domain_filter if self.domain_filter else ''}"
            }

        print(f"Found {len(examples)} examples to test{' for domain ' + self.domain_filter if self.domain_filter else ''}")
        if self.parallel:
            print(f"Running in parallel with {self.max_workers} workers")
        print("=" * 60)

        if self.parallel and len(examples) > 1:
            # Parallel execution
            self.results = run_examples_parallel(examples, self.repo_root, self.verbose, self.max_workers)
        else:
            # Sequential execution
            self.results = run_examples_sequential(examples, self.repo_root, self.verbose)

        # Generate summary
        summary = {
            "total": len(self.results),
            "passed": sum(1 for r in self.results if r["status"] == "passed"),
            "failed": sum(1 for r in self.results if r["status"] == "failed"),
            "error": sum(1 for r in self.results if r["status"] == "error"),
            "timeout": sum(1 for r in self.results if r["status"] == "timeout"),
            "total_execution_time": sum(r.get("execution_time", 0) for r in self.results if r.get("execution_time"))
        }

        summary["success_rate"] = summary["passed"] / summary["total"] if summary["total"] > 0 else 0

        return {
            "summary": summary,
            "results": self.results,
            "timestamp": time.time(),
            "command_line": sys.argv,
            "parallel_execution": self.parallel,
            "max_workers": self.max_workers
        }


def main():
    """Main function."""
    import argparse

    # Import os here to avoid circular imports
    import os

    parser = argparse.ArgumentParser(description="Test METAINFORMANT examples")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    parser.add_argument("--domain", "-d", help="Run examples for specific domain only")
    parser.add_argument("--continue-on-error", "-c", action="store_true", help="Continue testing even if examples fail")
    parser.add_argument("--junit-xml", action="store_true", help="Generate JUnit XML report for CI/CD systems")
    parser.add_argument("--html", action="store_true", help="Generate HTML test report for visualization")
    parser.add_argument("--parallel", "-p", action="store_true", help="Run examples in parallel for faster execution")
    parser.add_argument("--max-workers", type=int, help="Maximum number of parallel workers (default: 4)")

    args = parser.parse_args()

    # Create tester
    tester = ExampleTester(
        verbose=args.verbose,
        domain_filter=args.domain,
        continue_on_error=args.continue_on_error,
        parallel=args.parallel,
        max_workers=args.max_workers
    )

    try:
        # Run all examples
        results = tester.run_all_examples()

        # Save report
        output_dir = Path("output/examples")
        output_dir.mkdir(parents=True, exist_ok=True)

        save_json_report(results, output_dir)

        if args.junit_xml:
            generate_junit_xml(results, output_dir)

        if args.html:
            generate_html_report(results, output_dir)

        print_summary(results)

        # Exit with appropriate code
        summary = results["summary"]
        if summary["passed"] == summary["total"]:
            print("🎉 All examples passed!")
            sys.exit(0)
        else:
            print(f"❌ {summary['total'] - summary['passed']} examples failed")
            sys.exit(1)

    except Exception as e:
        print(f"Test runner failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
