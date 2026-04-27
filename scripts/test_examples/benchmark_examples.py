#!/usr/bin/env python3
"""Performance benchmarking script for METAINFORMANT examples.

This script establishes performance baselines and detects regressions in example execution times.

Usage:
    python scripts/benchmark_examples.py [--baseline] [--compare] [--threshold FLOAT] [--output DIR]

Arguments:
    --baseline: Establish new performance baseline
    --compare: Compare current performance against baseline
    --threshold: Performance regression threshold (default: 0.1 for 10% increase)
    --output: Output directory for benchmark results (default: output/examples/benchmarks)
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List


class ExampleBenchmarker:
    """Performance benchmarker for METAINFORMANT examples."""

    def __init__(self, output_dir: Path | None = None, runs: int = 3):
        self.output_dir = output_dir or Path("output/examples/benchmarks")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.baseline_file = self.output_dir / "baseline.json"
        self.runs = runs

    def establish_baseline(self) -> Dict[str, Any]:
        """Run examples multiple times and establish performance baseline."""
        print(f"Establishing performance baseline with {self.runs} runs per example...")

        all_results = []

        # Run examples multiple times using subprocess
        for run in range(self.runs):
            print(f"\nRun {run + 1}/{self.runs}")

            # Run test_examples.py script and capture output
            env = os.environ.copy()
            env["PYTHONPATH"] = str(Path.cwd() / "src")

            cmd = [sys.executable, "scripts/test_examples.py", "--continue-on-error", "--junit-xml", "--domain", "core"]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd(), env=env)

            if result.returncode != 0:
                print(f"❌ Test run {run + 1} failed: {result.stderr}")
                continue

            # Load results from the JSON output
            results_file = Path("output/examples/test_results.json")
            if results_file.exists():
                with open(results_file, "r") as f:
                    run_results = json.load(f)
                all_results.append(run_results)
            else:
                print(f"❌ No results file found for run {run + 1}")
                continue

        if not all_results:
            raise RuntimeError("No successful test runs completed")

        # Calculate baselines
        baseline = self._calculate_baseline(all_results)

        # Save baseline
        with open(self.baseline_file, "w") as f:
            json.dump(baseline, f, indent=2)

        print(f"Baseline established and saved to: {self.baseline_file}")
        return baseline

    def compare_performance(self, threshold: float = 0.1) -> Dict[str, Any]:
        """Compare current performance against baseline."""
        if not self.baseline_file.exists():
            print("❌ No baseline found. Run with --baseline first.")
            return {"error": "no_baseline"}

        print("Comparing current performance against baseline...")

        # Load baseline
        with open(self.baseline_file) as f:
            baseline = json.load(f)

        # Run current test using subprocess
        env = os.environ.copy()
        env["PYTHONPATH"] = str(Path.cwd() / "src")

        cmd = [sys.executable, "scripts/test_examples.py", "--continue-on-error", "--junit-xml"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd(), env=env)

        if result.returncode != 0:
            print(f"❌ Current test run failed: {result.stderr}")
            return {"error": "test_run_failed"}

        # Load current results
        results_file = Path("output/examples/test_results.json")
        if not results_file.exists():
            print("❌ No current results file found")
            return {"error": "no_results_file"}

        with open(results_file, "r") as f:
            current_results = json.load(f)

        # Compare results
        comparison = self._compare_results(baseline, current_results, threshold)

        # Save comparison
        comparison_file = self.output_dir / "comparison.json"
        with open(comparison_file, "w") as f:
            json.dump(comparison, f, indent=2)

        print(f"Comparison saved to: {comparison_file}")

        # Print summary
        self._print_comparison_summary(comparison)

        return comparison

    def _calculate_baseline(self, all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Calculate performance baseline from multiple runs."""
        baseline = {"timestamp": time.time(), "runs": self.runs, "examples": {}}

        # Group results by example
        example_results = {}
        for run_results in all_results:
            for result in run_results["results"]:
                example_path = result["example_path"]
                if example_path not in example_results:
                    example_results[example_path] = []
                example_results[example_path].append(result["execution_time"])

        # Calculate statistics for each example
        for example_path, times in example_results.items():
            baseline["examples"][example_path] = {
                "mean": statistics.mean(times),
                "median": statistics.median(times),
                "stdev": statistics.stdev(times) if len(times) > 1 else 0,
                "min": min(times),
                "max": max(times),
                "runs": len(times),
            }

        # Overall statistics
        all_times = [result["execution_time"] for run_results in all_results for result in run_results["results"]]
        baseline["overall"] = {
            "total_mean": statistics.mean(all_times),
            "total_median": statistics.median(all_times),
            "total_stdev": statistics.stdev(all_times) if len(all_times) > 1 else 0,
            "total_examples": len(example_results),
        }

        return baseline

    def _compare_results(self, baseline: Dict[str, Any], current: Dict[str, Any], threshold: float) -> Dict[str, Any]:
        """Compare current results against baseline."""
        comparison = {
            "timestamp": time.time(),
            "baseline_timestamp": baseline["timestamp"],
            "threshold": threshold,
            "examples": {},
            "regressions": [],
            "improvements": [],
        }

        # Compare each example
        for result in current["results"]:
            example_path = result["example_path"]
            current_time = result["execution_time"]

            if example_path in baseline["examples"]:
                baseline_stats = baseline["examples"][example_path]
                baseline_mean = baseline_stats["mean"]

                # Calculate change
                change = (current_time - baseline_mean) / baseline_mean
                is_regression = change > threshold
                is_improvement = change < -threshold

                comparison["examples"][example_path] = {
                    "current_time": current_time,
                    "baseline_mean": baseline_mean,
                    "change_percent": change * 100,
                    "is_regression": is_regression,
                    "is_improvement": is_improvement,
                    "status": "passed" if result["status"] == "passed" else "failed",
                }

                if is_regression:
                    comparison["regressions"].append(
                        {
                            "example": example_path,
                            "change_percent": change * 100,
                            "current_time": current_time,
                            "baseline_mean": baseline_mean,
                        }
                    )
                elif is_improvement:
                    comparison["improvements"].append(
                        {
                            "example": example_path,
                            "change_percent": change * 100,
                            "current_time": current_time,
                            "baseline_mean": baseline_mean,
                        }
                    )
            else:
                # New example
                comparison["examples"][example_path] = {
                    "current_time": current_time,
                    "baseline_mean": None,
                    "change_percent": None,
                    "is_regression": False,
                    "is_improvement": False,
                    "status": "passed" if result["status"] == "passed" else "failed",
                    "new_example": True,
                }

        return comparison

    def _print_comparison_summary(self, comparison: Dict[str, Any]) -> None:
        """Print a summary of the performance comparison."""
        regressions = len(comparison["regressions"])
        improvements = len(comparison["improvements"])
        threshold = comparison["threshold"] * 100

        print("\nPerformance Comparison Summary")
        print("=" * 50)
        print(f"Threshold: {threshold:.1f}%")
        print(f"Regressions (> +{threshold:.1f}%): {regressions}")
        print(f"Improvements (< -{threshold:.1f}%): {improvements}")

        if regressions > 0:
            print(f"\n❌ Regressions detected:")
            for reg in comparison["regressions"][:5]:  # Show top 5
                print(".1f")

            if len(comparison["regressions"]) > 5:
                print(f"  ... and {len(comparison['regressions']) - 5} more")

        if improvements > 0:
            print(f"\n✅ Performance improvements:")
            for imp in comparison["improvements"][:5]:  # Show top 5
                print(".1f")

        if regressions == 0 and improvements == 0:
            print("\n✅ No significant performance changes detected")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Benchmark METAINFORMANT examples performance")
    parser.add_argument("--baseline", action="store_true", help="Establish new performance baseline")
    parser.add_argument("--compare", action="store_true", help="Compare current performance against baseline")
    parser.add_argument(
        "--threshold", type=float, default=0.1, help="Performance regression threshold (default: 0.1 for 10%%)"
    )
    parser.add_argument("--output", type=Path, help="Output directory for benchmark results")
    parser.add_argument("--runs", type=int, default=3, help="Number of runs for baseline establishment (default: 3)")

    args = parser.parse_args()

    if not args.baseline and not args.compare:
        parser.error("Must specify either --baseline or --compare")

    # Create benchmarker
    benchmarker = ExampleBenchmarker(output_dir=args.output, runs=args.runs)

    try:
        if args.baseline:
            baseline = benchmarker.establish_baseline()
            print(f"✅ Baseline established for {len(baseline['examples'])} examples")
        elif args.compare:
            comparison = benchmarker.compare_performance(threshold=args.threshold)
            if comparison.get("error"):
                print(f"❌ {comparison['error']}")
                return 1

            regressions = len(comparison["regressions"])
            if regressions > 0:
                print(f"❌ {regressions} performance regressions detected")
                return 1
            else:
                print("✅ No performance regressions detected")
                return 0

    except Exception as e:
        print(f"❌ Benchmarking failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
