#!/usr/bin/env python3
"""Test runner for all METAINFORMANT examples.

This script discovers and runs all example files in the examples/ directory headlessly,
validates their outputs, and generates a comprehensive test report.

Usage:
    python scripts/test_examples.py [--verbose] [--domain DOMAIN] [--continue-on-error]

Arguments:
    --verbose: Enable verbose output
    --domain: Run examples for specific domain only (e.g., dna, rna, gwas)
    --continue-on-error: Continue testing even if individual examples fail
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

# Set headless environment before any matplotlib imports
os.environ["MPLBACKEND"] = "Agg"

# Add src to path for imports
import sys
from pathlib import Path
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root / "src"))

from metainformant.core import io


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

    def discover_examples(self) -> List[Path]:
        """Discover all example Python files."""
        if not self.examples_dir.exists():
            raise FileNotFoundError(f"Examples directory not found: {self.examples_dir}")

        all_examples = list(self.examples_dir.rglob("example_*.py"))

        if self.domain_filter:
            # Filter examples by domain
            filtered = []
            for example in all_examples:
                domain = example.parent.name
                if domain == self.domain_filter:
                    filtered.append(example)
            all_examples = filtered

        return sorted(all_examples)

    def run_example(self, example_path: Path) -> Dict[str, Any]:
        """Run a single example and capture results."""
        start_time = time.time()

        result = {
            "example_path": str(example_path.relative_to(self.examples_dir)),
            "domain": example_path.parent.name,
            "filename": example_path.name,
            "exit_code": None,
            "execution_time": None,
            "stdout": "",
            "stderr": "",
            "output_files_created": [],
            "json_outputs_valid": [],
            "json_outputs_invalid": [],
            "status": "running",
            "error_message": None,
            "error_details": None,
            "error_line": None,
            "error_file": None,
            "traceback": None
        }

        try:
            if self.verbose:
                print(f"Running {result['example_path']}...")

            # Run the example
            env = {**os.environ, "MPLBACKEND": "Agg"}  # Ensure headless
            env["PYTHONPATH"] = str(repo_root / "src")  # Add src to Python path

            process = subprocess.run(
                [sys.executable, str(example_path)],
                cwd=Path.cwd(),
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
                env=env
            )

            result["exit_code"] = process.returncode
            result["stdout"] = process.stdout
            result["stderr"] = process.stderr
            result["execution_time"] = time.time() - start_time

            # Parse error details from stderr
            self._parse_error_details(result)

            # Validate outputs
            output_files_created, json_outputs_valid, json_outputs_invalid = self.validate_outputs(example_path)

            result["output_files_created"] = output_files_created
            result["json_outputs_valid"] = json_outputs_valid
            result["json_outputs_invalid"] = json_outputs_invalid

            # Determine status
            if process.returncode == 0:
                # Check if expected outputs were created
                expected_outputs = self.get_expected_outputs(example_path)
                if expected_outputs and not output_files_created:
                    result["status"] = "failed"
                    result["error_message"] = f"Expected output files not created: {expected_outputs}"
                else:
                    result["status"] = "passed"
            else:
                result["status"] = "failed"
                if not result["error_message"]:
                    result["error_message"] = f"Non-zero exit code: {process.returncode}"

        except subprocess.TimeoutExpired:
            result["status"] = "timeout"
            result["error_message"] = "Example timed out after 5 minutes"
            result["execution_time"] = time.time() - start_time
        except Exception as e:
            result["status"] = "error"
            result["error_message"] = str(e)
            result["traceback"] = self._format_exception(e)
            result["execution_time"] = time.time() - start_time

        return result

    def get_expected_outputs(self, example_path: Path) -> List[str]:
        """Extract expected output files from example docstring."""
        try:
            with open(example_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Look for "Output:" section in docstring
            lines = content.split('\n')
            in_output_section = False
            outputs = []

            for line in lines:
                if line.strip().startswith('Output:'):
                    in_output_section = True
                    continue
                elif in_output_section and line.strip().startswith('"""'):
                    break
                elif in_output_section and line.strip():
                    # Extract output paths
                    output_line = line.strip()
                    if output_line.startswith('output/'):
                        outputs.append(output_line)

            return outputs

        except Exception:
            return []

    def validate_outputs(self, example_path: Path) -> tuple[List[str], List[str], List[str]]:
        """Validate that expected output files were created."""
        output_dir = Path("output/examples") / example_path.parent.name

        if not output_dir.exists():
            return [], [], []

        # Find all files in the output directory (created during this run)
        created_files = []
        json_valid = []
        json_invalid = []

        for file_path in output_dir.rglob("*"):
            if file_path.is_file():
                relative_path = str(file_path.relative_to(Path("output/examples")))
                created_files.append(relative_path)

                # Validate JSON files
                if file_path.suffix.lower() == '.json':
                    validation_result = self._validate_json_content(file_path, example_path.parent.name, example_path.stem)
                    if validation_result["valid"]:
                        json_valid.append(relative_path)
                    else:
                        json_invalid.append(f"{relative_path} ({validation_result['error']})")

        return created_files, json_valid, json_invalid

    def _validate_json_content(self, json_path: Path, domain: str, example_name: str) -> Dict[str, Any]:
        """Validate JSON content structure and data quality."""
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            # Basic structure validation
            if not isinstance(data, dict):
                return {"valid": False, "error": "Root must be a JSON object"}

            # Check for required fields
            required_fields = ["example", "domain"]
            for field in required_fields:
                if field not in data:
                    return {"valid": False, "error": f"Missing required field: {field}"}

            # Validate domain consistency
            if data.get("domain") != domain:
                return {"valid": False, "error": f"Domain mismatch: expected {domain}, got {data.get('domain')}"}

            # Validate results field exists
            if "results" not in data:
                return {"valid": False, "error": "Missing results field"}

            # Domain-specific validation
            domain_validation = self._validate_domain_content(data, domain, example_name)
            if not domain_validation["valid"]:
                return domain_validation

            return {"valid": True, "error": None}

        except json.JSONDecodeError as e:
            return {"valid": False, "error": f"Invalid JSON: {e}"}
        except Exception as e:
            return {"valid": False, "error": f"Validation error: {e}"}

    def _validate_domain_content(self, data: Dict[str, Any], domain: str, example_name: str) -> Dict[str, Any]:
        """Validate domain-specific content."""
        results = data.get("results", {})

        try:
            if domain == "dna":
                return self._validate_dna_content(results, example_name)
            elif domain == "gwas":
                return self._validate_gwas_content(results, example_name)
            elif domain == "ml":
                return self._validate_ml_content(results, example_name)
            elif domain == "core":
                return self._validate_core_content(results, example_name)
            else:
                # Generic validation for other domains
                if not isinstance(results, (dict, list)):
                    return {"valid": False, "error": "Results must be an object or array"}
                return {"valid": True, "error": None}

        except Exception as e:
            return {"valid": False, "error": f"Domain validation error: {e}"}

    def _validate_dna_content(self, results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
        """Validate DNA analysis results."""
        if example_name == "example_sequences":
            # Check for sequence data
            if not any(isinstance(v, dict) and "sequence" in v for v in results.values()):
                return {"valid": False, "error": "Missing sequence data in results"}

        elif example_name == "example_alignment":
            # Check for alignment data
            if not any("alignment" in str(v).lower() for v in results.values()):
                return {"valid": False, "error": "Missing alignment data in results"}

        return {"valid": True, "error": None}

    def _validate_gwas_content(self, results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
        """Validate GWAS analysis results."""
        if example_name == "example_association":
            # Check for statistical results
            if "p_values" not in results and not any("p_value" in str(k).lower() for k in results.keys()):
                return {"valid": False, "error": "Missing p-values in GWAS results"}

        return {"valid": True, "error": None}

    def _validate_ml_content(self, results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
        """Validate ML analysis results."""
        if example_name == "example_pipeline":
            # Check for performance metrics
            if "accuracy" not in results and "performance" not in results:
                return {"valid": False, "error": "Missing performance metrics in ML results"}

        return {"valid": True, "error": None}

    def _validate_core_content(self, results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
        """Validate core utility results."""
        # Core examples should have basic validation
        if not results:
            return {"valid": False, "error": "Empty results in core example"}

        return {"valid": True, "error": None}

    def run_all_examples(self) -> Dict[str, Any]:
        """Run all discovered examples and return results."""
        examples = self.discover_examples()

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
            self._run_examples_parallel(examples)
        else:
            # Sequential execution
            self._run_examples_sequential(examples)

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

    def _run_examples_sequential(self, examples: List[Path]) -> None:
        """Run examples sequentially."""
        for example_path in examples:
            result = self.run_example(example_path)
            self.results.append(result)

            # Print status
            status_icon = {"passed": "✓", "failed": "✗", "error": "✗", "timeout": "⏱"}.get(result["status"], "?")
            print(f"{status_icon} {result['example_path']} - {result['status']}")

            if self.verbose and result["error_message"]:
                print(f"    Error: {result['error_message']}")

            if not self.continue_on_error and result["status"] != "passed":
                print("Stopping due to failure (use --continue-on-error to continue)")
                break

    def _run_examples_parallel(self, examples: List[Path]) -> None:
        """Run examples in parallel using thread pools."""
        # For parallel execution, we need to be careful about output conflicts
        # Each example writes to its own domain subdirectory, so this should be safe

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all jobs
            future_to_example = {
                executor.submit(self.run_example, example_path): example_path
                for example_path in examples
            }

            # Process results as they complete
            for future in as_completed(future_to_example):
                example_path = future_to_example[future]
                try:
                    result = future.result()
                    self.results.append(result)

                    # Print status
                    status_icon = {"passed": "✓", "failed": "✗", "error": "✗", "timeout": "⏱"}.get(result["status"], "?")
                    print(f"{status_icon} {result['example_path']} - {result['status']}")

                    if self.verbose and result["error_message"]:
                        print(f"    Error: {result['error_message']}")

                except Exception as exc:
                    print(f"❌ {example_path.relative_to(self.examples_dir)} - error: {exc}")
                    # Add error result
                    self.results.append({
                        "example_path": str(example_path.relative_to(self.examples_dir)),
                        "domain": example_path.parent.name,
                        "filename": example_path.name,
                        "status": "error",
                        "error_message": str(exc),
                        "execution_time": 0,
                        "exit_code": -1,
                        "stdout": "",
                        "stderr": str(exc),
                        "output_files_created": [],
                        "json_outputs_valid": [],
                        "json_outputs_invalid": []
                    })

    def save_report(self, results: Dict[str, Any], junit_xml: bool = False, html: bool = False) -> None:
        """Save test results to file."""
        output_dir = Path("output/examples")
        output_dir.mkdir(parents=True, exist_ok=True)

        report_file = output_dir / "test_results.json"
        io.dump_json(results, report_file, indent=2)

        print(f"\nTest report saved to: {report_file}")

        # Generate JUnit XML for CI/CD systems
        if junit_xml:
            self._generate_junit_xml(results)

        # Generate HTML report
        if html:
            self._generate_html_report(results)

        # Print summary
        summary = results["summary"]
        print("\n" + "=" * 60)
        print("TEST SUMMARY")
        print("=" * 60)
        print(f"Total examples: {summary['total']}")
        print(f"Passed: {summary['passed']}")
        print(f"Failed: {summary['failed']}")
        print(f"Errors: {summary['error']}")
        print(f"Timeouts: {summary['timeout']}")
        print(".1f")
        print(".1f")
        print("=" * 60)

    def _generate_junit_xml(self, results: Dict[str, Any]) -> None:
        """Generate JUnit XML report for CI/CD systems."""
        import xml.etree.ElementTree as ET
        from xml.dom import minidom

        output_dir = Path("output/examples")
        junit_file = output_dir / "junit_report.xml"

        # Create root testsuites element
        testsuites = ET.Element("testsuites")
        testsuites.set("name", "METAINFORMANT Examples")
        testsuites.set("tests", str(results["summary"]["total"]))
        testsuites.set("failures", str(results["summary"]["failed"] + results["summary"]["error"]))
        testsuites.set("time", str(results["summary"]["total_execution_time"]))

        # Create testsuite for examples
        testsuite = ET.SubElement(testsuites, "testsuite")
        testsuite.set("name", "examples")
        testsuite.set("tests", str(results["summary"]["total"]))
        testsuite.set("failures", str(results["summary"]["failed"] + results["summary"]["error"]))
        testsuite.set("time", str(results["summary"]["total_execution_time"]))

        # Add test cases
        for result in results["results"]:
            testcase = ET.SubElement(testsuite, "testcase")
            testcase.set("name", result["example_path"])
            testcase.set("classname", f"examples.{result['domain']}")
            testcase.set("time", str(result.get("execution_time", 0)))

            if result["status"] != "passed":
                failure = ET.SubElement(testcase, "failure")
                failure.set("message", result.get("error_message", "Unknown error"))
                failure.set("type", result["status"])

                    # Add stderr if available
                if "stderr" in result and result["stderr"]:
                    failure.text = result["stderr"][:1000]  # Limit size

        # Write XML file
        rough_string = ET.tostring(testsuites, encoding="unicode")
        reparsed = minidom.parseString(rough_string)
        with open(junit_file, "w", encoding="utf-8") as f:
            f.write(reparsed.toprettyxml(indent="  "))

        print(f"JUnit XML report saved to: {junit_file}")

    def _parse_error_details(self, result: Dict[str, Any]) -> None:
        """Parse detailed error information from stderr."""
        stderr = result.get("stderr", "")
        if not stderr:
            return

        import re

        # Parse Python traceback for file and line info
        traceback_match = re.search(r'File "([^"]+)", line (\d+)', stderr)
        if traceback_match:
            result["error_file"] = traceback_match.group(1)
            result["error_line"] = int(traceback_match.group(2))

        # Extract main error message
        lines = stderr.strip().split('\n')
        for line in lines:
            if line.startswith('Traceback') or line.startswith('  File'):
                continue
            if line.strip() and not line.startswith(' ') and ':' in line:
                result["error_details"] = line.strip()
                break

        # Store full traceback
        if 'Traceback' in stderr:
            result["traceback"] = stderr

    def _format_exception(self, exc: Exception) -> str:
        """Format exception with traceback."""
        import traceback
        return ''.join(traceback.format_exception(type(exc), exc, exc.__traceback__))

    def _generate_html_report(self, results: Dict[str, Any]) -> None:
        """Generate HTML test report for visualization."""
        output_dir = Path("output/examples")
        html_file = output_dir / "test_report.html"

        summary = results["summary"]
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(results.get("timestamp", time.time())))

        # Group results by status
        passed = [r for r in results["results"] if r["status"] == "passed"]
        failed = [r for r in results["results"] if r["status"] == "failed"]
        errors = [r for r in results["results"] if r["status"] == "error"]
        timeouts = [r for r in results["results"] if r["status"] == "timeout"]

        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>METAINFORMANT Examples Test Report</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 8px 8px 0 0; }}
        .summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; padding: 30px; }}
        .metric {{ background: #f8f9fa; padding: 20px; border-radius: 8px; text-align: center; border-left: 4px solid #28a745; }}
        .metric.failed {{ border-left-color: #dc3545; }}
        .metric.error {{ border-left-color: #ffc107; }}
        .metric.timeout {{ border-left-color: #17a2b8; }}
        .metric h3 {{ margin: 0 0 10px 0; font-size: 2em; color: #333; }}
        .metric p {{ margin: 0; color: #666; font-size: 0.9em; }}
        .results {{ padding: 30px; }}
        .section {{ margin-bottom: 30px; }}
        .section h2 {{ color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
        .result-item {{ background: #f8f9fa; margin: 10px 0; padding: 15px; border-radius: 6px; border-left: 4px solid #28a745; }}
        .result-item.failed {{ border-left-color: #dc3545; background: #fff5f5; }}
        .result-item.error {{ border-left-color: #ffc107; background: #fffbf0; }}
        .result-item.timeout {{ border-left-color: #17a2b8; background: #f0f8ff; }}
        .result-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; }}
        .result-name {{ font-weight: bold; color: #333; }}
        .result-time {{ color: #666; font-size: 0.9em; }}
        .result-error {{ background: #f8f9fa; padding: 10px; border-radius: 4px; font-family: monospace; font-size: 0.8em; color: #dc3545; margin-top: 10px; white-space: pre-wrap; max-height: 200px; overflow-y: auto; }}
        .tabs {{ display: flex; border-bottom: 1px solid #dee2e6; margin-bottom: 20px; }}
        .tab {{ padding: 10px 20px; cursor: pointer; background: #f8f9fa; border: none; border-bottom: 2px solid transparent; }}
        .tab.active {{ background: white; border-bottom-color: #007bff; color: #007bff; font-weight: bold; }}
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
        .execution-info {{ background: #e9ecef; padding: 15px; border-radius: 6px; margin-bottom: 20px; }}
        .execution-info p {{ margin: 5px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>METAINFORMANT Examples Test Report</h1>
            <p>Generated on {timestamp}</p>
            <div class="execution-info">
                <p><strong>Total Examples:</strong> {summary['total']}</p>
                <p><strong>Execution Time:</strong> {summary['total_execution_time']:.2f}s</p>
                <p><strong>Parallel Execution:</strong> {results.get('parallel_execution', False)}</p>
{f"""                <p><strong>Workers:</strong> {results.get('max_workers', 'N/A')}</p>""" if results.get('parallel_execution') else ""}
            </div>
        </div>

        <div class="summary">
            <div class="metric">
                <h3>{summary['passed']}</h3>
                <p>Passed</p>
            </div>
            <div class="metric failed">
                <h3>{summary['failed']}</h3>
                <p>Failed</p>
            </div>
            <div class="metric error">
                <h3>{summary['error']}</h3>
                <p>Errors</p>
            </div>
            <div class="metric timeout">
                <h3>{summary['timeout']}</h3>
                <p>Timeouts</p>
            </div>
        </div>

        <div class="results">
            <div class="tabs">
                <button class="tab active" onclick="showTab('all')">All Results ({summary['total']})</button>
                {f'<button class="tab" onclick="showTab(\'failed\')">Failed ({len(failed)})</button>' if failed else ''}
                {f'<button class="tab" onclick="showTab(\'errors\')">Errors ({len(errors)})</button>' if errors else ''}
                {f'<button class="tab" onclick="showTab(\'passed\')">Passed ({len(passed)})</button>' if passed else ''}
            </div>

            <div id="all" class="tab-content active">
                <h2>All Results</h2>
                {self._generate_results_html(results['results'])}
            </div>

            {f'<div id="failed" class="tab-content"><h2>Failed Results</h2>{self._generate_results_html(failed)}</div>' if failed else ''}
            {f'<div id="errors" class="tab-content"><h2>Error Results</h2>{self._generate_results_html(errors)}</div>' if errors else ''}
            {f'<div id="passed" class="tab-content"><h2>Passed Results</h2>{self._generate_results_html(passed)}</div>' if passed else ''}
        </div>
    </div>

    <script>
        function showTab(tabName) {{
            // Hide all tab contents
            const contents = document.querySelectorAll('.tab-content');
            contents.forEach(content => content.classList.remove('active'));

            // Remove active class from all tabs
            const tabs = document.querySelectorAll('.tab');
            tabs.forEach(tab => tab.classList.remove('active'));

            // Show selected tab content
            document.getElementById(tabName).classList.add('active');

            // Add active class to clicked tab
            event.target.classList.add('active');
        }}
    </script>
</body>
</html>
"""

        with open(html_file, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f"HTML report saved to: {html_file}")

    def _generate_results_html(self, results: List[Dict[str, Any]]) -> str:
        """Generate HTML for a list of results."""
        if not results:
            return "<p>No results to display.</p>"

        html = ""
        for result in sorted(results, key=lambda x: (x["domain"], x["filename"])):
            status_class = result["status"]
            error_info = ""

            if result.get("error_message"):
                error_info = f'<div class="result-error"><strong>Error:</strong> {result["error_message"]}</div>'

            if result.get("traceback"):
                error_info += f'<div class="result-error"><strong>Traceback:</strong><br>{result["traceback"][:500]}...</div>'

            html += f"""
            <div class="result-item {status_class}">
                <div class="result-header">
                    <span class="result-name">{result['domain']}/{result['filename']}</span>
                    <span class="result-time">{result.get('execution_time', 0):.2f}s</span>
                </div>
                <div>
                    <strong>Status:</strong> {result['status']}
{f"""                    <br><strong>Exit Code:</strong> {result.get('exit_code', 'N/A')}""" if result.get('exit_code') is not None else ""}
{f"""                    <br><strong>Error Line:</strong> {result.get('error_line', 'N/A')}""" if result.get('error_line') else ""}
{f"""                    <br><strong>Files Created:</strong> {len(result.get('output_files_created', []))}""" if result.get('output_files_created') else ""}
                </div>
                {error_info}
            </div>
            """

        return html
        rough_string = ET.tostring(testsuites, encoding="unicode")
        reparsed = minidom.parseString(rough_string)
        with open(junit_file, "w", encoding="utf-8") as f:
            f.write(reparsed.toprettyxml(indent="  "))

        print(f"JUnit XML report saved to: {junit_file}")


def main():
    """Main function."""
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
        tester.save_report(results, junit_xml=args.junit_xml, html=args.html)

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
