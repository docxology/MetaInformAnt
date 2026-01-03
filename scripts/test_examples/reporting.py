#!/usr/bin/env python3
"""Reporting functionality for METAINFORMANT test runner."""

from __future__ import annotations

import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List
from xml.dom import minidom


def save_json_report(results: Dict[str, Any], output_dir: Path) -> None:
    """Save test results to JSON file.

    Args:
        results: Results dictionary
        output_dir: Output directory
    """
    report_file = output_dir / "test_results.json"
    from metainformant.core import io
    io.dump_json(results, report_file, indent=2)
    print(f"JSON report saved to: {report_file}")


def generate_junit_xml(results: Dict[str, Any], output_dir: Path) -> None:
    """Generate JUnit XML report for CI/CD systems.

    Args:
        results: Results dictionary
        output_dir: Output directory
    """
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


def generate_html_report(results: Dict[str, Any], output_dir: Path) -> None:
    """Generate HTML test report for visualization.

    Args:
        results: Results dictionary
        output_dir: Output directory
    """
    html_file = output_dir / "test_report.html"

    summary = results["summary"]
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(results.get("timestamp", time.time())))

    # Group results by status
    passed = [r for r in results["results"] if r["status"] == "passed"]
    failed = [r for r in results["results"] if r["status"] == "failed"]
    errors = [r for r in results["results"] if r["status"] == "error"]
    timeouts = [r for r in results["results"] if r["status"] == "timeout"]

    # Build HTML content
    html_parts = [
        '<!DOCTYPE html>',
        '<html lang="en">',
        '<head>',
        '    <meta charset="UTF-8">',
        '    <meta name="viewport" content="width=device-width, initial-scale=1.0">',
        '    <title>METAINFORMANT Examples Test Report</title>',
        '    <style>',
        '        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }',
        '        .container { max-width: 1200px; margin: 0 auto; background: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }',
        '        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 8px 8px 0 0; }',
        '        .summary { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; padding: 30px; }',
        '        .metric { background: #f8f9fa; padding: 20px; border-radius: 8px; text-align: center; border-left: 4px solid #28a745; }',
        '        .metric.failed { border-left-color: #dc3545; }',
        '        .metric.error { border-left-color: #ffc107; }',
        '        .metric.timeout { border-left-color: #17a2b8; }',
        '        .metric h3 { margin: 0 0 10px 0; font-size: 2em; color: #333; }',
        '        .metric p { margin: 0; color: #666; font-size: 0.9em; }',
        '        .results { padding: 30px; }',
        '        .section h2 { color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }',
        '        .result-item { background: #f8f9fa; margin: 10px 0; padding: 15px; border-radius: 6px; border-left: 4px solid #28a745; }',
        '        .result-item.failed { border-left-color: #dc3545; background: #fff5f5; }',
        '        .result-item.error { border-left-color: #ffc107; background: #fffbf0; }',
        '        .result-item.timeout { border-left-color: #17a2b8; background: #f0f8ff; }',
        '        .result-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; }',
        '        .result-name { font-weight: bold; color: #333; }',
        '        .result-time { color: #666; font-size: 0.9em; }',
        '        .result-error { background: #f8f9fa; padding: 10px; border-radius: 4px; font-family: monospace; font-size: 0.8em; color: #dc3545; margin-top: 10px; white-space: pre-wrap; max-height: 200px; overflow-y: auto; }',
        '        .tabs { display: flex; border-bottom: 1px solid #dee2e6; margin-bottom: 20px; }',
        '        .tab { padding: 10px 20px; cursor: pointer; background: #f8f9fa; border: none; border-bottom: 2px solid transparent; }',
        '        .tab.active { background: white; border-bottom-color: #007bff; color: #007bff; font-weight: bold; }',
        '        .tab-content { display: none; }',
        '        .tab-content.active { display: block; }',
        '        .execution-info { background: #e9ecef; padding: 15px; border-radius: 6px; margin-bottom: 20px; }',
        '        .execution-info p { margin: 5px 0; }',
        '    </style>',
        '</head>',
        '<body>',
        '    <div class="container">',
        '        <div class="header">',
        '            <h1>METAINFORMANT Examples Test Report</h1>',
        f'            <p>Generated on {timestamp}</p>',
        '            <div class="execution-info">',
        f'                <p><strong>Total Examples:</strong> {summary["total"]}</p>',
        f'                <p><strong>Execution Time:</strong> {summary["total_execution_time"]:.2f}s</p>',
        f'                <p><strong>Parallel Execution:</strong> {results.get("parallel_execution", False)}</p>',
    ]

    if results.get('parallel_execution'):
        html_parts.append(f'                <p><strong>Workers:</strong> {results.get("max_workers", "N/A")}</p>')

    html_parts.extend([
        '            </div>',
        '        </div>',
        '',
        '        <div class="summary">',
        '            <div class="metric">',
        f'                <h3>{summary["passed"]}</h3>',
        '                <p>Passed</p>',
        '            </div>',
        '            <div class="metric failed">',
        f'                <h3>{summary["failed"]}</h3>',
        '                <p>Failed</p>',
        '            </div>',
        '            <div class="metric error">',
        f'                <h3>{summary["error"]}</h3>',
        '                <p>Errors</p>',
        '            </div>',
        '            <div class="metric timeout">',
        f'                <h3>{summary["timeout"]}</h3>',
        '                <p>Timeouts</p>',
        '            </div>',
        '        </div>',
        '',
        '        <div class="results">',
        '            <div class="tabs">',
        f'                <button class="tab active" onclick="showTab(\'all\')">All Results ({summary["total"]})</button>',
    ])

    if failed:
        html_parts.append(f'                <button class="tab" onclick="showTab(\'failed\')">Failed ({len(failed)})</button>')
    if errors:
        html_parts.append(f'                <button class="tab" onclick="showTab(\'errors\')">Errors ({len(errors)})</button>')
    if passed:
        html_parts.append(f'                <button class="tab" onclick="showTab(\'passed\')">Passed ({len(passed)})</button>')

    html_parts.extend([
        '            </div>',
        '',
        '            <div id="all" class="tab-content active">',
        '                <h2>All Results</h2>',
        _generate_results_html(results['results']),
        '            </div>',
    ])

    if failed:
        html_parts.extend([
            '            <div id="failed" class="tab-content">',
            '                <h2>Failed Results</h2>',
            _generate_results_html(failed),
            '            </div>',
        ])

    if errors:
        html_parts.extend([
            '            <div id="errors" class="tab-content">',
            '                <h2>Error Results</h2>',
            _generate_results_html(errors),
            '            </div>',
        ])

    if passed:
        html_parts.extend([
            '            <div id="passed" class="tab-content">',
            '                <h2>Passed Results</h2>',
            _generate_results_html(passed),
            '            </div>',
        ])

    html_parts.extend([
        '        </div>',
        '    </div>',
        '',
        '    <script>',
        '        function showTab(tabName) {',
        '            // Hide all tab contents',
        '            const contents = document.querySelectorAll(\'.tab-content\');',
        '            contents.forEach(content => content.classList.remove(\'active\'));',
        '',
        '            // Remove active class from all tabs',
        '            const tabs = document.querySelectorAll(\'.tab\');',
        '            tabs.forEach(tab => tab.classList.remove(\'active\'));',
        '',
        '            // Show selected tab content',
        '            document.getElementById(tabName).classList.add(\'active\');',
        '',
        '            // Add active class to clicked tab',
        '            event.target.classList.add(\'active\');',
        '        }',
        '    </script>',
        '</body>',
        '</html>',
    ])

    html_content = '\n'.join(html_parts)

    with open(html_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"HTML report saved to: {html_file}")


def _generate_results_html(results: List[Dict[str, Any]]) -> str:
    """Generate HTML for a list of results.

    Args:
        results: List of result dictionaries

    Returns:
        HTML string
    """
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

        html_parts = [
            f'<div class="result-item {status_class}">',
            '<div class="result-header">',
            f'<span class="result-name">{result["domain"]}/{result["filename"]}</span>',
            f'<span class="result-time">{result.get("execution_time", 0):.2f}s</span>',
            '</div>',
            '<div>',
            f'<strong>Status:</strong> {result["status"]}'
        ]

        if result.get('exit_code') is not None:
            html_parts.append(f'<br><strong>Exit Code:</strong> {result.get("exit_code", "N/A")}')
        if result.get('error_line'):
            html_parts.append(f'<br><strong>Error Line:</strong> {result.get("error_line", "N/A")}')
        if result.get('output_files_created'):
            html_parts.append(f'<br><strong>Files Created:</strong> {len(result.get("output_files_created", []))}')

        html_parts.extend(['</div>', error_info, '</div>'])
        html += '\n'.join(html_parts)

    return html


def print_summary(results: Dict[str, Any]) -> None:
    """Print test summary to console.

    Args:
        results: Results dictionary
    """
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
