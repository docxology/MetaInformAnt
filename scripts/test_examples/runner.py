#!/usr/bin/env python3
"""Example execution functionality for METAINFORMANT test runner."""

from __future__ import annotations

import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

from .validator import validate_outputs


def run_example(example_path: Path, repo_root: Path, verbose: bool = False) -> Dict[str, Any]:
    """Run a single example and capture results.

    Args:
        example_path: Path to example file
        repo_root: Repository root path
        verbose: Whether to enable verbose output

    Returns:
        Dictionary with execution results
    """
    examples_dir = repo_root / "examples"
    start_time = time.time()

    result = {
        "example_path": str(example_path.relative_to(examples_dir)),
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
        "traceback": None,
    }

    try:
        if verbose:
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
            env=env,
        )

        result["exit_code"] = process.returncode
        result["stdout"] = process.stdout
        result["stderr"] = process.stderr
        result["execution_time"] = time.time() - start_time

        # Parse error details from stderr
        _parse_error_details(result)

        # Validate outputs
        output_files_created, json_outputs_valid, json_outputs_invalid = validate_outputs(example_path)

        result["output_files_created"] = output_files_created
        result["json_outputs_valid"] = json_outputs_valid
        result["json_outputs_invalid"] = json_outputs_invalid

        # Determine status
        if process.returncode == 0:
            # Check if expected outputs were created
            expected_outputs = _get_expected_outputs(example_path)
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
        result["traceback"] = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        result["execution_time"] = time.time() - start_time

    return result


def run_examples_sequential(examples: List[Path], repo_root: Path, verbose: bool = False) -> List[Dict[str, Any]]:
    """Run examples sequentially.

    Args:
        examples: List of example paths
        repo_root: Repository root path
        verbose: Whether to enable verbose output

    Returns:
        List of result dictionaries
    """
    results = []
    for example_path in examples:
        result = run_example(example_path, repo_root, verbose)
        results.append(result)

        # Print status
        status_icon = {"passed": "✓", "failed": "✗", "error": "✗", "timeout": "⏱"}.get(result["status"], "?")
        print(f"{status_icon} {result['example_path']} - {result['status']}")

        if verbose and result["error_message"]:
            print(f"    Error: {result['error_message']}")

    return results


def run_examples_parallel(
    examples: List[Path], repo_root: Path, verbose: bool = False, max_workers: int = 4
) -> List[Dict[str, Any]]:
    """Run examples in parallel using thread pools.

    Args:
        examples: List of example paths
        repo_root: Repository root path
        verbose: Whether to enable verbose output
        max_workers: Maximum number of worker threads

    Returns:
        List of result dictionaries
    """
    # For parallel execution, we need to be careful about output conflicts
    # Each example writes to its own domain subdirectory, so this should be safe

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_example = {
            executor.submit(run_example, example_path, repo_root, verbose): example_path for example_path in examples
        }

        # Process results as they complete
        for future in as_completed(future_to_example):
            example_path = future_to_example[future]
            try:
                result = future.result()
                results.append(result)

                # Print status
                status_icon = {"passed": "✓", "failed": "✗", "error": "✗", "timeout": "⏱"}.get(result["status"], "?")
                print(f"{status_icon} {result['example_path']} - {result['status']}")

                if verbose and result["error_message"]:
                    print(f"    Error: {result['error_message']}")

            except Exception as exc:
                print(f"❌ {example_path.relative_to(repo_root / 'examples')} - error: {exc}")
                # Add error result
                results.append(
                    {
                        "example_path": str(example_path.relative_to(repo_root / "examples")),
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
                        "json_outputs_invalid": [],
                    }
                )

    return results


def _get_expected_outputs(example_path: Path) -> List[str]:
    """Extract expected output files from example docstring.

    Args:
        example_path: Path to example file

    Returns:
        List of expected output patterns
    """
    try:
        with open(example_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Look for "Output:" section in docstring
        lines = content.split("\n")
        in_output_section = False
        outputs = []

        for line in lines:
            if line.strip().startswith("Output:"):
                in_output_section = True
                continue
            elif in_output_section and line.strip().startswith('"""'):
                break
            elif in_output_section and line.strip():
                # Extract output paths
                output_line = line.strip()
                if output_line.startswith("output/"):
                    outputs.append(output_line)

        return outputs

    except Exception:
        return []


def _parse_error_details(result: Dict[str, Any]) -> None:
    """Parse detailed error information from stderr.

    Args:
        result: Result dictionary to update
    """
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
    lines = stderr.strip().split("\n")
    for line in lines:
        if line.startswith("Traceback") or line.startswith("  File"):
            continue
        if line.strip() and not line.startswith(" ") and ":" in line:
            result["error_details"] = line.strip()
            break

    # Store full traceback
    if "Traceback" in stderr:
        result["traceback"] = stderr
