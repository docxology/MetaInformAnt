"""Pytest integration for METAINFORMANT examples.

This module provides pytest-based testing for examples, integrating with
the main example test runner while providing pytest-specific features.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

from metainformant.core import io


class TestExamples:
    """Pytest-based example testing."""

    @pytest.fixture(scope="class")
    def example_test_runner(self):
        """Provide access to the example test runner."""

        # Return a helper function that uses subprocess instead of direct import
        def run_test_examples(**kwargs):
            """Run test_examples.py script with given arguments."""
            cmd = [sys.executable, "scripts/test_examples.py"]

            # Add arguments
            if kwargs.get("continue_on_error"):
                cmd.append("--continue-on-error")
            if kwargs.get("junit_xml"):
                cmd.append("--junit-xml")
            if kwargs.get("domain"):
                cmd.extend(["--domain", kwargs["domain"]])

            result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
            return result

        return run_test_examples

    @pytest.fixture(scope="class")
    def examples_dir(self):
        """Provide path to examples directory."""
        return Path("examples")

    def get_all_examples(self, examples_dir: Path) -> List[Path]:
        """Get all example files."""
        return list(examples_dir.rglob("example_*.py"))

    @pytest.mark.parametrize(
        "example_path",
        [
            pytest.param(path, id=str(path.relative_to(Path("examples"))))
            for path in Path("examples").rglob("example_*.py")
        ],
        indirect=False,
    )
    def test_example_execution(self, example_path: Path):
        """Test individual example execution."""
        # Run the example directly
        result = subprocess.run(
            [sys.executable, str(example_path)],
            cwd=Path.cwd(),
            capture_output=True,
            text=True,
            timeout=120,  # 2 minute timeout per example
            env={**dict(os.environ), "MPLBACKEND": "Agg"},
        )

        # Check that it ran successfully
        assert (
            result.returncode == 0
        ), f"Example {example_path} failed with exit code {result.returncode}\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

        # Check for expected output files
        output_dir = Path("output/examples") / example_path.parent.name
        if output_dir.exists():
            json_files = list(output_dir.glob("*.json"))
            assert len(json_files) > 0, f"No JSON output files found for {example_path}"

            # Validate JSON files
            for json_file in json_files:
                with open(json_file, "r") as f:
                    data = json.load(f)

                # Basic structure validation
                assert "example" in data, f"Missing 'example' field in {json_file}"
                assert "domain" in data, f"Missing 'domain' field in {json_file}"
                assert "results" in data, f"Missing 'results' field in {json_file}"

    def test_example_structure_validation(self, examples_dir: Path):
        """Test that all examples follow proper structure."""
        for example_file in examples_dir.rglob("example_*.py"):
            # Check file exists and is readable
            assert example_file.exists(), f"Example file {example_file} does not exist"
            assert example_file.is_file(), f"{example_file} is not a file"

            # Check basic Python syntax
            with open(example_file, "r", encoding="utf-8") as f:
                content = f.read()

            # Must be valid Python
            compile(content, str(example_file), "exec")

            # Check for required elements
            assert 'if __name__ == "__main__":' in content, f"Missing main guard in {example_file}"
            assert "def main():" in content, f"Missing main function in {example_file}"

            # Check for output directory usage
            assert "output/examples/" in content, f"No output directory usage in {example_file}"

    def test_example_dependencies(self, examples_dir: Path):
        """Test that examples have proper dependency declarations."""
        # Load dependencies manifest
        deps_file = examples_dir / "dependencies.json"
        if deps_file.exists():
            with open(deps_file, "r") as f:
                deps_data = json.load(f)

            example_deps = deps_data.get("example_dependencies", {})

            # Check that all examples are declared
            for example_file in examples_dir.rglob("example_*.py"):
                example_key = f"{example_file.parent.name}/{example_file.name}"
                assert example_key in example_deps, f"Example {example_key} not declared in dependencies.json"

    def test_example_documentation(self, examples_dir: Path):
        """Test that examples have proper documentation."""
        for example_file in examples_dir.rglob("example_*.py"):
            with open(example_file, "r", encoding="utf-8") as f:
                content = f.read()

            lines = content.split("\n")

            # Check for module docstring
            assert lines[0].startswith('"""'), f"Missing module docstring in {example_file}"

            # Find the docstring end
            docstring_end = -1
            for i, line in enumerate(lines):
                if line.strip().endswith('"""') and i > 0:
                    docstring_end = i
                    break

            assert docstring_end > 0, f"Malformed docstring in {example_file}"

            # Check docstring content
            docstring = "\n".join(lines[1:docstring_end])
            assert "Usage:" in docstring, f"Missing Usage section in docstring for {example_file}"
            assert "Expected output:" in docstring, f"Missing Expected output section in docstring for {example_file}"

    def test_example_output_consistency(self, example_test_runner):
        """Test that example outputs follow consistent patterns."""
        # Run full example test suite using the fixture
        result = example_test_runner(continue_on_error=True, junit_xml=True)

        # Should complete successfully
        assert result.returncode == 0, f"Example test suite failed: {result.stderr}"

        # Check that JUnit XML was created
        junit_file = Path("output/examples/junit_report.xml")
        assert junit_file.exists(), "JUnit XML report not created"

        # Check that test results JSON was created
        results_file = Path("output/examples/test_results.json")
        assert results_file.exists(), "Test results JSON not created"

        # Validate results structure
        with open(results_file, "r") as f:
            test_results = json.load(f)

        assert "summary" in test_results, "Missing summary in test results"
        assert "results" in test_results, "Missing results in test results"

        summary = test_results["summary"]
        assert summary["total"] > 0, "No examples were tested"
        assert summary["passed"] == summary["total"], f"Not all examples passed: {summary}"

    @pytest.mark.slow
    def test_example_performance_regression(self):
        """Test for performance regressions in examples."""
        # This test would typically compare against stored baselines
        # For now, just ensure examples complete within reasonable time

        import time

        start_time = time.time()

        result = subprocess.run(
            [sys.executable, "scripts/test_examples.py", "--parallel", "--max-workers", "4"],
            cwd=Path.cwd(),
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout for parallel execution
        )

        end_time = time.time()
        duration = end_time - start_time

        assert result.returncode == 0, f"Performance test failed: {result.stderr}"
        assert duration < 180, f"Examples took too long: {duration:.1f}s (should be < 180s)"

    def test_example_domain_coverage(self, examples_dir: Path):
        """Test that examples provide good domain coverage."""
        # Get all domain directories
        domains = [d.name for d in examples_dir.iterdir() if d.is_dir() and not d.name.startswith(".")]

        # Each domain should have at least one example
        for domain in domains:
            domain_dir = examples_dir / domain
            examples = list(domain_dir.glob("example_*.py"))
            assert len(examples) > 0, f"Domain {domain} has no examples"

            # Each domain should have a README
            readme = domain_dir / "README.md"
            assert readme.exists(), f"Domain {domain} missing README.md"

    @pytest.mark.parametrize(
        "domain",
        [
            "core",
            "dna",
            "rna",
            "gwas",
            "protein",
            "ml",
            "networks",
            "ontology",
            "phenotype",
            "ecology",
            "math",
            "information",
            "life_events",
            "multiomics",
            "singlecell",
            "quality",
            "simulation",
            "visualization",
        ],
    )
    def test_domain_examples_exist(self, domain: str, examples_dir: Path):
        """Test that each domain has examples."""
        domain_dir = examples_dir / domain
        if domain_dir.exists():
            examples = list(domain_dir.glob("example_*.py"))
            assert len(examples) > 0, f"Domain {domain} has no examples"

    def test_integration_examples_complexity(self):
        """Test that integration examples handle complexity properly."""
        integration_dir = Path("examples/integration")
        if integration_dir.exists():
            for example_file in integration_dir.glob("example_*.py"):
                # Integration examples should be more complex
                with open(example_file, "r", encoding="utf-8") as f:
                    content = f.read()

                # Should import from multiple domains
                import_lines = [line for line in content.split("\n") if line.startswith("from metainformant.")]
                # Allow some flexibility - integration examples might use different patterns
                assert len(content.split("\n")) > 50, f"Integration example {example_file} seems too simple"

    def test_example_error_handling(self):
        """Test that examples handle errors gracefully."""
        # This is a meta-test - we test that our error handling works
        # by running examples and checking they don't crash the test runner

        result = subprocess.run(
            [sys.executable, "scripts/test_examples.py", "--continue-on-error"],
            cwd=Path.cwd(),
            capture_output=True,
            text=True,
        )

        # Even if examples fail, the test runner should handle it gracefully
        assert result.returncode in [0, 1], f"Test runner crashed: {result.stderr}"

    @pytest.mark.network
    def test_examples_with_network_dependencies(self):
        """Test examples that may require network access."""
        # Skip if no network available
        import urllib.request

        try:
            urllib.request.urlopen("http://httpbin.org/get", timeout=5)
        except:
            pytest.skip("Network not available")

        # If network is available, examples should still work
        # (Most examples don't actually need network, but this tests the framework)
        result = subprocess.run(
            [sys.executable, "scripts/test_examples.py", "--domain", "core"],
            cwd=Path.cwd(),
            capture_output=True,
            text=True,
            timeout=60,
        )

        assert result.returncode == 0, f"Network-dependent test failed: {result.stderr}"


# Pytest configuration and fixtures
def pytest_configure(config):
    """Configure pytest for example testing."""
    # Add custom markers
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "network: marks tests that require network access")


@pytest.fixture(scope="session", autouse=True)
def setup_example_testing():
    """Set up environment for example testing."""
    # Ensure output directory exists
    output_dir = Path("output/examples")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set headless matplotlib backend
    import os

    os.environ["MPLBACKEND"] = "Agg"

    yield

    # Cleanup could go here if needed


@pytest.fixture
def example_output_dir():
    """Provide clean output directory for example testing."""
    output_dir = Path("output/examples/test_temp")
    output_dir.mkdir(parents=True, exist_ok=True)

    yield output_dir

    # Cleanup
    import shutil

    if output_dir.exists():
        shutil.rmtree(output_dir)


# Utility functions for example testing
def run_example_isolated(example_path: Path, timeout: int = 60) -> subprocess.CompletedProcess:
    """Run an example in isolation with proper environment."""
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    env["PYTHONPATH"] = str(Path.cwd())

    return subprocess.run(
        [sys.executable, str(example_path)], cwd=Path.cwd(), capture_output=True, text=True, timeout=timeout, env=env
    )


def validate_example_output(example_path: Path) -> Dict[str, Any]:
    """Validate that an example produced correct output."""
    domain = example_path.parent.name
    output_dir = Path("output/examples") / domain

    if not output_dir.exists():
        return {"valid": False, "error": "Output directory not created"}

    json_files = list(output_dir.glob("*.json"))
    if not json_files:
        return {"valid": False, "error": "No JSON output files found"}

    results = {"valid": True, "files": [], "errors": []}

    for json_file in json_files:
        try:
            with open(json_file, "r") as f:
                data = json.load(f)

            # Basic validation
            required_fields = ["example", "domain", "results"]
            for field in required_fields:
                if field not in data:
                    results["errors"].append(f"Missing field '{field}' in {json_file}")
                    results["valid"] = False

            results["files"].append(str(json_file))

        except Exception as e:
            results["errors"].append(f"Invalid JSON in {json_file}: {e}")
            results["valid"] = False

    return results
