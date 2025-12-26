"""Pytest configuration and fixtures for METAINFORMANT examples testing."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Generator

import pytest

import metainformant


@pytest.fixture(scope="session")
def metainformant_version() -> str:
    """Provide METAINFORMANT version for testing."""
    return metainformant.__version__


@pytest.fixture(scope="session")
def examples_base_dir() -> Path:
    """Provide the base examples directory."""
    return Path("examples")


@pytest.fixture(scope="session")
def example_domains(examples_base_dir: Path) -> list[str]:
    """Provide list of available example domains."""
    return [d.name for d in examples_base_dir.iterdir()
            if d.is_dir() and not d.name.startswith('.')]


@pytest.fixture(params=[
    "core", "dna", "rna", "gwas", "protein", "epigenome", "ontology",
    "phenotype", "ecology", "math", "information", "life_events",
    "multiomics", "singlecell", "quality", "networks", "ml",
    "simulation", "visualization"
])
def domain(request, example_domains) -> str:
    """Parametrize tests across all domains."""
    domain = request.param
    if domain not in example_domains:
        pytest.skip(f"Domain {domain} not available")
    return domain


@pytest.fixture
def domain_examples_dir(examples_base_dir: Path, domain: str) -> Path:
    """Provide the directory for a specific domain's examples."""
    domain_dir = examples_base_dir / domain
    if not domain_dir.exists():
        pytest.skip(f"Domain directory {domain_dir} does not exist")
    return domain_dir


@pytest.fixture
def domain_examples(domain_examples_dir: Path) -> list[Path]:
    """Provide list of example files for a domain."""
    return list(domain_examples_dir.glob("example_*.py"))


@pytest.fixture
def example_file(request, domain_examples) -> Path:
    """Provide a specific example file for testing."""
    if not domain_examples:
        pytest.skip("No examples available in domain")

    # Allow parametrization by example name
    example_name = getattr(request, 'param', None)
    if example_name:
        example_path = None
        for ex in domain_examples:
            if ex.name == f"example_{example_name}.py":
                example_path = ex
                break
        if not example_path:
            pytest.skip(f"Example {example_name} not found")
        return example_path
    else:
        # Return first example as default
        return domain_examples[0]


@pytest.fixture
def example_output_dir(example_file: Path) -> Generator[Path, None, None]:
    """Provide a clean output directory for example testing."""
    domain = example_file.parent.name
    output_dir = Path("output/examples") / domain

    # Create if doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Clean any existing files for this example
    example_name = example_file.stem
    for file_path in output_dir.glob(f"{example_name}*"):
        if file_path.is_file():
            file_path.unlink()

    yield output_dir


@pytest.fixture
def temp_example_dir() -> Generator[Path, None, None]:
    """Provide a temporary directory for example testing."""
    with tempfile.TemporaryDirectory(prefix="metainformant_examples_") as temp_dir:
        yield Path(temp_dir)


@pytest.fixture
def mock_example_environment(monkeypatch, temp_example_dir: Path) -> Dict[str, Any]:
    """Set up a mock environment for example testing."""
    # Mock output directory
    mock_output = temp_example_dir / "output" / "examples"
    mock_output.mkdir(parents=True)

    # Set environment variables
    monkeypatch.setenv("MPLBACKEND", "Agg")
    monkeypatch.setattr("metainformant.core.paths.PLATFORM", "test")

    return {
        "output_dir": mock_output,
        "temp_dir": temp_example_dir,
        "env_vars": dict(os.environ)
    }


@pytest.fixture(scope="session")
def example_dependencies() -> Dict[str, Any]:
    """Load example dependencies manifest."""
    deps_file = Path("examples/dependencies.json")
    if deps_file.exists():
        import json
        with open(deps_file, 'r') as f:
            return json.load(f)
    return {}


@pytest.fixture
def example_dependencies_for_domain(example_dependencies: Dict[str, Any], domain: str) -> Dict[str, Any]:
    """Provide dependencies for a specific domain."""
    return example_dependencies.get("domain_dependencies", {}).get(domain, {})


@pytest.fixture
def check_example_dependencies(example_dependencies: Dict[str, Any]):
    """Provide a function to check if example dependencies are available."""

    def _check_deps(example_path: str) -> Dict[str, Any]:
        """Check dependencies for a specific example."""
        if example_path in example_dependencies.get("example_dependencies", {}):
            deps_config = example_dependencies["example_dependencies"][example_path]

            # Check required dependencies
            required = deps_config.get("required", [])
            optional = deps_config.get("optional", [])

            missing_required = []
            missing_optional = []

            for dep in required:
                if not _is_package_available(dep):
                    missing_required.append(dep)

            for dep in optional:
                if not _is_package_available(dep):
                    missing_optional.append(dep)

            return {
                "available": len(missing_required) == 0,
                "missing_required": missing_required,
                "missing_optional": missing_optional
            }
        else:
            return {"available": True, "missing_required": [], "missing_optional": []}

    def _is_package_available(package_name: str) -> bool:
        """Check if a Python package is available."""
        try:
            if "." in package_name:
                # Handle submodules
                module_parts = package_name.split(".")
                import importlib
                module = importlib.import_module(module_parts[0])
                for part in module_parts[1:]:
                    module = getattr(module, part)
                return True
            else:
                import importlib
                importlib.import_module(package_name)
                return True
        except ImportError:
            return False

    return _check_deps


@pytest.fixture
def skip_if_missing_dependencies(check_example_dependencies):
    """Skip test if required dependencies are missing."""

    def _skip_if_missing(example_path: str):
        deps_status = check_example_dependencies(example_path)
        if not deps_status["available"]:
            missing = deps_status["missing_required"]
            pytest.skip(f"Missing required dependencies: {', '.join(missing)}")

    return _skip_if_missing


@pytest.fixture
def performance_tracker() -> Dict[str, Any]:
    """Provide a performance tracker for examples."""
    import time

    tracker = {
        "start_times": {},
        "end_times": {},
        "durations": {}
    }

    def start_tracking(name: str):
        tracker["start_times"][name] = time.time()

    def stop_tracking(name: str):
        if name in tracker["start_times"]:
            tracker["end_times"][name] = time.time()
            tracker["durations"][name] = tracker["end_times"][name] - tracker["start_times"][name]

    def get_duration(name: str) -> float:
        return tracker["durations"].get(name, 0.0)

    tracker["start"] = start_tracking
    tracker["stop"] = stop_tracking
    tracker["get_duration"] = get_duration

    return tracker


@pytest.fixture
def example_validator():
    """Provide an example validator for testing."""

    def validate_example_structure(example_path: Path) -> Dict[str, Any]:
        """Validate basic example structure."""
        issues = []

        if not example_path.exists():
            return {"valid": False, "issues": ["File does not exist"]}

        try:
            with open(example_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Check syntax
            compile(content, str(example_path), 'exec')

            # Check required elements
            if 'if __name__ == "__main__":' not in content:
                issues.append("Missing main guard")

            if 'def main():' not in content:
                issues.append("Missing main function")

            if 'output/examples/' not in content:
                issues.append("No output directory usage")

            return {
                "valid": len(issues) == 0,
                "issues": issues
            }

        except SyntaxError as e:
            return {"valid": False, "issues": [f"Syntax error: {e}"]}
        except Exception as e:
            return {"valid": False, "issues": [f"Validation error: {e}"]}

    return validate_example_structure


# Markers for conditional testing
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers",
        "requires_amalgkit: mark test as requiring amalgkit CLI tool"
    )
    config.addinivalue_line(
        "markers",
        "requires_network: mark test as requiring network access"
    )
    config.addinivalue_line(
        "markers",
        "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers",
        "example_domain(domain): mark test as specific to an example domain"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection based on available dependencies."""
    # Check for amalgkit
    amalgkit_available = _check_external_tool("amalgkit")

    for item in items:
        # Skip amalgkit tests if not available
        if "requires_amalgkit" in item.keywords and not amalgkit_available:
            item.add_marker(pytest.mark.skip(reason="amalgkit not available"))

        # Skip network tests if network not available
        if "requires_network" in item.keywords:
            if not _check_network_available():
                item.add_marker(pytest.mark.skip(reason="network not available"))


def _check_external_tool(tool_name: str) -> bool:
    """Check if an external tool is available."""
    import subprocess

    try:
        result = subprocess.run(
            ["which", tool_name],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except:
        return False


def _check_network_available() -> bool:
    """Check if network is available."""
    import urllib.request

    try:
        urllib.request.urlopen('http://httpbin.org/get', timeout=5)
        return True
    except:
        return False
