"""Common test utilities and helpers for METAINFORMANT.

This module provides shared utilities for test development, following the
STRICT NO-MOCKING policy. All utilities use real implementations and
graceful degradation for missing dependencies.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pytest


# External tool checking utilities
def check_external_tool_available(tool_name: str) -> bool:
    """Check if an external CLI tool is available on PATH."""
    return shutil.which(tool_name) is not None


def require_external_tool(tool_name: str, reason: str = "") -> None:
    """Require an external tool for a test, skip if not available."""
    if not check_external_tool_available(tool_name):
        reason_msg = f" - {reason}" if reason else ""
        pytest.skip(f"External tool '{tool_name}' not available on PATH{reason_msg}")


def get_available_external_tools(tool_names: List[str]) -> List[str]:
    """Get list of available external tools from a list of tool names."""
    return [tool for tool in tool_names if check_external_tool_available(tool)]


# Network connectivity utilities
def check_network_connectivity(url: str = "https://httpbin.org/status/200", timeout: int = 5) -> bool:
    """Check if network connectivity is available."""
    try:
        import requests

        response = requests.get(url, timeout=timeout)
        return response.status_code == 200
    except Exception:
        return False


def require_network_connectivity(reason: str = "") -> None:
    """Require network connectivity for a test, skip if not available."""
    if not check_network_connectivity():
        reason_msg = f" - {reason}" if reason else ""
        pytest.skip(f"Network connectivity not available{reason_msg}")


# Test data generation helpers
def generate_sample_dna_sequence(length: int = 100, gc_content: float = 0.5) -> str:
    """Generate a sample DNA sequence with specified GC content."""
    import random

    bases = ["A", "T", "G", "C"]
    # Adjust base probabilities to achieve target GC content
    gc_prob = gc_content
    at_prob = (1 - gc_content) / 2

    weights = [at_prob, at_prob, gc_prob / 2, gc_prob / 2]  # A, T, G, C

    sequence = "".join(random.choices(bases, weights=weights, k=length))
    return sequence


def generate_sample_protein_sequence(length: int = 50) -> str:
    """Generate a sample protein sequence."""
    import random

    # Common amino acids with realistic frequencies
    amino_acids = "ARNDCQEGHILKMFPSTWYV"
    weights = [1] * len(amino_acids)  # Equal weights for simplicity

    sequence = "".join(random.choices(list(amino_acids), weights=weights, k=length))
    return sequence


def create_sample_fasta_file(output_path: Path, sequences: Dict[str, str]) -> Path:
    """Create a sample FASTA file with the given sequences."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")

    return output_path


def create_sample_fastq_file(output_path: Path, reads: Dict[str, str], quality: str = "I" * 50) -> Path:
    """Create a sample FASTQ file with the given reads."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for read_id, sequence in reads.items():
            f.write(f"@{read_id}\n{sequence}\n+\n{quality[:len(sequence)]}\n")

    return output_path


# Environment variable management
class EnvironmentVarManager:
    """Context manager for temporary environment variable management."""

    def __init__(self, **env_vars):
        self.env_vars = env_vars
        self.original_values = {}

    def __enter__(self):
        # Store original values
        for var_name in self.env_vars:
            self.original_values[var_name] = os.environ.get(var_name)

        # Set new values
        for var_name, var_value in self.env_vars.items():
            if var_value is None:
                os.environ.pop(var_name, None)
            else:
                os.environ[var_name] = str(var_value)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore original values
        for var_name, original_value in self.original_values.items():
            if original_value is None:
                os.environ.pop(var_name, None)
            else:
                os.environ[var_name] = original_value


def with_env_vars(**env_vars):
    """Decorator to run a test function with specific environment variables."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            with EnvironmentVarManager(**env_vars):
                return func(*args, **kwargs)

        return wrapper

    return decorator


# Output directory management
def ensure_output_directory(subdir: str = "") -> Path:
    """Ensure output directory exists and return its path."""
    from metainformant.core import io

    output_path = Path("output")
    if subdir:
        output_path = output_path / subdir

    return io.ensure_directory(output_path)


def create_isolated_output_dir(prefix: str = "test_") -> Path:
    """Create an isolated output directory for a test."""
    return ensure_output_directory(f"{prefix}{os.getpid()}")


# Dependency verification helpers
def check_python_module_available(module_name: str) -> bool:
    """Check if a Python module is available for import."""
    try:
        __import__(module_name)
        return True
    except ImportError:
        return False


def require_python_module(module_name: str, reason: str = "") -> None:
    """Require a Python module for a test, skip if not available."""
    if not check_python_module_available(module_name):
        reason_msg = f" - {reason}" if reason else ""
        pytest.skip(f"Python module '{module_name}' not available{reason_msg}")


def check_uv_dependency_group(group_name: str) -> bool:
    """Check if a uv dependency group is available."""
    # This is a simplified check - in practice, we'd check if uv sync was run for this group
    # For now, just check if uv is available
    return check_external_tool_available("uv")


# File and directory utilities
def create_temporary_file(content: str = "", suffix: str = ".txt") -> Path:
    """Create a temporary file with optional content."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    try:
        if content:
            with os.fdopen(fd, "w") as f:
                f.write(content)
        else:
            os.close(fd)
        return Path(path)
    except Exception:
        os.close(fd)
        raise


def create_temporary_directory() -> Path:
    """Create a temporary directory."""
    return Path(tempfile.mkdtemp())


# Test timing and performance utilities
def measure_execution_time(func, *args, **kwargs) -> tuple[Any, float]:
    """Measure execution time of a function."""
    import time

    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()

    return result, end_time - start_time


# Data validation helpers
def validate_fasta_file(file_path: Path) -> bool:
    """Validate that a file is properly formatted FASTA."""
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        if not lines or not lines[0].startswith(">"):
            return False

        # Check that sequence lines contain valid nucleotides
        valid_bases = set("ATCGNUatcgnu-")
        for line in lines[1:]:
            line = line.strip()
            if line.startswith(">"):
                continue
            if not all(c in valid_bases for c in line):
                return False

        return True
    except Exception:
        return False


def validate_fastq_file(file_path: Path) -> bool:
    """Validate that a file is properly formatted FASTQ."""
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        if len(lines) % 4 != 0:
            return False

        for i in range(0, len(lines), 4):
            # Check FASTQ format: @seq, sequence, +, quality
            if not lines[i].startswith("@"):
                return False
            if not lines[i + 2].startswith("+"):
                return False
            # Quality string should be same length as sequence
            if len(lines[i + 1].strip()) != len(lines[i + 3].strip()):
                return False

        return True
    except Exception:
        return False


# Test result analysis helpers
def count_test_functions_in_file(file_path: Path) -> int:
    """Count the number of test functions in a Python file."""
    try:
        with open(file_path, "r") as f:
            content = f.read()

        # Simple regex-based counting (could be improved)
        import re

        test_functions = re.findall(r"def test_\w+", content)
        return len(test_functions)
    except Exception:
        return 0


def analyze_test_coverage(test_results: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze test coverage results."""
    # This would integrate with coverage.py results
    # For now, return a basic structure
    return {
        "total_statements": 0,
        "covered_statements": 0,
        "coverage_percentage": 0.0,
        "missing_lines": [],
    }


# Platform detection utilities
def get_platform_info() -> Dict[str, str]:
    """Get information about the current platform."""
    import platform

    return {
        "system": platform.system(),
        "release": platform.release(),
        "version": platform.version(),
        "machine": platform.machine(),
        "processor": platform.processor(),
    }


def is_ci_environment() -> bool:
    """Check if running in a CI environment."""
    ci_indicators = [
        "CI",
        "CONTINUOUS_INTEGRATION",
        "TRAVIS",
        "APPVEYOR",
        "CIRCLECI",
        "GITHUB_ACTIONS",
        "GITLAB_CI",
        "JENKINS",
    ]

    return any(os.environ.get(indicator) for indicator in ci_indicators)
