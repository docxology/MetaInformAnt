"""Tests for RNA environment and dependency checking functions.

This module tests all environment checking functions following NO_MOCKING_POLICY.
All tests use real implementations and real CLI tools.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from metainformant.rna import environment


class TestEnvironmentChecks:
    """Test individual environment checking functions."""

    def test_check_amalgkit_returns_tuple(self):
        """Test that check_amalgkit returns a tuple."""
        ok, msg = environment.check_amalgkit()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_check_sra_toolkit_returns_tuple(self):
        """Test that check_sra_toolkit returns a tuple."""
        ok, msg = environment.check_sra_toolkit()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_check_kallisto_returns_tuple(self):
        """Test that check_kallisto returns a tuple."""
        ok, msg = environment.check_kallisto()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_check_metainformant_returns_tuple(self):
        """Test that check_metainformant returns a tuple."""
        ok, msg = environment.check_metainformant()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        # Should always be available since we're running tests
        assert ok is True

    def test_check_virtual_env_returns_tuple(self):
        """Test that check_virtual_env returns a tuple."""
        ok, msg = environment.check_virtual_env()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_check_rscript_returns_tuple(self):
        """Test that check_rscript returns a tuple."""
        ok, msg = environment.check_rscript()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_check_dependencies_returns_dict(self):
        """Test that check_dependencies returns a dictionary."""
        deps = environment.check_dependencies()
        assert isinstance(deps, dict)
        assert len(deps) > 0

        # Verify structure
        expected_keys = {
            "virtual_env",
            "metainformant",
            "amalgkit",
            "sra_toolkit",
            "kallisto",
            "rscript",
        }
        assert set(deps.keys()) == expected_keys

        # Verify each value is a tuple
        for key, value in deps.items():
            assert isinstance(value, tuple)
            assert len(value) == 2
            assert isinstance(value[0], bool)
            assert isinstance(value[1], str)

    def test_validate_environment_returns_dict(self):
        """Test that validate_environment returns a dictionary."""
        result = environment.validate_environment()
        assert isinstance(result, dict)

        # Verify structure
        assert "all_passed" in result
        assert "dependencies" in result
        assert "recommendations" in result

        assert isinstance(result["all_passed"], bool)
        assert isinstance(result["dependencies"], dict)
        assert isinstance(result["recommendations"], list)

        # Verify dependencies match check_dependencies
        deps = environment.check_dependencies()
        assert result["dependencies"] == deps

    def test_validate_environment_all_passed_logic(self):
        """Test that all_passed reflects dependency availability."""
        result = environment.validate_environment()
        deps = result["dependencies"]

        # all_passed should be True only if all dependencies are available
        expected_all_passed = all(available for available, _ in deps.values())
        assert result["all_passed"] == expected_all_passed

    def test_validate_environment_recommendations(self):
        """Test that recommendations are generated for missing dependencies."""
        result = environment.validate_environment()
        deps = result["dependencies"]

        # Count missing dependencies
        missing_count = sum(1 for available, _ in deps.values() if not available)

        # Should have recommendations for missing dependencies
        assert len(result["recommendations"]) >= missing_count

        # Verify recommendations are strings
        for rec in result["recommendations"]:
            assert isinstance(rec, str)
            assert len(rec) > 0

    def test_check_metainformant_always_available(self):
        """Test that metainformant check always succeeds in test environment."""
        ok, msg = environment.check_metainformant()
        assert ok is True
        assert "v" in msg.lower() or "unknown" in msg.lower()

    def test_check_virtual_env_detection(self):
        """Test virtual environment detection logic."""
        ok, msg = environment.check_virtual_env()

        # Check if we're actually in a venv
        in_venv = os.environ.get("VIRTUAL_ENV") is not None or sys.prefix != sys.base_prefix

        # The check should match reality
        assert ok == in_venv
        assert len(msg) > 0


class TestEnvironmentErrorHandling:
    """Test error handling in environment checks."""

    def test_check_amalgkit_handles_missing_tool(self):
        """Test that check_amalgkit handles missing tool gracefully."""
        ok, msg = environment.check_amalgkit()
        # Should return False with error message if not available
        if not ok:
            assert len(msg) > 0
            assert "not found" in msg.lower() or "error" in msg.lower()

    def test_check_sra_toolkit_handles_missing_tool(self):
        """Test that check_sra_toolkit handles missing tool gracefully."""
        ok, msg = environment.check_sra_toolkit()
        if not ok:
            assert len(msg) > 0
            assert "not found" in msg.lower() or "error" in msg.lower()

    def test_check_kallisto_handles_missing_tool(self):
        """Test that check_kallisto handles missing tool gracefully."""
        ok, msg = environment.check_kallisto()
        if not ok:
            assert len(msg) > 0
            assert "not found" in msg.lower() or "error" in msg.lower()

    def test_check_rscript_handles_missing_tool(self):
        """Test that check_rscript handles missing tool gracefully."""
        ok, msg = environment.check_rscript()
        if not ok:
            assert len(msg) > 0
            assert "not found" in msg.lower() or "error" in msg.lower()


class TestEnvironmentDocumentation:
    """Test that environment functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all environment checking functions have docstrings."""
        functions = [
            environment.check_amalgkit,
            environment.check_sra_toolkit,
            environment.check_kallisto,
            environment.check_metainformant,
            environment.check_virtual_env,
            environment.check_rscript,
            environment.check_dependencies,
            environment.validate_environment,
        ]

        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0
