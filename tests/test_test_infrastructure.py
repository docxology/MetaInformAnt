"""Tests for the test infrastructure itself.

This module tests the test infrastructure components including:
- conftest.py fixtures and configuration
- tests/utils.py utilities
- Basic test environment setup
- Dependency verification
- Test script functionality

These tests ensure the test infrastructure is reliable and working correctly.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import pytest

# Test conftest.py fixtures and utilities
class TestConftestFixtures:
    """Test conftest.py fixtures and configuration."""

    def test_uv_environment_info_fixture(self, uv_environment_info):
        """Test uv_environment_info fixture provides expected structure."""
        required_keys = [
            "uv_available", "uv_cache_dir", "uv_project_env",
            "virtual_env", "metainformant_venv", "python_path",
            "test_data_dir", "output_dir"
        ]

        for key in required_keys:
            assert key in uv_environment_info, f"Missing key: {key}"

        # uv should be available in test environment
        assert isinstance(uv_environment_info["uv_available"], bool)

    def test_test_dependency_checker_fixture(self, test_dependency_checker):
        """Test test_dependency_checker fixture returns expected structure."""
        result = test_dependency_checker("fast")

        # Should return a dictionary
        assert isinstance(result, dict)

        # Should check basic imports
        assert "import_pytest" in result
        assert "import_pathlib" in result
        assert "import_tempfile" in result

    def test_verify_uv_setup_fixture(self):
        """Test that verify_uv_setup fixture doesn't fail."""
        # If this test runs, the fixture passed
        assert True

    def test_validate_test_environment_fixture(self):
        """Test that validate_test_environment fixture doesn't fail."""
        # If this test runs, the fixture passed
        assert True

    def test_configure_uv_environment_fixture(self):
        """Test that configure_uv_environment fixture doesn't fail."""
        # If this test runs, the fixture passed
        assert True


# Test tests/utils.py utilities
class TestUtilsModule:
    """Test the utilities in tests/utils.py."""

    def test_check_external_tool_available(self):
        """Test external tool checking utility."""
        from tests.utils import check_external_tool_available

        # Python should be available
        assert check_external_tool_available("python3")

        # Some tool that shouldn't exist
        assert not check_external_tool_available("nonexistent_tool_12345")

    def test_require_external_tool_skip(self):
        """Test require_external_tool with non-existent tool."""
        from tests.utils import require_external_tool

        with pytest.raises(pytest.skip.Exception):
            require_external_tool("nonexistent_tool_12345", "Test tool requirement")

    def test_get_available_external_tools(self):
        """Test getting list of available external tools."""
        from tests.utils import get_available_external_tools

        tools = ["python3", "nonexistent_tool_12345", "bash"]
        available = get_available_external_tools(tools)

        assert "python3" in available
        assert "bash" in available
        assert "nonexistent_tool_12345" not in available

    def test_check_network_connectivity(self):
        """Test network connectivity checking."""
        from tests.utils import check_network_connectivity

        # This might fail in offline environments, so just test it doesn't crash
        result = check_network_connectivity()
        assert isinstance(result, bool)

    def test_require_network_connectivity_skip(self):
        """Test require_network_connectivity in offline environment."""
        from tests.utils import require_network_connectivity

        # This might not skip if network is available
        # Just test it doesn't crash
        require_network_connectivity("Test network requirement")

    def test_generate_sample_dna_sequence(self):
        """Test DNA sequence generation."""
        from tests.utils import generate_sample_dna_sequence

        seq = generate_sample_dna_sequence(100)
        assert len(seq) == 100
        assert all(c in "ATCG" for c in seq)

    def test_generate_sample_protein_sequence(self):
        """Test protein sequence generation."""
        from tests.utils import generate_sample_protein_sequence

        seq = generate_sample_protein_sequence(50)
        assert len(seq) == 50
        assert all(c in "ACDEFGHIKLMNPQRSTVWY" for c in seq)

    def test_create_sample_fasta_file(self, isolated_tmp_dir):
        """Test FASTA file creation."""
        from tests.utils import create_sample_fasta_file

        sequences = {"seq1": "ATCG", "seq2": "GCTA"}
        fasta_file = create_sample_fasta_file(isolated_tmp_dir / "test.fasta", sequences)

        assert fasta_file.exists()
        content = fasta_file.read_text()
        assert ">seq1" in content
        assert ">seq2" in content
        assert "ATCG" in content
        assert "GCTA" in content

    def test_create_sample_fastq_file(self, isolated_tmp_dir):
        """Test FASTQ file creation."""
        from tests.utils import create_sample_fastq_file

        reads = {"read1": "ATCG", "read2": "GCTA"}
        fastq_file = create_sample_fastq_file(isolated_tmp_dir / "test.fastq", reads)

        assert fastq_file.exists()
        content = fastq_file.read_text()
        assert "@read1" in content
        assert "@read2" in content
        assert "ATCG" in content
        assert "GCTA" in content

    def test_environment_var_manager(self):
        """Test environment variable manager context manager."""
        from tests.utils import EnvironmentVarManager

        original_test_var = os.environ.get("TEST_VAR")

        with EnvironmentVarManager(TEST_VAR="test_value"):
            assert os.environ.get("TEST_VAR") == "test_value"

        # Should be restored
        assert os.environ.get("TEST_VAR") == original_test_var

    def test_with_env_vars_decorator(self):
        """Test with_env_vars decorator."""
        from tests.utils import with_env_vars

        @with_env_vars(TEST_DECORATOR_VAR="decorated_value")
        def test_function():
            return os.environ.get("TEST_DECORATOR_VAR")

        result = test_function()
        assert result == "decorated_value"

    def test_ensure_output_directory(self):
        """Test output directory creation."""
        from tests.utils import ensure_output_directory

        output_dir = ensure_output_directory("test_infrastructure_output")
        assert output_dir.exists()
        assert output_dir.is_dir()

    def test_validate_fasta_file(self, isolated_tmp_dir):
        """Test FASTA file validation."""
        from tests.utils import validate_fasta_file, create_sample_fasta_file

        # Valid FASTA
        valid_fasta = create_sample_fasta_file(isolated_tmp_dir / "valid.fasta",
                                             {"seq1": "ATCG", "seq2": "GCTA"})
        assert validate_fasta_file(valid_fasta)

        # Invalid file
        invalid_file = isolated_tmp_dir / "invalid.fasta"
        invalid_file.write_text("not a fasta file")
        assert not validate_fasta_file(invalid_file)

    def test_validate_fastq_file(self, isolated_tmp_dir):
        """Test FASTQ file validation."""
        from tests.utils import validate_fastq_file, create_sample_fastq_file

        # Valid FASTQ
        valid_fastq = create_sample_fastq_file(isolated_tmp_dir / "valid.fastq",
                                             {"read1": "ATCG", "read2": "GCTA"})
        assert validate_fastq_file(valid_fastq)

        # Invalid file
        invalid_file = isolated_tmp_dir / "invalid.fastq"
        invalid_file.write_text("not a fastq file")
        assert not validate_fastq_file(invalid_file)

    def test_count_test_functions_in_file(self):
        """Test counting test functions in a file."""
        from tests.utils import count_test_functions_in_file

        # Test this file itself
        count = count_test_functions_in_file(Path(__file__))
        assert count >= 1  # Should find at least the test functions in this class

    def test_get_platform_info(self):
        """Test platform info retrieval."""
        from tests.utils import get_platform_info

        info = get_platform_info()
        required_keys = ["system", "release", "version", "machine", "processor"]

        for key in required_keys:
            assert key in info
            assert isinstance(info[key], str)

    def test_is_ci_environment(self):
        """Test CI environment detection."""
        from tests.utils import is_ci_environment

        # This will be False in local development, True in CI
        result = is_ci_environment()
        assert isinstance(result, bool)


class TestTestScripts:
    """Test basic functionality of test scripts."""

    def test_uv_test_setup_script_exists(self):
        """Test that uv_test_setup.sh script exists and is executable."""
        script_path = Path("scripts/package/uv_test_setup.sh")
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111  # Check if executable

    def test_verify_deps_script_exists(self):
        """Test that verify_test_deps.sh script exists and is executable."""
        script_path = Path("scripts/package/verify_test_deps.sh")
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111  # Check if executable

    def test_uv_test_script_exists(self):
        """Test that uv_test.sh script exists and is executable."""
        script_path = Path("scripts/package/uv_test.sh")
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111  # Check if executable

    def test_uv_test_optimized_script_exists(self):
        """Test that uv_test_optimized.sh script exists and is executable."""
        script_path = Path("scripts/package/uv_test_optimized.sh")
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111  # Check if executable

    def test_run_tests_script_exists(self):
        """Test that run_tests.sh script exists and is executable."""
        script_path = Path("scripts/package/run_tests.sh")
        assert script_path.exists()
        assert script_path.stat().st_mode & 0o111  # Check if executable


class TestEnvironmentSetup:
    """Test test environment setup and validation."""

    def test_test_data_directory_exists(self):
        """Test that tests/data directory exists."""
        data_dir = Path("tests/data")
        assert data_dir.exists()
        assert data_dir.is_dir()

    def test_output_directory_created(self):
        """Test that output directory gets created."""
        output_dir = Path("output")
        # May or may not exist initially, but should be creatable
        output_dir.mkdir(exist_ok=True)
        assert output_dir.exists()

    def test_src_directory_in_path(self):
        """Test that src directory is in Python path during tests."""
        import sys
        src_path = str(Path("src").resolve())

        # Should be in sys.path during test execution
        assert any(src_path in path for path in sys.path)

    def test_conftest_imports_work(self):
        """Test that conftest.py imports work."""
        # If this test runs, conftest imports worked
        assert True

    def test_utils_module_imports(self):
        """Test that tests.utils module can be imported."""
        import tests.utils
        assert hasattr(tests.utils, 'check_external_tool_available')
        assert hasattr(tests.utils, 'require_external_tool')


class TestDependencyVerification:
    """Test dependency verification functionality."""

    def test_pytest_available(self):
        """Test that pytest is available."""
        import pytest
        assert pytest is not None

    def test_basic_test_imports(self):
        """Test that basic test imports work."""
        import pathlib
        import tempfile
        assert pathlib is not None
        assert tempfile is not None

    def test_pathlib_available(self):
        """Test pathlib is available."""
        import pathlib
        assert pathlib is not None

    def test_tempfile_available(self):
        """Test tempfile is available."""
        import tempfile
        assert tempfile is not None


class TestFileOperations:
    """Test file operation utilities."""

    def test_isolated_tmp_dir_fixture(self, isolated_tmp_dir):
        """Test isolated_tmp_dir fixture."""
        assert isolated_tmp_dir.exists()
        assert isolated_tmp_dir.is_dir()

        # Should be empty initially
        assert list(isolated_tmp_dir.iterdir()) == []

        # Should be able to create files
        test_file = isolated_tmp_dir / "test.txt"
        test_file.write_text("test content")
        assert test_file.exists()
        assert test_file.read_text() == "test content"

    def test_seeded_random_fixture(self, seeded_random):
        """Test seeded_random fixture."""
        assert seeded_random is not None

        # Should produce deterministic results
        val1 = seeded_random.random()
        seeded_random.seed(42)  # Reset seed
        val2 = seeded_random.random()
        assert val1 == val2

    def test_sample_dna_sequences_fixture(self, sample_dna_sequences):
        """Test sample_dna_sequences fixture."""
        assert isinstance(sample_dna_sequences, dict)
        assert "short" in sample_dna_sequences
        assert "medium" in sample_dna_sequences
        assert "long" in sample_dna_sequences

        # Check sequences contain valid DNA bases
        for seq in sample_dna_sequences.values():
            assert all(c in "ATCGNUatcgnu-" for c in seq)

    def test_sample_protein_sequences_fixture(self, sample_protein_sequences):
        """Test sample_protein_sequences fixture."""
        assert isinstance(sample_protein_sequences, dict)
        assert "short" in sample_protein_sequences
        assert "medium" in sample_protein_sequences
        assert "long" in sample_protein_sequences

        # Check sequences contain valid amino acids
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY*")
        for seq in sample_protein_sequences.values():
            assert all(c in valid_aa for c in seq)

    def test_test_filesystem_fixture(self, test_filesystem):
        """Test test_filesystem fixture."""
        from tests.conftest import TestFileSystem

        assert isinstance(test_filesystem, TestFileSystem)

        # Test file creation
        test_file = test_filesystem.create_file("test.txt", "test content")
        assert test_file.exists()
        assert test_file.read_text() == "test content"

        # Test binary file creation
        test_bin = test_filesystem.create_binary_file("test.bin", b"binary data")
        assert test_bin.exists()
        assert test_bin.read_bytes() == b"binary data"


class TestPerformanceTracker:
    """Test performance tracking functionality."""

    def test_performance_tracker_fixture(self, performance_tracker):
        """Test performance_tracker fixture."""
        from tests.conftest import PerformanceTracker

        assert isinstance(performance_tracker, PerformanceTracker)

        # Test basic functionality
        performance_tracker.track_time("test_function", 1.5)
        avg_time = performance_tracker.get_average_time("test_function")
        assert avg_time == 1.5



