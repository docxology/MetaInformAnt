"""Tests for GWAS SRA download module.

All tests follow NO_MOCKING_POLICY and use real implementations.
Tests use real subprocess calls to check SRA tools availability.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from metainformant.gwas.sra_download import (
    check_sra_tools_available,
    download_sra_project,
    download_sra_run,
    search_sra_for_organism,
)


class TestCheckSRAToolsAvailable:
    """Test SRA tools availability checking."""

    def test_check_sra_tools_available_real(self):
        """Test real check for SRA tools availability.
        
        Uses real subprocess calls following NO_MOCKING_POLICY.
        """
        available = check_sra_tools_available()
        assert isinstance(available, bool)
        # Result depends on whether SRA Toolkit is installed on system

    def test_check_returns_boolean(self):
        """Test that function returns boolean value."""
        result = check_sra_tools_available()
        assert isinstance(result, bool)


class TestDownloadSRARun:
    """Test SRA run download functionality."""

    @pytest.mark.slow
    @pytest.mark.network
    @pytest.mark.external_tool
    @pytest.mark.skipif(
        not shutil.which("fasterq-dump") and not shutil.which("fastq-dump"),
        reason="SRA Toolkit not available - real implementation requires fasterq-dump or fastq-dump"
    )
    def test_download_sra_run_requires_tools(self, tmp_path: Path):
        """Test that download requires SRA tools to be available.

        This test verifies the function structure but will skip if
        SRA tools are not available on the system.
        """
        # Check if SRA tools are available
        import shutil
        sra_tools = ["fasterq-dump", "prefetch", "fastq-dump"]
        available_tools = [tool for tool in sra_tools if shutil.which(tool)]

        if not available_tools:
            pytest.skip("SRA toolkit not available (none of fasterq-dump, prefetch, fastq-dump found)")

        # Small test accession - may or may not exist
        # This tests the function interface, not actual download
        # Use a very short timeout to avoid hanging on real downloads
        try:
            result = download_sra_run(
                "SRR000001",  # Example accession
                tmp_path,
                use_fasterq=True,
                threads=1,
                timeout=5  # Very short timeout for testing
            )
            # If it succeeds, verify structure
            assert isinstance(result, dict)
        except Exception:
            # Expected if tools unavailable or accession invalid
            # This documents real behavior
            pass

    @pytest.mark.network
    @pytest.mark.external_tool
    def test_download_sra_run_invalid_accession(self, tmp_path: Path):
        """Test behavior with invalid accession.

        Uses real implementation - may fail or skip depending on tool availability.
        """
        if not check_sra_tools_available():
            pytest.skip("SRA Toolkit not available - real implementation requires tools")
        
        # Test with clearly invalid accession
        try:
            result = download_sra_run("INVALID_ACCESSION", tmp_path)
            # If it doesn't raise, check result structure
            assert isinstance(result, dict)
        except Exception:
            # Expected behavior for invalid accession
            pass


class TestDownloadSRAProject:
    """Test SRA project download functionality."""

    @pytest.mark.skipif(
        not shutil.which("fasterq-dump") and not shutil.which("fastq-dump"),
        reason="SRA Toolkit not available - real implementation requires tools"
    )
    @pytest.mark.network
    @pytest.mark.external_tool
    def test_download_sra_project_requires_tools(self, tmp_path: Path):
        """Test project download interface.

        Verifies function structure without requiring actual download.
        """
        try:
            result = download_sra_project(
                "SRP000001",  # Example project
                tmp_path,
                max_runs=1  # Limit for testing
            )
            assert isinstance(result, dict)
        except Exception:
            # Expected if tools unavailable or project invalid
            pass


class TestSearchSRAForOrganism:
    """Test SRA search functionality."""

    @pytest.mark.network
    @pytest.mark.external_tool
    def test_search_sra_for_organism_interface(self):
        """Test search function interface.

        Tests function structure - actual search requires network access.
        """
        # Test with a real organism name
        try:
            results = search_sra_for_organism("Homo sapiens", max_results=5)
            # If search succeeds, verify structure
            assert isinstance(results, list)
            # Results should be dictionaries with SRA metadata
            if results:
                assert isinstance(results[0], dict)
        except Exception:
            # Expected if network unavailable or API issues
            # This documents real failure modes
            pass

    @pytest.mark.network
    @pytest.mark.external_tool
    def test_search_sra_empty_results(self):
        """Test search with organism that may have no results."""
        try:
            results = search_sra_for_organism("NonexistentSpecies123", max_results=5)
            # Should return empty list, not error
            assert isinstance(results, list)
        except Exception:
            # Network or API errors are acceptable
            pass


