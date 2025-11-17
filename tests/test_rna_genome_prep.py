"""Tests for RNA genome preparation functions.

This module tests genome preparation functionality following NO_MOCKING_POLICY.
All tests use real file operations.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.rna import genome_prep


class TestFindRnaFastaInGenomeDir:
    """Test find_rna_fasta_in_genome_dir function."""

    def test_find_rna_fasta_not_found(self, tmp_path: Path):
        """Test find_rna_fasta_in_genome_dir when file not found."""
        result = genome_prep.find_rna_fasta_in_genome_dir(tmp_path, "GCF_123456789")
        assert result is None

    def test_find_rna_fasta_api_pattern(self, tmp_path: Path):
        """Test find_rna_fasta_in_genome_dir finds API extraction pattern."""
        accession = "GCF_123456789"
        api_dir = tmp_path / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession
        api_dir.mkdir(parents=True)
        rna_file = api_dir / "rna.fna"
        rna_file.write_text(">test\nATCG")
        
        result = genome_prep.find_rna_fasta_in_genome_dir(tmp_path, accession)
        assert result == rna_file

    def test_find_rna_fasta_cli_pattern(self, tmp_path: Path):
        """Test find_rna_fasta_in_genome_dir finds CLI extraction pattern."""
        accession = "GCF_123456789"
        cli_dir = tmp_path / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession
        cli_dir.mkdir(parents=True)
        rna_file = cli_dir / "rna.fna"
        rna_file.write_text(">test\nATCG")
        
        result = genome_prep.find_rna_fasta_in_genome_dir(tmp_path, accession)
        assert result == rna_file

    def test_find_rna_fasta_direct_pattern(self, tmp_path: Path):
        """Test find_rna_fasta_in_genome_dir finds direct file pattern."""
        accession = "GCF_123456789"
        rna_file = tmp_path / "rna.fna"
        rna_file.write_text(">test\nATCG")
        
        result = genome_prep.find_rna_fasta_in_genome_dir(tmp_path, accession)
        assert result == rna_file

    def test_find_rna_fasta_prefers_unzipped(self, tmp_path: Path):
        """Test find_rna_fasta_in_genome_dir prefers unzipped files."""
        accession = "GCF_123456789"
        api_dir = tmp_path / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession
        api_dir.mkdir(parents=True)
        
        # Create both zipped and unzipped
        rna_gz = api_dir / "rna.fna.gz"
        rna_unzipped = api_dir / "rna.fna"
        rna_unzipped.write_text(">test\nATCG")
        
        with gzip.open(rna_gz, "wt") as f:
            f.write(">test\nATCG")
        
        result = genome_prep.find_rna_fasta_in_genome_dir(tmp_path, accession)
        assert result == rna_unzipped


class TestGetExpectedIndexPath:
    """Test get_expected_index_path function."""

    def test_get_expected_index_path(self, tmp_path: Path):
        """Test get_expected_index_path returns correct path."""
        work_dir = tmp_path / "work"
        species_name = "Test_species"
        
        expected = genome_prep.get_expected_index_path(work_dir, species_name)
        assert isinstance(expected, Path)
        assert expected.name.endswith(".idx")
        assert species_name in str(expected)


class TestGenomePrepDocumentation:
    """Test that genome_prep functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all genome_prep functions have docstrings."""
        functions = [
            genome_prep.find_rna_fasta_in_genome_dir,
            genome_prep.download_rna_fasta_from_ftp,
            genome_prep.download_cds_fasta_from_ftp,
            genome_prep.extract_transcripts_from_gff,
            genome_prep.prepare_transcriptome_for_kallisto,
            genome_prep.build_kallisto_index,
            genome_prep.get_expected_index_path,
            genome_prep.prepare_genome_for_quantification,
            genome_prep.verify_genome_status,
            genome_prep.orchestrate_genome_setup,
        ]
        
        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0


class TestGenomePrepErrorHandling:
    """Test error handling in genome_prep functions."""

    def test_find_rna_fasta_handles_missing_dir(self):
        """Test find_rna_fasta_in_genome_dir handles missing directory."""
        missing_dir = Path("/nonexistent/path/to/genome")
        result = genome_prep.find_rna_fasta_in_genome_dir(missing_dir, "GCF_123")
        assert result is None

    def test_get_expected_index_path_creates_path(self, tmp_path: Path):
        """Test get_expected_index_path creates path structure."""
        work_dir = tmp_path / "work"
        species_name = "Test_species"
        
        expected = genome_prep.get_expected_index_path(work_dir, species_name)
        # Path should be valid even if directory doesn't exist
        assert isinstance(expected, Path)
        assert expected.parent.name == "index" or "index" in str(expected)

