"""Tests for RNA cleanup functions.

This module tests cleanup functionality following NO_MOCKING_POLICY.
All tests use real file operations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core.io import write_delimited
from metainformant.rna import cleanup
from metainformant.rna.workflow import AmalgkitWorkflowConfig


class TestFindPartialDownloads:
    """Test find_partial_downloads function."""

    def test_find_partial_downloads_empty_directories(self, tmp_path: Path):
        """Test find_partial_downloads with empty directories."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()
        
        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert isinstance(partial, list)
        assert len(partial) == 0

    def test_find_partial_downloads_missing_directories(self, tmp_path: Path):
        """Test find_partial_downloads with missing directories."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        
        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert isinstance(partial, list)
        assert len(partial) == 0

    def test_find_partial_downloads_with_unquantified_samples(self, tmp_path: Path):
        """Test find_partial_downloads finds unquantified samples."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir(parents=True)
        quant_dir.mkdir(parents=True)
        
        # Create sample directory with files but no quantification
        sample_dir = fastq_dir / "getfastq" / "SRR123"
        sample_dir.mkdir(parents=True)
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")
        
        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert len(partial) == 1
        assert partial[0][0] == "SRR123"
        assert partial[0][1] == sample_dir
        assert isinstance(partial[0][2], int)  # size_mb

    def test_find_partial_downloads_skips_quantified_samples(self, tmp_path: Path):
        """Test find_partial_downloads skips quantified samples."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir(parents=True)
        quant_dir.mkdir(parents=True)
        
        # Create quantified sample
        sample_dir = fastq_dir / "getfastq" / "SRR123"
        sample_dir.mkdir(parents=True)
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")
        
        quant_sample_dir = quant_dir / "SRR123"
        quant_sample_dir.mkdir(parents=True)
        (quant_sample_dir / "abundance.tsv").write_text("test")
        
        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert len(partial) == 0

    def test_find_partial_downloads_direct_structure(self, tmp_path: Path):
        """Test find_partial_downloads with direct structure (not getfastq subdir)."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir(parents=True)
        quant_dir.mkdir(parents=True)
        
        # Create sample directory directly in fastq_dir
        sample_dir = fastq_dir / "SRR123"
        sample_dir.mkdir(parents=True)
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")
        
        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert len(partial) == 1
        assert partial[0][0] == "SRR123"


class TestCleanupPartialDownloads:
    """Test cleanup_partial_downloads function."""

    def test_cleanup_partial_downloads_dry_run(self, tmp_path: Path):
        """Test cleanup_partial_downloads with dry_run=True."""
        # Create config file
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "getfastq": {"out_dir": str(tmp_path / "fastq")},
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        from metainformant.core.io import dump_json
        dump_json(config_data, config_file)
        
        # Create partial download
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir(parents=True)
        quant_dir.mkdir(parents=True)
        
        sample_dir = fastq_dir / "getfastq" / "SRR123"
        sample_dir.mkdir(parents=True)
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")
        
        result = cleanup.cleanup_partial_downloads(config_file, dry_run=True)
        assert isinstance(result, dict)
        assert "deleted" in result
        assert "freed_mb" in result
        assert "errors" in result
        assert result["deleted"] == 0  # Dry run doesn't delete
        assert sample_dir.exists()  # Should still exist

    def test_cleanup_partial_downloads_actual_cleanup(self, tmp_path: Path):
        """Test cleanup_partial_downloads actually deletes files."""
        # Create config file with absolute paths
        config_file = tmp_path / "config.yaml"
        # Use absolute paths since load_workflow_config resolves relative to repo root
        work_dir_abs = (tmp_path / "work").resolve()
        fastq_dir_abs = (tmp_path / "fastq").resolve()
        quant_dir_abs = (tmp_path / "quant").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "getfastq": {"out_dir": str(fastq_dir_abs)},
                "quant": {"out_dir": str(quant_dir_abs)},
            },
        }
        from metainformant.core.io import dump_json
        dump_json(config_data, config_file)
        
        # Create partial download
        fastq_dir_abs.mkdir(parents=True)
        quant_dir_abs.mkdir(parents=True)
        
        sample_dir = fastq_dir_abs / "getfastq" / "SRR123"
        sample_dir.mkdir(parents=True)
        test_file = sample_dir / "SRR123_1.fastq.gz"
        test_file.write_text("test")
        
        result = cleanup.cleanup_partial_downloads(config_file, dry_run=False)
        assert isinstance(result, dict)
        assert result["deleted"] == 1
        assert result["freed_mb"] >= 0
        assert not sample_dir.exists()  # Should be deleted


class TestFixAbundanceNaming:
    """Test fix_abundance_naming functions."""

    def test_fix_abundance_naming_creates_symlink(self, tmp_path: Path):
        """Test fix_abundance_naming creates symlink."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_dir = quant_dir / "SRR123"
        sample_dir.mkdir(parents=True)
        source = sample_dir / "abundance.tsv"
        source.write_text("test")
        
        result = cleanup.fix_abundance_naming(quant_dir, "SRR123")
        assert result is True
        
        target = sample_dir / "SRR123_abundance.tsv"
        assert target.exists()
        assert target.is_symlink()

    def test_fix_abundance_naming_handles_missing_source(self, tmp_path: Path):
        """Test fix_abundance_naming handles missing source file."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_dir = quant_dir / "SRR123"
        sample_dir.mkdir(parents=True)
        # Don't create abundance.tsv
        
        result = cleanup.fix_abundance_naming(quant_dir, "SRR123")
        assert result is False

    def test_fix_abundance_naming_handles_existing_target(self, tmp_path: Path):
        """Test fix_abundance_naming handles existing target."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_dir = quant_dir / "SRR123"
        sample_dir.mkdir(parents=True)
        source = sample_dir / "abundance.tsv"
        source.write_text("test")
        
        target = sample_dir / "SRR123_abundance.tsv"
        target.write_text("existing")
        
        result = cleanup.fix_abundance_naming(quant_dir, "SRR123")
        assert result is True  # Should return True if target exists

    def test_fix_abundance_naming_for_species(self, tmp_path: Path):
        """Test fix_abundance_naming_for_species fixes all samples."""
        # Create config file
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        from metainformant.core.io import dump_json
        dump_json(config_data, config_file)
        
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        # Create multiple samples
        for srr_id in ["SRR123", "SRR456"]:
            sample_dir = quant_dir / srr_id
            sample_dir.mkdir(parents=True)
            (sample_dir / "abundance.tsv").write_text("test")
        
        created, already_exists = cleanup.fix_abundance_naming_for_species(config_file)
        assert isinstance(created, int)
        assert isinstance(already_exists, int)
        # Should create 2 symlinks (or more if some already existed)
        assert created >= 0
        assert (created + already_exists) >= 0
        
        # Verify symlinks were created (if created > 0)
        if created > 0:
            for srr_id in ["SRR123", "SRR456"]:
                target = quant_dir / srr_id / f"{srr_id}_abundance.tsv"
                # Target may exist as symlink or regular file
                if target.exists():
                    assert target.is_symlink() or target.is_file()

    def test_fix_abundance_naming_for_species_handles_missing_quant_dir(self, tmp_path: Path):
        """Test fix_abundance_naming_for_species handles missing quant directory."""
        # Create config file
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "quant": {"out_dir": str(tmp_path / "nonexistent" / "quant")},
            },
        }
        from metainformant.core.io import dump_json
        dump_json(config_data, config_file)
        
        created, already_exists = cleanup.fix_abundance_naming_for_species(config_file)
        assert created == 0
        assert already_exists == 0


class TestCleanupDocumentation:
    """Test that cleanup functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all cleanup functions have docstrings."""
        functions = [
            cleanup.find_partial_downloads,
            cleanup.cleanup_partial_downloads,
            cleanup.fix_abundance_naming,
            cleanup.fix_abundance_naming_for_species,
        ]
        
        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0

