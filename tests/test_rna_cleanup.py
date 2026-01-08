"""Tests for RNA cleanup functions.

This module tests cleanup functionality following NO_MOCKING_POLICY.
All tests use real file operations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna import cleanup


class TestFindPartialDownloads:
    """Test find_partial_downloads function."""

    def test_empty_directories(self, tmp_path: Path):
        """Test with empty directories."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert partial == []

    def test_missing_directories(self, tmp_path: Path):
        """Test with missing directories."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"

        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert partial == []

    def test_finds_unquantified(self, tmp_path: Path):
        """Test finds samples without quantification."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        sample_dir = fastq_dir / "SRR123"
        sample_dir.mkdir()
        (sample_dir / "SRR123_1.fastq.gz").write_text("test data")

        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert len(partial) == 1
        assert partial[0][0] == "SRR123"

    def test_skips_quantified(self, tmp_path: Path):
        """Test skips already quantified samples."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        # Sample with FASTQ
        sample_fastq = fastq_dir / "SRR123"
        sample_fastq.mkdir()
        (sample_fastq / "SRR123_1.fastq.gz").write_text("test")

        # Sample quantified
        sample_quant = quant_dir / "SRR123"
        sample_quant.mkdir()
        (sample_quant / "abundance.tsv").write_text("test")

        partial = cleanup.find_partial_downloads(fastq_dir, quant_dir)
        assert partial == []


class TestCleanupPartialDownloads:
    """Test cleanup_partial_downloads function."""

    def test_dry_run(self, tmp_path: Path):
        """Test dry run doesn't delete."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        sample_dir = fastq_dir / "SRR123"
        sample_dir.mkdir()
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")

        result = cleanup.cleanup_partial_downloads(tmp_path, dry_run=True)
        assert result["deleted"] == 0
        assert sample_dir.exists()

    def test_actual_cleanup(self, tmp_path: Path):
        """Test actual deletion."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        sample_dir = fastq_dir / "SRR123"
        sample_dir.mkdir()
        (sample_dir / "SRR123_1.fastq.gz").write_text("test")

        result = cleanup.cleanup_partial_downloads(tmp_path, dry_run=False)
        assert result["deleted"] == 1
        assert not sample_dir.exists()


class TestFixAbundanceNaming:
    """Test fix_abundance_naming function."""

    def test_creates_symlink(self, tmp_path: Path):
        """Test symlink creation."""
        quant_dir = tmp_path / "quant"
        sample_dir = quant_dir / "SRR123"
        sample_dir.mkdir(parents=True)
        (sample_dir / "abundance.tsv").write_text("test")

        result = cleanup.fix_abundance_naming(quant_dir, "SRR123")
        assert result is True
        assert (sample_dir / "SRR123_abundance.tsv").is_symlink()

    def test_handles_missing_source(self, tmp_path: Path):
        """Test returns False if source missing."""
        quant_dir = tmp_path / "quant"
        sample_dir = quant_dir / "SRR123"
        sample_dir.mkdir(parents=True)

        result = cleanup.fix_abundance_naming(quant_dir, "SRR123")
        assert result is False


class TestCleanupUnquantifiedSamples:
    """Test cleanup_unquantified_samples function."""

    def test_removes_unquantified(self, tmp_path: Path):
        """Test removing unquantified samples."""
        fastq_dir = tmp_path / "fastq"
        quant_dir = tmp_path / "quant"
        fastq_dir.mkdir()
        quant_dir.mkdir()

        # Quantified sample (keep)
        (quant_dir / "S1").mkdir()
        (quant_dir / "S1" / "abundance.tsv").write_text("data")
        (fastq_dir / "S1_1.fastq").write_text("data")

        # Unquantified sample (remove)
        (fastq_dir / "S2_1.fastq").write_text("data")

        cleaned = cleanup.cleanup_unquantified_samples(tmp_path)
        assert "S2" in cleaned
        assert not (fastq_dir / "S2_1.fastq").exists()
        assert (fastq_dir / "S1_1.fastq").exists()


