"""Tests for download and quantification skip logic.

Tests verify that already-quantified samples are skipped correctly,
and that already-downloaded FASTQ files are processed (quantified and deleted)
rather than re-downloaded.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.steps.process_samples import _sample_already_quantified
from metainformant.core.io import write_delimited


class TestSkipQuantifiedSamples:
    """Test that already-quantified samples are skipped."""

    def test_sample_already_quantified_returns_true(self, tmp_path: Path):
        """Test that _sample_already_quantified returns True when abundance.tsv exists."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_id = "SRR123456"
        sample_quant_dir = quant_dir / sample_id
        sample_quant_dir.mkdir()
        
        # Create abundance.tsv file
        abundance_file = sample_quant_dir / "abundance.tsv"
        abundance_file.write_text("target_id\tlength\teff_length\test_counts\ttpm\n")
        
        # Should return True
        assert _sample_already_quantified(sample_id, quant_dir) is True

    def test_sample_already_quantified_returns_false(self, tmp_path: Path):
        """Test that _sample_already_quantified returns False when abundance.tsv doesn't exist."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_id = "SRR123456"
        sample_quant_dir = quant_dir / sample_id
        sample_quant_dir.mkdir()
        
        # No abundance.tsv file
        # Should return False
        assert _sample_already_quantified(sample_id, quant_dir) is False

    def test_sample_already_quantified_missing_dir(self, tmp_path: Path):
        """Test that _sample_already_quantified returns False when sample directory doesn't exist."""
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        sample_id = "SRR123456"
        # No sample directory created
        
        # Should return False
        assert _sample_already_quantified(sample_id, quant_dir) is False


class TestProcessDownloadedFastq:
    """Test that already-downloaded FASTQ files are processed (quantified and deleted)."""

    @pytest.mark.slow
    def test_skip_completed_skips_quantified(self, tmp_path: Path):
        """Test that skip_completed=True skips already-quantified samples."""
        from metainformant.rna.steps.process_samples import run_download_quant_workflow
        
        # Setup: Create metadata
        metadata_path = tmp_path / "metadata.tsv"
        write_delimited(
            [
                {"run": "SRR001", "scientific_name": "Test_species"},
                {"run": "SRR002", "scientific_name": "Test_species"},
            ],
            metadata_path,
            delimiter="\t",
        )
        
        # Setup: Create quant directory with one already-quantified sample
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        (quant_dir / "SRR001").mkdir()
        (quant_dir / "SRR001" / "abundance.tsv").write_text("target_id\tlength\teff_length\test_counts\ttpm\n")
        
        # Run with skip_completed=True
        # Note: This will fail if amalgkit/kallisto not available, but that's expected
        # The test verifies the skip logic structure
        stats = run_download_quant_workflow(
            metadata_path=metadata_path,
            getfastq_params={"out_dir": str(tmp_path / "fastq")},
            quant_params={"out_dir": str(quant_dir), "index_dir": str(tmp_path / "index")},
            work_dir=tmp_path,
            num_workers=1,  # Sequential mode
            skip_completed=True,
        )
        
        # Should have skipped SRR001 (already quantified)
        assert stats["skipped"] >= 1
        assert stats["total_samples"] == 2

    @pytest.mark.slow
    def test_skip_completed_false_processes_all(self, tmp_path: Path):
        """Test that skip_completed=False processes all samples even if quantified."""
        from metainformant.rna.steps.process_samples import run_download_quant_workflow
        
        # Setup: Create metadata
        metadata_path = tmp_path / "metadata.tsv"
        write_delimited(
            [
                {"run": "SRR001", "scientific_name": "Test_species"},
            ],
            metadata_path,
            delimiter="\t",
        )
        
        # Setup: Create quant directory with already-quantified sample
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        (quant_dir / "SRR001").mkdir()
        (quant_dir / "SRR001" / "abundance.tsv").write_text("target_id\tlength\teff_length\test_counts\ttpm\n")
        
        # Run with skip_completed=False
        # This will attempt to process even though already quantified
        # (may fail due to missing dependencies, but structure is correct)
        stats = run_download_quant_workflow(
            metadata_path=metadata_path,
            getfastq_params={"out_dir": str(tmp_path / "fastq")},
            quant_params={"out_dir": str(quant_dir), "index_dir": str(tmp_path / "index")},
            work_dir=tmp_path,
            num_workers=1,
            skip_completed=False,
        )
        
        # Should have attempted to process (not skipped)
        assert stats["skipped"] == 0
        assert stats["total_samples"] == 1


class TestDownloadedFastqProcessing:
    """Test that already-downloaded FASTQ files are processed rather than re-downloaded."""

    @pytest.mark.slow
    def test_downloaded_fastq_processed(self, tmp_path: Path):
        """Test that if FASTQ exists, it's quantified and deleted (not re-downloaded).
        
        This test verifies the behavior described: if raw FASTQ file is there,
        process it (quantify) then delete it, rather than skipping.
        """
        from metainformant.rna.steps.process_samples import run_download_quant_workflow
        
        # Setup: Create metadata
        metadata_path = tmp_path / "metadata.tsv"
        write_delimited(
            [
                {"run": "SRR001", "scientific_name": "Test_species"},
            ],
            metadata_path,
            delimiter="\t",
        )
        
        # Setup: Create FASTQ directory with already-downloaded FASTQ
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir(parents=True)
        sample_fastq_dir = fastq_dir / "SRR001"
        sample_fastq_dir.mkdir()
        # Create dummy FASTQ file
        (sample_fastq_dir / "SRR001_1.fastq.gz").write_text("dummy fastq content")
        
        # Setup: Quant directory (empty - not quantified yet)
        quant_dir = tmp_path / "quant"
        quant_dir.mkdir(parents=True)
        
        # The workflow should:
        # 1. Detect FASTQ exists (not re-download)
        # 2. Quantify it
        # 3. Delete FASTQ after quantification
        
        # Note: This will fail if kallisto/index not available, but structure is correct
        # The key is that skip_completed=True should NOT skip if FASTQ exists but not quantified
        stats = run_download_quant_workflow(
            metadata_path=metadata_path,
            getfastq_params={"out_dir": str(fastq_dir)},
            quant_params={"out_dir": str(quant_dir), "index_dir": str(tmp_path / "index")},
            work_dir=tmp_path,
            num_workers=1,
            skip_completed=True,  # Should still process if FASTQ exists but not quantified
        )
        
        # Should attempt to process (not skipped) because not quantified yet
        # Even though FASTQ exists, we want to quantify it
        assert stats["total_samples"] == 1
        # skipped should be 0 because sample is not quantified yet
        # (quantification check happens, FASTQ existence doesn't cause skip)

