"""Tests for RNA pipeline module.

All tests follow real-implementation policy and use real implementations.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.engine.pipeline import RNAPipelineConfig, summarize_finalize_tables


class TestRNAPipelineConfig:
    """Test RNAPipelineConfig dataclass."""

    def test_config_creation(self, tmp_path: Path):
        """Test creating pipeline configuration."""
        config = RNAPipelineConfig(work_dir=tmp_path / "work")
        assert config.work_dir == tmp_path / "work"
        assert isinstance(config.work_dir, Path)

    def test_config_work_dir_path(self, tmp_path: Path):
        """Test that work_dir can be any Path."""
        work_dir = tmp_path / "pipeline" / "output"
        config = RNAPipelineConfig(work_dir=work_dir)
        assert config.work_dir == work_dir


class TestSummarizeFinalizeTables:
    """Test summarize_finalize_tables function."""

    def test_summarize_empty_directory(self, tmp_path: Path):
        """Test summarizing empty directory."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        counts = summarize_finalize_tables(empty_dir)
        assert counts["_total"] == 0
        assert "_error" not in counts

    def test_summarize_nonexistent_directory(self, tmp_path: Path):
        """Test summarizing non-existent directory."""
        nonexistent = tmp_path / "does_not_exist"
        counts = summarize_finalize_tables(nonexistent)
        assert counts["_total"] == 0
        assert "_error" in counts

    def test_summarize_with_tsv_files(self, tmp_path: Path):
        """Test summarizing directory with TSV files."""
        finalize_dir = tmp_path / "finalize"
        finalize_dir.mkdir()

        # Create some TSV files
        (finalize_dir / "summary.tsv").write_text("test")
        (finalize_dir / "metadata.tsv").write_text("test")
        (finalize_dir / "summary.tsv").write_text("test")  # Duplicate name

        # Create subdirectory with TSV
        subdir = finalize_dir / "subdir"
        subdir.mkdir()
        (subdir / "data.tsv").write_text("test")

        counts = summarize_finalize_tables(finalize_dir)
        assert counts["summary.tsv"] == 1  # Counts by filename, not path
        assert counts["metadata.tsv"] == 1
        assert counts["data.tsv"] == 1

    def test_summarize_ignores_non_tsv(self, tmp_path: Path):
        """Test that non-TSV files are ignored."""
        finalize_dir = tmp_path / "finalize"
        finalize_dir.mkdir()

        (finalize_dir / "data.tsv").write_text("test")
        (finalize_dir / "data.txt").write_text("test")
        (finalize_dir / "data.csv").write_text("test")

        counts = summarize_finalize_tables(finalize_dir)
        assert "data.tsv" in counts
        assert "data.txt" not in counts
        assert "data.csv" not in counts
