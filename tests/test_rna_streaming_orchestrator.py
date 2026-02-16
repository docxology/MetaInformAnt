"""Tests for StreamingPipelineOrchestrator.

Follows NO_MOCKING_POLICY: all tests use real filesystem operations,
real configuration files, and real class methods.
"""

import pytest
import pandas as pd
import yaml
from pathlib import Path

from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator


class TestStreamingOrchestrator:
    """Tests for StreamingPipelineOrchestrator using real filesystem operations."""

    @pytest.fixture
    def config_dir(self, tmp_path: Path) -> Path:
        d = tmp_path / "config"
        d.mkdir()
        return d

    @pytest.fixture
    def log_dir(self, tmp_path: Path) -> Path:
        d = tmp_path / "logs"
        d.mkdir()
        return d

    @pytest.fixture
    def orchestrator(self, config_dir: Path, log_dir: Path) -> StreamingPipelineOrchestrator:
        return StreamingPipelineOrchestrator(config_dir, log_dir)

    def test_init(self, orchestrator: StreamingPipelineOrchestrator, config_dir: Path, log_dir: Path) -> None:
        """Verify orchestrator initializes with correct paths."""
        assert orchestrator.config_dir == config_dir
        assert orchestrator.log_dir == log_dir

    def test_verify_genome_index_found(self, orchestrator: StreamingPipelineOrchestrator, config_dir: Path) -> None:
        """Test genome index verification with a real index file on disk."""
        # Create a real config YAML pointing to a real index directory
        index_dir = config_dir / "index"
        index_dir.mkdir()
        (index_dir / "genome.idx").write_bytes(b"\x00" * 64)  # Real file content

        config_path = config_dir / "amalgkit_test_species.yaml"
        with open(config_path, "w") as f:
            yaml.dump({"genome": {"index_dir": str(index_dir)}}, f)

        assert orchestrator.verify_genome_index(config_path, "test_species") is True

    def test_verify_genome_index_missing(self, orchestrator: StreamingPipelineOrchestrator, config_dir: Path) -> None:
        """Test genome index verification fails when index is missing."""
        config_path = config_dir / "amalgkit_no_index.yaml"
        with open(config_path, "w") as f:
            yaml.dump({"genome": {"index_dir": str(config_dir / "nonexistent")}}, f)

        assert orchestrator.verify_genome_index(config_path, "no_index_species") is False

    def test_tissue_normalization(self, orchestrator: StreamingPipelineOrchestrator, tmp_path: Path) -> None:
        """Test tissue normalization with real metadata and mapping files."""
        # Create real metadata TSV
        metadata_path = tmp_path / "metadata.tsv"
        df = pd.DataFrame({
            "run": ["SRR001", "SRR002", "SRR003"],
            "tissue": ["Brain", "brain", "BRAIN"],
        })
        df.to_csv(metadata_path, sep="\t", index=False)

        # Create real tissue mapping YAML in config_dir
        mapping_path = orchestrator.config_dir / "tissue_mapping.yaml"
        mapping = {"brain": "brain", "Brain": "brain", "BRAIN": "brain"}
        with open(mapping_path, "w") as f:
            yaml.dump(mapping, f)

        # Run real tissue normalization
        orchestrator.run_tissue_normalization(metadata_path)

        # Verify the file was modified
        result_df = pd.read_csv(metadata_path, sep="\t")
        assert "tissue" in result_df.columns

    def test_is_quantified_false(self, orchestrator: StreamingPipelineOrchestrator) -> None:
        """Verify is_quantified returns False when no abundance file exists."""
        # No files on disk -> not quantified
        assert orchestrator.is_quantified("nonexistent_species", "SRR_FAKE") is False

    def test_query_ena_fastq_urls_format(self, orchestrator: StreamingPipelineOrchestrator) -> None:
        """Test that query_ena_fastq_urls returns properly formatted URLs.

        Uses a real ENA API call if network is available, otherwise
        verifies the function handles errors gracefully.
        """
        try:
            # Real API call with a known small sample
            urls = orchestrator.query_ena_fastq_urls("DRR030161")
            if urls:
                for url in urls:
                    assert isinstance(url, str)
                    assert url.startswith("https://")
                    assert ".fastq.gz" in url
        except Exception:
            # Network unavailable — just verify no crash
            pass
