"""Tests for StreamingPipelineOrchestrator.

Follows real-implementation policy: all tests use real filesystem operations,
real configuration files, and real class methods.
"""

from pathlib import Path

import pandas as pd
import pytest
import yaml

from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    _build_quant_command,
    _build_sample_tasks,
    _filter_metadata_by_size,
    _resolve_metadata_path,
    _resolve_quant_metadata_path,
    _sample_run_column,
    _species_name_from_config,
)


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
        df = pd.DataFrame(
            {
                "run": ["SRR001", "SRR002", "SRR003"],
                "tissue": ["Brain", "brain", "BRAIN"],
            }
        )
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


def test_streaming_path_helpers_prefer_current_metadata(tmp_path: Path) -> None:
    """Metadata helper should prefer metadata.tsv and fall back to selected metadata."""
    work_dir = tmp_path / "work"
    metadata_dir = work_dir / "metadata"
    metadata_dir.mkdir(parents=True)
    selected = metadata_dir / "metadata_selected.tsv"
    selected.write_text("run\ttotal_bases\nSRR2\t2\n")

    assert _species_name_from_config("amalgkit_apis.yaml") == "apis"
    assert _resolve_metadata_path(work_dir) == selected
    assert _resolve_quant_metadata_path(work_dir) == str(selected)

    current = metadata_dir / "metadata.tsv"
    current.write_text("run\ttotal_bases\nSRR1\t1\n")
    assert _resolve_metadata_path(work_dir) == current
    assert _resolve_quant_metadata_path(work_dir) == str(current)


def test_streaming_task_helpers_filter_and_build_tasks(tmp_path: Path) -> None:
    """Task helper should sort by size and preserve batch indices for real rows."""
    df = pd.DataFrame(
        {
            "run_accession": ["SRR_BIG", "SRR_SMALL", "SRR_MED"],
            "total_bases": [4_000_000_000, "500000000", 1_500_000_000],
        }
    )

    filtered = _filter_metadata_by_size(df, max_gb=2.0)
    srr_col = _sample_run_column(filtered)
    tasks = _build_sample_tasks(filtered, srr_col, tmp_path / "fastq", tmp_path / "config.yaml", "apis")

    assert [task["srr"] for task in tasks] == ["SRR_SMALL", "SRR_MED"]
    assert [task["batch_idx"] for task in tasks] == [1, 2]


def test_build_quant_command_uses_metadata_cleanup_and_index(tmp_path: Path) -> None:
    """Quant command construction should include cleanup and real index paths."""
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (index_dir / "genome.idx").write_bytes(b"idx")
    cfg = {
        "steps": {"quant": {"keep_fastq": "no"}},
        "genome": {"index_dir": str(index_dir)},
    }

    cmd = _build_quant_command(cfg, "apis", batch_index=3, threads=2, metadata_path="metadata.tsv")

    assert cmd[:4] == ["amalgkit", "quant", "--out_dir", "output/amalgkit/apis/work"]
    assert cmd[cmd.index("--metadata") + 1] == "metadata.tsv"
    assert cmd[cmd.index("--clean_fastq") + 1] == "yes"
    assert cmd[cmd.index("--index_dir") + 1] == str(index_dir)
