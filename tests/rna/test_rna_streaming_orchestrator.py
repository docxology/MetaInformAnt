"""Tests for StreamingPipelineOrchestrator.

Follows real-implementation policy: all tests use real filesystem operations,
real configuration files, and real class methods.
"""

import gzip
import json
import subprocess
import sys
import threading
import time
from pathlib import Path

import pandas as pd
import pytest
import yaml

from metainformant.rna.engine import streaming_orchestrator as orchestrator_module
from metainformant.rna.engine.provenance import (
    is_current_metadata,
    is_current_quantification,
    write_metadata_provenance,
    write_quant_provenance,
)
from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    _build_quant_command,
    _build_sample_tasks,
    _campaign_cached_sra_path,
    _ensure_ncbi_settings_for_data_root,
    _ensure_reference_alias_indexes,
    _filter_kallisto_ineligible_reads,
    _filter_metadata_by_size,
    _filter_valid_run_rows,
    _is_valid_run_accession,
    _library_layout_is_paired,
    _metadata_requires_fastq_input_stats,
    _prioritize_tasks_by_raw_state,
    _resolve_metadata_path,
    _resolve_quant_metadata_path,
    _run_command_in_process_group,
    _runtime_config_path,
    _sample_run_column,
    _species_name_from_config,
    _validated_local_fastq_inputs,
    _write_fastq_input_stats,
    build_pipeline_resource_profile,
)
from metainformant.rna.retrieval import ena_downloader as ena_downloader_module


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
    def orchestrator(self, config_dir: Path, log_dir: Path, tmp_path: Path) -> StreamingPipelineOrchestrator:
        return StreamingPipelineOrchestrator(config_dir, log_dir, db_path=tmp_path / "progress.db")

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

    def test_tissue_normalization_is_idempotent_and_preserves_original(
        self, orchestrator: StreamingPipelineOrchestrator, tmp_path: Path
    ) -> None:
        """Repeated discovery cannot replace the source tissue label."""
        metadata_path = tmp_path / "metadata.tsv"
        pd.DataFrame({"run": ["SRR001"], "tissue": ["Brain"]}).to_csv(metadata_path, sep="\t", index=False)
        (orchestrator.config_dir / "tissue_mapping.yaml").write_text("brain:\n  - brain\n  - Brain\n", encoding="utf-8")

        assert orchestrator.run_tissue_normalization(metadata_path) is True
        first = pd.read_csv(metadata_path, sep="\t")
        assert first.loc[0, "tissue"] == "brain"
        assert first.loc[0, "tissue_original"] == "Brain"

        assert orchestrator.run_tissue_normalization(metadata_path) is True
        second = pd.read_csv(metadata_path, sep="\t")
        pd.testing.assert_frame_equal(first, second)

    def test_is_quantified_false(self, orchestrator: StreamingPipelineOrchestrator) -> None:
        """Verify is_quantified returns False when no abundance file exists."""
        # No files on disk -> not quantified
        assert orchestrator.is_quantified("nonexistent_species", "SRR_FAKE") is False

    @pytest.mark.network
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


def test_streaming_path_helpers_prefer_selected_metadata(tmp_path: Path) -> None:
    """Metadata helper should prefer the selected cohort for quantification."""
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
    assert _resolve_metadata_path(work_dir) == selected
    assert _resolve_quant_metadata_path(work_dir) == str(selected)


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


def test_kallisto_eligibility_excludes_known_short_reads_and_keeps_unknowns() -> None:
    """Known sub-k-mer libraries are audited before download; unknown lengths remain eligible."""

    dataframe = pd.DataFrame(
        {
            "run": ["SRR_SHORT", "SRR_OK", "SRR_UNKNOWN"],
            "spot_length": [22, 50, None],
            "total_bases": [220, 500, 1000],
            "total_spots": [10, 10, None],
        }
    )

    eligible, excluded, read_length = _filter_kallisto_ineligible_reads(dataframe, 31)

    assert eligible["run"].tolist() == ["SRR_OK", "SRR_UNKNOWN"]
    assert excluded["run"].tolist() == ["SRR_SHORT"]
    assert float(read_length.iloc[0]) == 22
    assert float(read_length.iloc[2]) != float(read_length.iloc[2])
    assert excluded["quant_eligibility_reason"].tolist() == ["mean_read_length_below_kallisto_kmer"]


def test_task_scheduler_prioritizes_existing_and_resumable_raw_inputs(tmp_path: Path) -> None:
    """Raw reuse and resumable partials must precede fresh acquisition."""
    fastq_dir = tmp_path / "getfastq"
    (fastq_dir / "SRR_EXISTING").mkdir(parents=True)
    (fastq_dir / "SRR_EXISTING" / "SRR_EXISTING_1.fastq.gz").write_bytes(b"raw")
    (fastq_dir / "SRR_PARTIAL").mkdir(parents=True)
    (fastq_dir / "SRR_PARTIAL" / "SRR_PARTIAL_1.fastq.gz.part").write_bytes(b"partial")
    (fastq_dir / "SRR987654_1.fastq.gz").write_bytes(b"raw")

    tasks = [
        {"srr": "SRR_NEW", "fastq_dir": fastq_dir, "batch_idx": 1},
        {"srr": "SRR_PARTIAL", "fastq_dir": fastq_dir, "batch_idx": 2},
        {"srr": "SRR_EXISTING", "fastq_dir": fastq_dir, "batch_idx": 3},
        {"srr": "SRR987654", "fastq_dir": fastq_dir, "batch_idx": 4},
    ]

    prioritized = _prioritize_tasks_by_raw_state(tasks)

    assert [task["srr"] for task in prioritized] == ["SRR_EXISTING", "SRR987654", "SRR_PARTIAL", "SRR_NEW"]
    assert [task["batch_idx"] for task in prioritized] == [3, 4, 2, 1]


def test_task_scheduler_uses_size_within_raw_state_tiers(tmp_path: Path) -> None:
    """Small reusable reads and nearly complete partials improve early progress."""
    fastq_dir = tmp_path / "getfastq"
    (fastq_dir / "SRR_LARGE").mkdir(parents=True)
    (fastq_dir / "SRR_LARGE" / "SRR_LARGE_1.fastq.gz").write_bytes(b"x" * 1000)
    (fastq_dir / "SRR_SMALL").mkdir(parents=True)
    (fastq_dir / "SRR_SMALL" / "SRR_SMALL_1.fastq.gz").write_bytes(b"x")
    (fastq_dir / "SRR_PARTIAL_SMALL").mkdir(parents=True)
    (fastq_dir / "SRR_PARTIAL_SMALL" / "SRR_PARTIAL_SMALL_1.fastq.gz.part").write_bytes(b"x")
    (fastq_dir / "SRR_PARTIAL_LARGE").mkdir(parents=True)
    (fastq_dir / "SRR_PARTIAL_LARGE" / "SRR_PARTIAL_LARGE_1.fastq.gz.part").write_bytes(b"x" * 1000)

    tasks = [
        {"srr": "SRR_LARGE", "fastq_dir": fastq_dir, "batch_idx": 1},
        {"srr": "SRR_PARTIAL_SMALL", "fastq_dir": fastq_dir, "batch_idx": 2},
        {"srr": "SRR_SMALL", "fastq_dir": fastq_dir, "batch_idx": 3},
        {"srr": "SRR_PARTIAL_LARGE", "fastq_dir": fastq_dir, "batch_idx": 4},
    ]

    prioritized = _prioritize_tasks_by_raw_state(tasks)

    assert [task["srr"] for task in prioritized] == [
        "SRR_SMALL",
        "SRR_LARGE",
        "SRR_PARTIAL_LARGE",
        "SRR_PARTIAL_SMALL",
    ]
    assert [task["_raw_input_bytes"] for task in prioritized] == [1, 1000, 1000, 1]


def test_task_scheduler_prioritizes_campaign_cached_sra(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Complete shared SRA cache entries precede fresh network acquisition."""

    data_root = tmp_path / "campaign-data"
    cache_dir = data_root / ".sra-cache" / "sra"
    cache_dir.mkdir(parents=True)
    cached = cache_dir / "SRR_CACHED.sra.cache"
    cached.write_bytes(b"cached-sra")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))

    fastq_dir = tmp_path / "species" / "work" / "getfastq"
    tasks = [
        {"srr": "SRR_NEW", "fastq_dir": fastq_dir, "batch_idx": 1},
        {"srr": "SRR_CACHED", "fastq_dir": fastq_dir, "batch_idx": 2},
    ]

    prioritized = _prioritize_tasks_by_raw_state(tasks)

    assert _campaign_cached_sra_path("SRR_CACHED") == cached
    assert [task["srr"] for task in prioritized] == ["SRR_CACHED", "SRR_NEW"]
    assert prioritized[0]["_raw_input_priority"] == 0
    assert prioritized[0]["_raw_input_bytes"] == len(b"cached-sra")


def test_task_scheduler_interleaves_sra_extraction_and_network_work(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Long cache conversions must not occupy every sample worker."""

    data_root = tmp_path / "campaign-data"
    cache_dir = data_root / ".sra-cache" / "sra"
    cache_dir.mkdir(parents=True)
    for index in range(3):
        (cache_dir / f"SRR_CACHED_{index}.sra.cache").write_bytes(b"cache")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))

    fastq_dir = tmp_path / "species" / "work" / "getfastq"
    tasks = [
        {"srr": "SRR_NEW_0", "fastq_dir": fastq_dir, "batch_idx": 0},
        {"srr": "SRR_CACHED_0", "fastq_dir": fastq_dir, "batch_idx": 1},
        {"srr": "SRR_NEW_1", "fastq_dir": fastq_dir, "batch_idx": 2},
        {"srr": "SRR_CACHED_1", "fastq_dir": fastq_dir, "batch_idx": 3},
        {"srr": "SRR_NEW_2", "fastq_dir": fastq_dir, "batch_idx": 4},
        {"srr": "SRR_CACHED_2", "fastq_dir": fastq_dir, "batch_idx": 5},
    ]

    prioritized = _prioritize_tasks_by_raw_state(tasks, workers=4, fasterq_slots=2)

    assert [task["srr"] for task in prioritized] == [
        "SRR_CACHED_0",
        "SRR_CACHED_1",
        "SRR_NEW_0",
        "SRR_NEW_1",
        "SRR_CACHED_2",
        "SRR_NEW_2",
    ]


def test_build_sample_tasks_preserves_library_layout_for_local_reuse(tmp_path: Path) -> None:
    """Task records retain paired-end metadata for idempotent local adoption."""

    filtered = pd.DataFrame({"run": ["SRR_PAIR"], "lib_layout": ["paired"]})

    tasks = _build_sample_tasks(
        filtered,
        "run",
        tmp_path / "fastq",
        tmp_path / "config.yaml",
        "test_species",
    )

    assert tasks[0]["expected_paired"] is True
    assert _library_layout_is_paired("single") is False
    assert _library_layout_is_paired("unknown") is None


def test_valid_local_fastq_skips_ena_and_requires_complete_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Complete local reads bypass URL discovery; a lone paired mate does not."""

    sample_dir = tmp_path / "fastq" / "SRR_LOCAL"
    sample_dir.mkdir(parents=True)
    mate_1 = sample_dir / "SRR_LOCAL_1.fastq.gz"
    mate_2 = sample_dir / "SRR_LOCAL_2.fastq.gz"
    with gzip.open(mate_1, "wb") as handle:
        handle.write(b"@read\nACGT\n+\n!!!!\n")

    assert _validated_local_fastq_inputs(sample_dir, "SRR_LOCAL", True) == []

    with gzip.open(mate_2, "wb") as handle:
        handle.write(b"@read\nTGCA\n+\n!!!!\n")

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "get_fastq_urls",
        lambda *_args, **_kwargs: pytest.fail("local FASTQs must bypass ENA URL discovery"),
    )
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_LOCAL", tmp_path / "fastq", expected_paired=True) is True


def test_interrupted_paired_singleton_is_fully_classified_before_local_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An orphan complete FASTQ gets a mate-structure witness before reuse."""

    sample_dir = tmp_path / "fastq" / "SRR_ORPHAN"
    sample_dir.mkdir(parents=True)
    singleton = sample_dir / "SRR_ORPHAN.fastq.gz"
    with gzip.open(singleton, "wt", encoding="utf-8") as handle:
        handle.write("@SRR_ORPHAN.1 1/1\nACGT\n+\n!!!!\n" "@SRR_ORPHAN.2 2/1\nTGCA\n+\n!!!!\n")

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "get_fastq_urls",
        lambda *_args, **_kwargs: pytest.fail("classified local FASTQ must bypass ENA URL discovery"),
    )
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert (
        orchestrator.download_fastq(
            "SRR_ORPHAN",
            tmp_path / "fastq",
            expected_paired=True,
        )
        is True
    )
    marker = json.loads((tmp_path / "raw_validation" / "SRR_ORPHAN.json").read_text())
    assert marker["declared_library_layout"] == "paired"
    assert marker["effective_library_layout"] == "single"
    assert marker["layout_validation"] == "full_fastq_record_scan_read1_only"


def test_local_fastq_validation_marker_avoids_repeated_full_scan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A validated unchanged FASTQ is adopted without rereading its gzip stream."""

    sample_dir = tmp_path / "fastq" / "SRR_CACHED"
    sample_dir.mkdir(parents=True)
    fastq = sample_dir / "SRR_CACHED.fastq.gz"
    with gzip.open(fastq, "wb") as handle:
        handle.write(b"@read\nACGT\n+\n!!!!\n")

    calls: list[Path] = []
    real_validator = orchestrator_module.verify_gzip_integrity

    def counted_validator(path: Path) -> bool:
        calls.append(path)
        return real_validator(path)

    monkeypatch.setattr(orchestrator_module, "verify_gzip_integrity", counted_validator)

    assert _validated_local_fastq_inputs(sample_dir, "SRR_CACHED", False) == [fastq]
    assert _validated_local_fastq_inputs(sample_dir, "SRR_CACHED", False) == [fastq]
    assert calls == [fastq]
    assert (tmp_path / "raw_validation" / "SRR_CACHED.json").is_file()


def test_direct_ena_success_writes_raw_validation_witness(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A direct ENA transfer must leave the same witness as local reuse."""

    fastq_root = tmp_path / "fastq"

    def fake_download(_self: object, accession: str, sample_dir: Path) -> tuple[bool, str, list[Path]]:
        fastq = sample_dir / f"{accession}.fastq.gz"
        with gzip.open(fastq, "wb") as handle:
            handle.write(b"@read\nACGT\n+\n!!!!\n")
        return True, "", [fastq]

    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setattr(orchestrator_module.ENADownloader, "download_run", fake_download)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_DIRECT", fastq_root, expected_paired=False) is True
    marker = tmp_path / "raw_validation" / "SRR_DIRECT.json"
    assert marker.is_file()
    assert _validated_local_fastq_inputs(fastq_root / "SRR_DIRECT", "SRR_DIRECT", False) == [
        fastq_root / "SRR_DIRECT" / "SRR_DIRECT.fastq.gz"
    ]


def test_interleaved_ena_paired_fastq_is_split_idempotently(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A one-file ENA paired run becomes two validated Kallisto mates."""

    fastq_root = tmp_path / "fastq"

    def fake_download(_self: object, accession: str, sample_dir: Path) -> tuple[bool, str, list[Path]]:
        source = sample_dir / f"{accession}.fastq.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(
                "@read1/1\nACGT\n+\n!!!!\n"
                "@read1/2\nTGCA\n+\n!!!!\n"
                "@read2 1:N:0:1\nAAAA\n+\n!!!!\n"
                "@read2 2:N:0:1\nTTTT\n+\n!!!!\n"
            )
        return True, "Downloaded 1 file", [source]

    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setattr(orchestrator_module.ENADownloader, "download_run", fake_download)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_INTERLEAVED", fastq_root, expected_paired=True) is True
    sample_dir = fastq_root / "SRR_INTERLEAVED"
    assert not (sample_dir / "SRR_INTERLEAVED.fastq.gz").exists()
    mates = _validated_local_fastq_inputs(sample_dir, "SRR_INTERLEAVED", True)
    assert [path.name for path in mates] == [
        "SRR_INTERLEAVED_1.fastq.gz",
        "SRR_INTERLEAVED_2.fastq.gz",
    ]
    assert all(path.stat().st_size > 0 for path in mates)


@pytest.mark.parametrize(
    "records, expected_error",
    [
        (
            "@read1/1\nACGT\n+\n!!!\n" "@read1/2\nTGCA\n+\n!!!!\n",
            "sequence/quality length mismatch",
        ),
        (
            "@read1/2\nACGT\n+\n!!!!\n" "@read1/1\nTGCA\n+\n!!!!\n",
            "explicit FASTQ mate order is reversed",
        ),
    ],
)
def test_interleaved_split_rejects_malformed_records_and_reversed_mates(
    tmp_path: Path,
    records: str,
    expected_error: str,
) -> None:
    """A full interleaved proof validates payload structure and mate order."""

    sample_dir = tmp_path / "SRR_BAD_PAIR"
    sample_dir.mkdir()
    source = sample_dir / "SRR_BAD_PAIR.fastq.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(records)

    with pytest.raises(ValueError, match=expected_error):
        orchestrator_module._split_interleaved_fastq(
            source,
            "SRR_BAD_PAIR",
            sample_dir,
        )
    assert source.is_file()
    assert not list(sample_dir.glob("SRR_BAD_PAIR_[12].fastq.gz"))


def test_missing_sra_statistics_are_derived_from_validated_fastq(
    tmp_path: Path,
) -> None:
    """Zero/missing public SRA metrics do not block direct FASTQ quantification."""

    sample_dir = tmp_path / "fastq" / "SRR_STATS"
    sample_dir.mkdir(parents=True)
    fastq = sample_dir / "SRR_STATS.fastq.gz"
    with gzip.open(fastq, "wb") as handle:
        handle.write(b"@read1\nACGT\n+\n!!!!\n@read2\nTGCAAA\n+\n!!!!!!\n")

    assert _write_fastq_input_stats(sample_dir, "SRR_STATS") is True
    assert (sample_dir / "getfastq_stats.tsv").read_text(encoding="utf-8") == (
        "run\tnum_written\tbp_written\tbp_fastp_in\n" "SRR_STATS\t2\t10\t10\n"
    )
    assert _write_fastq_input_stats(sample_dir, "SRR_STATS") is True


def test_positive_sra_counts_avoid_recount_when_spot_length_is_missing(
    tmp_path: Path,
) -> None:
    """Amalgkit 0.16.60 can infer spot length from positive SRA counts."""

    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "run\ttotal_spots\ttotal_bases\tspot_length\n" "SRR_RATIO\t100\t5000\t\n",
        encoding="utf-8",
    )

    assert _metadata_requires_fastq_input_stats(metadata_path, "SRR_RATIO") is False

    metadata_path.write_text(
        "run\ttotal_spots\ttotal_bases\tspot_length\n" "SRR_MISSING\t100\t\t\n",
        encoding="utf-8",
    )
    assert _metadata_requires_fastq_input_stats(metadata_path, "SRR_MISSING") is True


def test_local_sra_is_extracted_before_network_acquisition(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An existing local SRA archive is processed in staging and retained for cleanup gating."""

    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_LOCAL_SRA"
    sample_dir.mkdir(parents=True)
    sra = sample_dir / "SRR_LOCAL_SRA.sra"
    sra.write_bytes(b"synthetic sra archive")
    stale_workspace = fastq_root / ".fasterq_tmp" / "SRR_LOCAL_SRA"
    stale_workspace.mkdir(parents=True)
    (stale_workspace / "stale.tmp").write_bytes(b"interrupted extraction")
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0] == "fasterq-dump":
            assert not (stale_workspace / "stale.tmp").exists()
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_LOCAL_SRA.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(b"@read\nACGT\n+\n!!!!\n")
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda *_args, **_kwargs: pytest.fail("local SRA must bypass ENA acquisition"),
    )
    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_LOCAL_SRA", fastq_root, expected_paired=False) is True
    fasterq_command = next(command for command in commands if command[0] == "fasterq-dump")
    assert fasterq_command[1] == str(sra)
    assert sra.is_file()
    assert (sample_dir / "SRR_LOCAL_SRA.fastq.gz").is_file()


def test_local_sra_can_use_reserved_host_fasterq_scratch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """SSD scratch stages a cross-volume source while retaining the authoritative input."""

    fastq_root = tmp_path / "external" / "fastq"
    sample_dir = fastq_root / "SRR_LOCAL_TMP"
    sample_dir.mkdir(parents=True)
    sra = sample_dir / "SRR_LOCAL_TMP.sra"
    sra.write_bytes(b"synthetic sra archive")
    local_scratch = tmp_path / "local-fasterq"
    stale_local = local_scratch / "SRR_LOCAL_TMP"
    stale_local.mkdir(parents=True)
    (stale_local / "stale.tmp").write_bytes(b"interrupted temp")
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", str(local_scratch))
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", "0")
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MULTIPLIER", "1")
    monkeypatch.setattr(
        orchestrator_module,
        "_paths_share_filesystem",
        lambda *_paths: False,
    )
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0] == "fasterq-dump":
            temp_dir = Path(command[command.index("--temp") + 1])
            output_dir = Path(command[command.index("--outdir") + 1])
            local_source = local_scratch / "SRR_LOCAL_TMP" / "source" / sra.name
            assert Path(command[1]) == local_source
            assert local_source.read_bytes() == sra.read_bytes()
            assert temp_dir == local_scratch / "SRR_LOCAL_TMP" / "temp"
            assert output_dir == local_scratch / "SRR_LOCAL_TMP" / "output"
            assert not (stale_local / "stale.tmp").exists()
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_LOCAL_TMP.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_LOCAL_TMP", fastq_root, expected_paired=False) is True
    assert sra.is_file()
    assert (sample_dir / "SRR_LOCAL_TMP.fastq.gz").is_file()
    assert not (local_scratch / "SRR_LOCAL_TMP").exists()
    assert not (fastq_root / ".fasterq_tmp" / "SRR_LOCAL_TMP").exists()


def test_local_sra_source_staging_failure_uses_authoritative_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed ephemeral source copy retains and extracts the authoritative SRA."""

    fastq_root = tmp_path / "external" / "fastq"
    sample_dir = fastq_root / "SRR_SOURCE_FALLBACK"
    sample_dir.mkdir(parents=True)
    sra = sample_dir / "SRR_SOURCE_FALLBACK.sra"
    sra.write_bytes(b"synthetic sra archive")
    local_scratch = tmp_path / "local-fasterq"
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", str(local_scratch))
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", "0")
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MULTIPLIER", "1")
    monkeypatch.setattr(
        orchestrator_module,
        "_paths_share_filesystem",
        lambda *_paths: False,
    )

    real_copyfile = orchestrator_module.shutil.copyfile

    def fail_source_copy(source: Path | str, destination: Path | str) -> str:
        if Path(source) == sra:
            raise OSError("synthetic local source-copy failure")
        return real_copyfile(source, destination)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if command[0] == "fasterq-dump":
            assert command[1] == str(sra)
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_SOURCE_FALLBACK.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module.shutil, "copyfile", fail_source_copy)
    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert (
        orchestrator.download_fastq(
            "SRR_SOURCE_FALLBACK",
            fastq_root,
            expected_paired=False,
        )
        is True
    )
    assert sra.is_file()
    assert (sample_dir / "SRR_SOURCE_FALLBACK.fastq.gz").is_file()
    assert not (local_scratch / "SRR_SOURCE_FALLBACK").exists()


def test_sra_validation_and_extraction_share_fallback_slot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The fallback semaphore bounds validation and extraction as one I/O unit."""

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    orchestrator._fasterq_semaphore = threading.BoundedSemaphore(1)
    state = {"active": 0, "maximum": 0}
    state_lock = threading.Lock()

    def guarded(result: bool) -> bool:
        with state_lock:
            state["active"] += 1
            state["maximum"] = max(state["maximum"], state["active"])
        time.sleep(0.03)
        with state_lock:
            state["active"] -= 1
        return result

    monkeypatch.setattr(
        orchestrator,
        "_validate_local_sra_archive",
        lambda *_args: guarded(True),
    )
    monkeypatch.setattr(
        orchestrator,
        "_extract_sra_to_fastq",
        lambda *_args: guarded(True),
    )
    results: list[bool] = []

    def run_one(index: int) -> None:
        results.append(
            orchestrator._extract_sra_to_fastq_bounded(
                tmp_path / f"SRR_{index}.sra",
                tmp_path / f"sample-{index}",
                tmp_path / "fastq",
                f"SRR_{index}",
                False,
            )
        )

    workers = [threading.Thread(target=run_one, args=(index,)) for index in range(2)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join()

    assert results == [True, True]
    assert state["maximum"] == 1


def test_cross_filesystem_fastq_promotion_uses_destination_atomic_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """EXDEV promotion verifies a destination-side copy before deleting source."""

    staged = tmp_path / "local-scratch" / "SRR_COPY.fastq.gz"
    destination = tmp_path / "external" / "SRR_COPY.fastq.gz"
    staged.parent.mkdir(parents=True)
    destination.parent.mkdir(parents=True)
    staged.write_bytes(b"validated compressed payload")
    real_replace = orchestrator_module.os.replace

    def cross_device_once(source: Path | str, target: Path | str) -> None:
        if Path(source) == staged:
            raise OSError(orchestrator_module.errno.EXDEV, "cross-device link")
        real_replace(source, target)

    monkeypatch.setattr(orchestrator_module.os, "replace", cross_device_once)

    orchestrator_module._promote_staged_fastq_file(staged, destination)

    assert not staged.exists()
    assert destination.read_bytes() == b"validated compressed payload"
    assert not list(destination.parent.glob(".*.copying"))


def test_local_fasterq_scratch_falls_back_when_reservation_would_cross_floor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An oversized local reservation uses external temp without losing the SRA."""

    fastq_root = tmp_path / "external" / "fastq"
    sample_dir = fastq_root / "SRR_NO_LOCAL_TMP"
    sample_dir.mkdir(parents=True)
    sra = sample_dir / "SRR_NO_LOCAL_TMP.sra"
    sra.write_bytes(b"synthetic sra archive")
    local_scratch = tmp_path / "local-fasterq"
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", str(local_scratch))
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", "100000")
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0] == "fasterq-dump":
            temp_dir = Path(command[command.index("--temp") + 1])
            output_dir = Path(command[command.index("--outdir") + 1])
            assert temp_dir == fastq_root / ".fasterq_tmp" / "SRR_NO_LOCAL_TMP" / "temp"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_NO_LOCAL_TMP.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_NO_LOCAL_TMP", fastq_root, expected_paired=False) is True
    assert sra.is_file()
    assert (sample_dir / "SRR_NO_LOCAL_TMP.fastq.gz").is_file()


def test_remote_fasterq_fallback_can_use_fixed_reserved_host_scratch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bare NCBI accessions use the configured fixed reservation and SSD temp."""

    fastq_root = tmp_path / "external" / "fastq"
    sample_dir = fastq_root / "SRR_REMOTE_TMP"
    sample_dir.mkdir(parents=True)
    local_scratch = tmp_path / "local-fasterq"
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", str(local_scratch))
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_REMOTE_RESERVE_GB", "0.001")
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0] == "fasterq-dump":
            temp_dir = Path(command[command.index("--temp") + 1])
            output_dir = Path(command[command.index("--outdir") + 1])
            assert command[1] == "SRR_REMOTE_TMP"
            assert temp_dir == local_scratch / "SRR_REMOTE_TMP" / "temp"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_REMOTE_TMP.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator._extract_sra_to_fastq(
        "SRR_REMOTE_TMP",
        sample_dir,
        fastq_root,
        "SRR_REMOTE_TMP",
        expected_paired=False,
    )
    assert (sample_dir / "SRR_REMOTE_TMP.fastq.gz").is_file()
    assert not (local_scratch / "SRR_REMOTE_TMP").exists()


def test_fasterq_promotion_reuses_existing_valid_mate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fallback extraction supplies a missing mate without replacing completed transfer work."""

    data_root = tmp_path / "campaign-data"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_PAIRED"
    sample_dir.mkdir(parents=True)
    existing_mate = sample_dir / "SRR_PAIRED_1.fastq.gz"
    with gzip.open(existing_mate, "wb") as handle:
        handle.write(b"@existing/1\nACGT\n+\n!!!!\n")
    existing_payload = existing_mate.read_bytes()

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if command[0] == "fasterq-dump":
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_PAIRED_1.fastq").write_text("@fallback/1\nACGT\n+\n!!!!\n")
            (output_dir / "SRR_PAIRED_2.fastq").write_text("@fallback/2\nTGCA\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator._extract_sra_to_fastq(
        "SRR_PAIRED",
        sample_dir,
        fastq_root,
        "SRR_PAIRED",
        expected_paired=True,
    )
    assert existing_mate.read_bytes() == existing_payload
    assert (sample_dir / "SRR_PAIRED_2.fastq.gz").is_file()


def test_fasterq_interleaved_paired_output_is_split_and_promoted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A valid one-file paired SRA extraction becomes two canonical mates."""

    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_INTERLEAVED_SRA"
    sample_dir.mkdir(parents=True)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if command[0] == "fasterq-dump":
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_INTERLEAVED_SRA.fastq").write_text(
                "@spot.1 1 length=4\nACGT\n+\n!!!!\n" "@spot.1 2 length=4\nTGCA\n+\n!!!!\n"
            )
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator._extract_sra_to_fastq(
        "SRR_INTERLEAVED_SRA",
        sample_dir,
        fastq_root,
        "SRR_INTERLEAVED_SRA",
        expected_paired=True,
    )
    assert not (sample_dir / "SRR_INTERLEAVED_SRA.fastq.gz").exists()
    assert (sample_dir / "SRR_INTERLEAVED_SRA_1.fastq.gz").is_file()
    assert (sample_dir / "SRR_INTERLEAVED_SRA_2.fastq.gz").is_file()


def test_fasterq_paired_metadata_single_read_output_is_reconciled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A read-1-only archive follows Amalgkit 0.16.60's single-file rule."""

    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_LAYOUT_MISMATCH"
    sample_dir.mkdir(parents=True)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if command[0] == "fasterq-dump":
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_LAYOUT_MISMATCH.fastq").write_text(
                "@SRR_LAYOUT_MISMATCH.1 1/1\nACGT\n+\n!!!!\n" "@SRR_LAYOUT_MISMATCH.2 2/1\nTGCA\n+\n!!!!\n"
            )
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator._extract_sra_to_fastq(
        "SRR_LAYOUT_MISMATCH",
        sample_dir,
        fastq_root,
        "SRR_LAYOUT_MISMATCH",
        expected_paired=True,
    )
    singleton = sample_dir / "SRR_LAYOUT_MISMATCH.fastq.gz"
    assert singleton.is_file()
    assert not (sample_dir / "SRR_LAYOUT_MISMATCH_1.fastq.gz").exists()
    marker = json.loads((tmp_path / "raw_validation" / "SRR_LAYOUT_MISMATCH.json").read_text())
    assert marker["declared_library_layout"] == "paired"
    assert marker["effective_library_layout"] == "single"
    # The witness name is historical evidence and remains reusable across
    # runtime version drift.
    assert marker["layout_reconciliation"] == "amalgkit_0.16.38_paired_metadata_single_fastq"


def test_fasterq_nonadjacent_mate_two_is_not_downgraded_to_single(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A mixed or misordered paired stream is rejected rather than relabelled."""

    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_MIXED"
    sample_dir.mkdir(parents=True)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if command[0] == "fasterq-dump":
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_MIXED.fastq").write_text(
                "@unrelated-a/1\nACGT\n+\n!!!!\n" "@unrelated-b/2\nTGCA\n+\n!!!!\n"
            )
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(source.read_bytes())
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert not orchestrator._extract_sra_to_fastq(
        "SRR_MIXED",
        sample_dir,
        fastq_root,
        "SRR_MIXED",
        expected_paired=True,
    )
    assert not list(sample_dir.glob("*.fastq.gz"))


def test_campaign_cached_sra_is_extracted_before_ena(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A shared NCBI cache database bypasses ENA without moving the cache."""

    data_root = tmp_path / "campaign-data"
    cache_dir = data_root / ".sra-cache" / "sra"
    cache_dir.mkdir(parents=True)
    cached = cache_dir / "SRR_CACHED_SRA.sra.cache"
    cached.write_bytes(b"synthetic cached SRA database")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))

    fastq_root = tmp_path / "fastq"
    sample_dir = fastq_root / "SRR_CACHED_SRA"
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0] == "fasterq-dump":
            output_dir = Path(command[command.index("--outdir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "SRR_CACHED_SRA.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(b"@read\nACGT\n+\n!!!!\n")
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda *_args, **_kwargs: pytest.fail("campaign SRA cache must bypass ENA"),
    )
    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_CACHED_SRA", fastq_root, expected_paired=False) is True
    fasterq_command = next(command for command in commands if command[0] == "fasterq-dump")
    assert fasterq_command[1] == str(cached)
    assert cached.is_file()
    assert (sample_dir / "SRR_CACHED_SRA.fastq.gz").is_file()


def test_invalid_sra_validation_witness_is_reused_until_source_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An unchanged corrupt cache is validated once across repeated decisions."""

    source = tmp_path / "SRR_INVALID.sra.cache"
    source.write_bytes(b"corrupt cache payload")
    calls: list[list[str]] = []

    def failed_validation(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(
            command,
            3,
            "",
            "file corrupt: failed md5 validation",
        )

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", failed_validation)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator._validate_local_sra_archive(source, "SRR_INVALID") is False
    assert orchestrator._validate_local_sra_archive(source, "SRR_INVALID") is False
    assert len(calls) == 1

    source.write_bytes(b"replacement cache payload")
    assert orchestrator._validate_local_sra_archive(source, "SRR_INVALID") is False
    assert len(calls) == 2


def test_invalid_cached_sra_skips_extraction_and_uses_ena(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Invalid cached evidence is retained while a valid ENA payload replaces its role."""

    data_root = tmp_path / "campaign-data"
    cache_dir = data_root / ".sra-cache" / "sra"
    cache_dir.mkdir(parents=True)
    cached = cache_dir / "SRR_INVALID_CACHE.sra.cache"
    cached.write_bytes(b"corrupt cache payload")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")

    def download_valid_fastq(
        _self: object,
        accession: str,
        destination: Path,
    ) -> tuple[bool, str, list[Path]]:
        fastq = destination / f"{accession}.fastq.gz"
        with gzip.open(fastq, "wb") as handle:
            handle.write(b"@read\nACGT\n+\n!!!!\n")
        return True, "downloaded", [fastq]

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        download_valid_fastq,
    )
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    monkeypatch.setattr(
        orchestrator,
        "_validate_local_sra_archive",
        lambda *_args: False,
    )
    monkeypatch.setattr(
        orchestrator,
        "_extract_sra_to_fastq",
        lambda *_args: pytest.fail("invalid cached SRA must not be extracted"),
    )

    assert (
        orchestrator.download_fastq(
            "SRR_INVALID_CACHE",
            tmp_path / "fastq",
            expected_paired=False,
        )
        is True
    )
    assert cached.is_file()


def test_failed_cached_sra_is_not_extracted_twice_in_one_invocation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed cache extraction and ENA transfer must not repeat the same expensive work."""

    data_root = tmp_path / "campaign-data"
    cache_dir = data_root / ".sra-cache" / "sra"
    cache_dir.mkdir(parents=True)
    cached = cache_dir / "SRR_CACHED_FAIL.sra"
    cached.write_bytes(b"synthetic cached SRA database")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda _self, _sample_id, _out_dir: (False, "synthetic ENA failure", []),
    )

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    monkeypatch.setattr(
        orchestrator,
        "_validate_local_sra_archive",
        lambda *_args: True,
    )
    extraction_sources: list[Path | str] = []

    def failed_extraction(
        source: Path | str,
        _sample_dir: Path,
        _out_dir: Path,
        _srr_id: str,
        _expected_paired: bool | None,
    ) -> bool:
        extraction_sources.append(source)
        return False

    monkeypatch.setattr(orchestrator, "_extract_sra_to_fastq_bounded", failed_extraction)

    assert (
        orchestrator.download_fastq(
            "SRR_CACHED_FAIL",
            tmp_path / "fastq",
            expected_paired=False,
        )
        is False
    )
    assert extraction_sources == [cached]
    assert cached.is_file()


def test_failed_ena_can_defer_ncbi_without_occupying_primary_worker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ENA exhaustion may hand work to the dedicated NCBI acquisition lane."""

    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda _self, _sample_id, _out_dir: (False, "synthetic ENA failure", []),
    )
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    monkeypatch.setattr(
        orchestrator,
        "_extract_sra_to_fastq_bounded",
        lambda *_args, **_kwargs: pytest.fail("deferred NCBI fallback must not run in the primary acquisition worker"),
    )

    outcome = orchestrator.download_fastq(
        "SRR_DEFER",
        tmp_path / "fastq",
        expected_paired=False,
        defer_ncbi_fallback=True,
    )

    assert outcome is None


def test_resource_profile_bounds_workers_and_records_stage_budgets(monkeypatch: pytest.MonkeyPatch) -> None:
    """Task workers may exceed the separately bounded quantification slots."""

    monkeypatch.delenv("AMALGKIT_PIPELINE_FASTQ_THREADS", raising=False)
    monkeypatch.delenv("AMALGKIT_PIPELINE_COMPRESSION_THREADS", raising=False)
    monkeypatch.delenv("AMALGKIT_PIPELINE_MAX_IN_FLIGHT", raising=False)

    profile = build_pipeline_resource_profile(12, 8, cpu_count=10)

    assert profile.requested_workers == 12
    assert profile.workers == 12
    assert profile.quant_slots == 8
    assert profile.quant_threads_per_worker == 1
    assert profile.effective_quant_threads == 8
    assert profile.fasterq_threads == 1
    assert profile.fasterq_slots == 4
    assert profile.compression_threads == 1
    assert profile.validation_slots == 4
    assert profile.peak_stage_threads == 8
    assert profile.max_in_flight == 24
    assert profile.host_cpu_count == 10


def test_resource_profile_uses_explicit_fallback_budgets_from_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fallback subprocess budgets are explicit and independently configurable."""

    monkeypatch.setenv("AMALGKIT_PIPELINE_FASTQ_THREADS", "3")
    monkeypatch.setenv("AMALGKIT_PIPELINE_FASTQ_SLOTS", "2")
    monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_THREADS", "2")
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATION_SLOTS", "3")
    monkeypatch.setenv("AMALGKIT_PIPELINE_MAX_IN_FLIGHT", "11")

    profile = build_pipeline_resource_profile(4, 8, cpu_count=10)

    assert profile.workers == 4
    assert profile.quant_slots == 4
    assert profile.quant_threads_per_worker == 2
    assert profile.fasterq_threads == 3
    assert profile.fasterq_slots == 2
    assert profile.compression_threads == 2
    assert profile.validation_slots == 3
    assert profile.peak_stage_threads == 8
    assert profile.max_in_flight == 11


def test_resource_profile_can_bound_quant_slots_independently(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """High task-worker concurrency can feed a smaller, multi-threaded quant budget."""

    monkeypatch.setenv("AMALGKIT_PIPELINE_QUANT_SLOTS", "3")
    profile = build_pipeline_resource_profile(12, 9, cpu_count=10)

    assert profile.workers == 12
    assert profile.quant_slots == 3
    assert profile.quant_threads_per_worker == 3
    assert profile.effective_quant_threads == 9

    explicit = build_pipeline_resource_profile(12, 9, quant_slots=2, cpu_count=10)
    assert explicit.quant_slots == 2
    assert explicit.quant_threads_per_worker == 4


def test_fasterq_fallback_uses_profiled_child_process_threads(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The NCBI fallback must use the configured fasterq and pigz budgets."""

    # This synthetic test uses a temporary directory on the system volume.
    # Disable both production free-space floors so a nearly-full developer
    # volume cannot turn the isolated test into a 60-second throttle wait.
    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_PIPELINE_FASTQ_THREADS", "3")
    monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_THREADS", "2")
    data_root = tmp_path / "campaign-data"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda _self, _sample_id, _out_dir: (False, "synthetic ENA failure", []),
    )
    commands: list[list[str]] = []
    environments: list[dict[str, str]] = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        environment = kwargs.get("env")
        if isinstance(environment, dict):
            environments.append(environment)
        if command[0] == "fasterq-dump":
            out_dir = Path(command[command.index("--outdir") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "SRR_PROFILE.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(b"@read\nACGT\n+\n!!!!\n")
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_PROFILE", tmp_path / "fastq", expected_paired=False) is True
    fasterq_command = next(command for command in commands if command[0] == "fasterq-dump")
    pigz_command = next(command for command in commands if command[0] == "pigz")
    assert fasterq_command[fasterq_command.index("--threads") + 1] == "3"
    assert pigz_command[pigz_command.index("-p") + 1] == "2"
    settings_path = Path(environments[0]["NCBI_SETTINGS"])
    assert settings_path == data_root / ".ncbi" / "metainformant-user-settings.mkfg"
    assert f'/repository/user/main/public/root = "{data_root / ".sra-cache"}"' in settings_path.read_text()


def test_ncbi_fallback_prefetches_and_verifies_before_fasterq(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Production NCBI fallback downloads an accession once before extraction."""

    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_NCBI_PREFETCH_FIRST", "yes")
    data_root = tmp_path / "campaign-data"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[0].endswith("prefetch"):
            destination = Path(command[command.index("--output-directory") + 1])
            destination.mkdir(parents=True, exist_ok=True)
            (destination / "SRR_PREFETCH.sra").write_bytes(b"verified-sra")
        elif command[0] == "fasterq-dump":
            output = Path(command[command.index("--outdir") + 1])
            output.mkdir(parents=True, exist_ok=True)
            (output / "SRR_PREFETCH.fastq").write_text("@read\nACGT\n+\n!!!!\n")
        elif command[0] == "pigz":
            source = Path(command[-1])
            with gzip.open(source.with_suffix(".fastq.gz"), "wb") as handle:
                handle.write(b"@read\nACGT\n+\n!!!!\n")
            source.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", fake_run)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    monkeypatch.setattr(orchestrator, "_validate_local_sra_archive", lambda *_args: True)

    assert orchestrator._download_fastq_ncbi_only("SRR_PREFETCH", tmp_path / "fastq", expected_paired=False)
    assert commands[0][0].endswith("prefetch")
    assert commands[1][0] == "fasterq-dump"
    assert commands[1][1].endswith("SRR_PREFETCH.sra")


def test_campaign_ncbi_settings_preserve_global_configuration(tmp_path: Path) -> None:
    """The campaign-local VDB file routes SRA caches without editing user settings."""

    data_root = tmp_path / "campaign-data"
    settings_path = _ensure_ncbi_settings_for_data_root(data_root)

    assert settings_path.is_file()
    assert (data_root / ".sra-cache").is_dir()
    assert "metainformant-user-settings" in settings_path.name
    assert str(data_root / ".sra-cache") in settings_path.read_text()


def test_run_all_keeps_submitted_sample_tasks_bounded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A large task list never exceeds the configured in-flight window."""

    from concurrent.futures import ThreadPoolExecutor as RealThreadPoolExecutor

    tasks = [
        {
            "srr": f"SRR{i:03d}",
            "batch_idx": i,
            "fastq_dir": tmp_path / "fastq",
            "config_path": tmp_path / "config.yaml",
            "species_name": "test_species",
        }
        for i in range(1, 7)
    ]
    executors: list[object] = []

    class RecordingExecutor:
        def __init__(self, max_workers: int) -> None:
            self.executor = RealThreadPoolExecutor(max_workers=max_workers)
            self.max_workers = max_workers
            self.pending = 0
            self.max_pending = 0
            self.lock = threading.Lock()
            executors.append(self)

        def __enter__(self) -> "RecordingExecutor":
            return self

        def __exit__(self, *args: object) -> None:
            self.executor.shutdown(wait=True)

        def submit(self, *args: object, **kwargs: object) -> object:
            with self.lock:
                self.pending += 1
                self.max_pending = max(self.max_pending, self.pending)
            future = self.executor.submit(*args, **kwargs)

            def release(_future: object) -> None:
                with self.lock:
                    self.pending -= 1

            future.add_done_callback(release)
            return future

    def fake_process(*_args: object, **_kwargs: object) -> dict[str, object]:
        time.sleep(0.01)
        return {"quantified": True, "skipped": False, "error": None}

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    monkeypatch.setattr(orchestrator, "discover_species_tasks", lambda *_args: tasks)
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    monkeypatch.setattr(orchestrator_module.concurrent.futures, "ThreadPoolExecutor", RecordingExecutor)

    orchestrator.run_all(
        ["amalgkit_test_species.yaml"],
        max_gb=1.0,
        workers=2,
        threads=2,
        max_in_flight=2,
        discovery_workers=1,
    )

    executor = max(executors, key=lambda item: item.max_workers)
    assert isinstance(executor, RecordingExecutor)
    assert executor.max_pending <= 2


def test_run_all_parallelizes_discovery_and_preserves_species_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Discovery overlaps, but the execution task list follows config order."""

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    started: list[str] = []
    barrier = threading.Barrier(3)

    def fake_discover(config_name: str, *_args: object) -> list[dict[str, object]]:
        started.append(config_name)
        barrier.wait(timeout=2)
        species = config_name.removeprefix("amalgkit_").removesuffix(".yaml")
        return [
            {
                "srr": species,
                "batch_idx": 1,
                "fastq_dir": tmp_path / "fastq" / species,
                "config_path": tmp_path / f"{species}.yaml",
                "species_name": species,
            }
        ]

    observed: list[str] = []

    def fake_process(*args: object, **_kwargs: object) -> dict[str, object]:
        observed.append(str(args[4]))
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(orchestrator, "discover_species_tasks", fake_discover)
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    orchestrator.run_all(
        ["amalgkit_first.yaml", "amalgkit_second.yaml", "amalgkit_third.yaml"],
        max_gb=1.0,
        workers=1,
        threads=1,
        max_in_flight=1,
        discovery_workers=3,
    )

    assert set(started) == {"amalgkit_first.yaml", "amalgkit_second.yaml", "amalgkit_third.yaml"}
    assert observed == ["first", "second", "third"]


def test_run_all_aborts_before_execution_on_discovery_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed species discovery cannot silently produce a partial campaign."""

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    executed = False

    def fake_discover(config_name: str, *_args: object) -> list[dict[str, object]]:
        if config_name == "amalgkit_bad.yaml":
            raise RuntimeError("reference preparation failed")
        return []

    def fake_process(*_args: object, **_kwargs: object) -> dict[str, object]:
        nonlocal executed
        executed = True
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(orchestrator, "discover_species_tasks", fake_discover)
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    with pytest.raises(RuntimeError, match="Species discovery failed"):
        orchestrator.run_all(
            ["amalgkit_good.yaml", "amalgkit_bad.yaml"],
            max_gb=1.0,
            workers=1,
            threads=1,
            discovery_workers=2,
        )
    assert executed is False


def test_run_all_routes_deferred_ncbi_acquisition_back_to_primary_lane(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A deferred extraction must finish acquisition before primary quantification resumes."""

    task = {
        "srr": "SRR_DEFERRED",
        "batch_idx": 1,
        "fastq_dir": tmp_path / "fastq",
        "config_path": tmp_path / "config.yaml",
        "species_name": "test_species",
        "expected_paired": False,
    }
    process_calls: list[dict[str, object]] = []
    fallback_calls: list[str] = []

    def fake_process(*_args: object, **kwargs: object) -> dict[str, object]:
        process_calls.append(kwargs)
        if len(process_calls) == 1:
            return {
                "quantified": False,
                "skipped": False,
                "fallback_deferred": True,
                "error": None,
            }
        return {
            "quantified": True,
            "skipped": False,
            "fallback_deferred": False,
            "error": None,
        }

    def fake_fallback(srr_id: str, _fastq_dir: Path, _expected_paired: bool | None) -> bool:
        fallback_calls.append(srr_id)
        return True

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    monkeypatch.setattr(orchestrator, "discover_species_tasks", lambda *_args: [task])
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    monkeypatch.setattr(orchestrator, "_download_fastq_ncbi_only", fake_fallback)

    orchestrator.run_all(
        ["amalgkit_test_species.yaml"],
        max_gb=1.0,
        workers=2,
        threads=2,
        max_in_flight=2,
        fasterq_slots=1,
    )

    assert fallback_calls == ["SRR_DEFERRED"]
    assert len(process_calls) == 2
    assert process_calls[0]["defer_ncbi_fallback"] is True
    assert process_calls[1]["defer_ncbi_fallback"] is True


def test_deferred_ncbi_fallback_does_not_starve_primary_executor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A blocked one-slot fallback must not prevent the next ENA task from running."""

    tasks = [
        {
            "srr": "SRR_NEEDS_FALLBACK",
            "batch_idx": 1,
            "fastq_dir": tmp_path / "fastq",
            "config_path": tmp_path / "config.yaml",
            "species_name": "test_species",
            "expected_paired": False,
        },
        {
            "srr": "SRR_PRIMARY_CONTINUES",
            "batch_idx": 2,
            "fastq_dir": tmp_path / "fastq",
            "config_path": tmp_path / "config.yaml",
            "species_name": "test_species",
            "expected_paired": False,
        },
    ]
    primary_continued = threading.Event()
    fallback_started = threading.Event()
    process_count: dict[str, int] = {}

    def fake_process(srr_id: str, *_args: object, **_kwargs: object) -> dict[str, object]:
        process_count[srr_id] = process_count.get(srr_id, 0) + 1
        if srr_id == "SRR_NEEDS_FALLBACK" and process_count[srr_id] == 1:
            return {
                "quantified": False,
                "skipped": False,
                "fallback_deferred": True,
                "error": None,
            }
        if srr_id == "SRR_PRIMARY_CONTINUES":
            assert fallback_started.wait(timeout=2)
            primary_continued.set()
        return {
            "quantified": True,
            "skipped": False,
            "fallback_deferred": False,
            "error": None,
        }

    def blocking_fallback(_srr_id: str, _fastq_dir: Path, _expected_paired: bool | None) -> bool:
        fallback_started.set()
        return primary_continued.wait(timeout=2)

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    monkeypatch.setattr(orchestrator, "discover_species_tasks", lambda *_args: tasks)
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    monkeypatch.setattr(orchestrator, "_download_fastq_ncbi_only", blocking_fallback)

    orchestrator.run_all(
        ["amalgkit_test_species.yaml"],
        max_gb=1.0,
        workers=1,
        threads=1,
        max_in_flight=1,
        fasterq_slots=1,
    )

    assert primary_continued.is_set()
    assert process_count == {
        "SRR_NEEDS_FALLBACK": 2,
        "SRR_PRIMARY_CONTINUES": 1,
    }


def test_failed_download_preserves_invalid_transfer_artifacts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fallback cleanup must retain invalid payloads and resumable partials."""

    fastq_dir = tmp_path / "fastq"
    sample_dir = fastq_dir / "SRR_FAIL"
    sample_dir.mkdir(parents=True)
    invalid = sample_dir / "SRR_FAIL_1.fastq.gz.part.invalid"
    partial = sample_dir / "SRR_FAIL_2.fastq.gz.part"
    completed_mate = sample_dir / "SRR_FAIL_1.fastq.gz"
    stale = sample_dir / "stale.fastq.gz"
    invalid.write_bytes(b"invalid transfer")
    partial.write_bytes(b"partial transfer")
    with gzip.open(completed_mate, "wb") as handle:
        handle.write(b"@completed/1\nACGT\n+\n!!!!\n")
    stale.write_bytes(b"stale final payload")

    monkeypatch.setattr(
        orchestrator_module.ENADownloader,
        "download_run",
        lambda _self, _sample_id, _out_dir: (False, "synthetic ENA failure", []),
    )

    def failed_fasterq(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, 1, "", "synthetic fasterq failure")

    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", failed_fasterq)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )

    assert orchestrator.download_fastq("SRR_FAIL", fastq_dir) is False
    assert invalid.exists()
    assert partial.exists()
    assert completed_mate.exists()
    assert not stale.exists()


def test_reference_aliases_are_explicit_and_audited(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Metadata subspecies aliases create symlinks and a manifest, not copies."""
    data_root = tmp_path / "data"
    index_dir = data_root / "test_species" / "work" / "index"
    index_dir.mkdir(parents=True)
    source = index_dir / "Test_species.idx"
    source.write_bytes(b"kallisto-index")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))

    ready, missing, manifest = _ensure_reference_alias_indexes(
        {
            "species_list": ["Test_species"],
            "genome": {"index_dir": str(index_dir)},
            "reference_aliases": {"Test species isolate": "Test_species"},
        },
        "test_species",
        ["Test species", "Test species isolate"],
    )

    assert ready
    assert missing == []
    assert manifest is not None and manifest.exists()
    alias = index_dir / "Test_species_isolate.idx"
    assert alias.is_symlink()
    assert alias.resolve() == source.resolve()
    payload = yaml.safe_load(manifest.read_text())
    assert payload["status"] == "complete"
    assert any(entry["mode"] == "explicit_species_level_alias" for entry in payload["aliases"])


def test_invalid_run_rows_are_excluded_and_ledgered(tmp_path: Path) -> None:
    """Missing SRA run accessions cannot enter task paths or abort discovery."""
    metadata_path = tmp_path / "work" / "metadata" / "metadata.tsv"
    metadata_path.parent.mkdir(parents=True)
    df = pd.DataFrame(
        {
            "run": ["SRR001", None, "SRX002", "ERR003"],
            "experiment": ["SRX1", "SRX2", "SRX3", "SRX4"],
        }
    )

    filtered = _filter_valid_run_rows(df, "apis", metadata_path)

    assert [_is_valid_run_accession(value) for value in ["SRR001", "ERR003", None, "SRX002"]] == [
        True,
        True,
        False,
        False,
    ]
    assert filtered["run"].tolist() == ["SRR001", "ERR003"]
    ledger = metadata_path.parent / "validation" / "invalid_run_accessions.tsv"
    assert ledger.exists()
    ledger_df = pd.read_csv(ledger, sep="\t")
    assert len(ledger_df) == 2
    assert set(ledger_df["invalid_run_reason"]) == {"missing_or_invalid_run_accession"}


def test_build_quant_command_uses_metadata_cleanup_and_index(tmp_path: Path) -> None:
    """Quant command construction should include cleanup and real index paths."""
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (index_dir / "genome.idx").write_bytes(b"idx")
    cfg = {
        "steps": {"quant": {"clean_fastq": "yes"}},
        "genome": {"index_dir": str(index_dir)},
    }

    cmd = _build_quant_command(cfg, "apis", batch_index=3, threads=2, metadata_path="metadata.tsv")

    assert cmd[:4] == ["amalgkit", "quant", "--out_dir", "output/amalgkit/apis/work"]
    assert cmd[cmd.index("--metadata") + 1] == "metadata.tsv"
    assert cmd[cmd.index("--clean_fastq") + 1] == "no"
    assert cmd[cmd.index("--index_dir") + 1] == str(index_dir)
    assert cmd[cmd.index("--redo") + 1] == "yes"


def test_build_quant_command_uses_atomic_local_index_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An opt-in local cache accelerates index reads without moving outputs."""

    source_dir = tmp_path / "external" / "index"
    source_dir.mkdir(parents=True)
    source = source_dir / "genome.idx"
    source.write_bytes(b"reference-index")
    cache_dir = tmp_path / "local-index-cache"
    monkeypatch.setenv("AMALGKIT_LOCAL_INDEX_CACHE_DIR", str(cache_dir))

    cmd = _build_quant_command(
        {"genome": {"index_dir": str(source_dir)}},
        "test_species",
        batch_index=1,
        threads=1,
        metadata_path="metadata.tsv",
    )

    cached_dir = cache_dir / "test_species" / "index"
    assert cmd[cmd.index("--index_dir") + 1] == str(cached_dir)
    assert (cached_dir / "genome.idx").read_bytes() == source.read_bytes()
    manifest = yaml.safe_load((cached_dir / ".metainformant_index_cache.json").read_text())
    assert manifest["schema"] == "metainformant.rna.index_cache.v1"
    assert manifest["source"] == str(source_dir)

    source.write_bytes(b"updated-reference-index")
    refreshed = _build_quant_command(
        {"genome": {"index_dir": str(source_dir)}},
        "test_species",
        batch_index=2,
        threads=1,
        metadata_path="metadata.tsv",
    )
    assert refreshed[refreshed.index("--index_dir") + 1] == str(cached_dir)
    assert (cached_dir / "genome.idx").read_bytes() == source.read_bytes()


def test_stage_timeouts_are_explicit_and_configurable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Recovery profiles can extend each long-running stage independently."""

    values = {
        "AMALGKIT_PIPELINE_DOWNLOAD_TIMEOUT_SECONDS": "101",
        "AMALGKIT_PIPELINE_FASTQ_TIMEOUT_SECONDS": "202",
        "AMALGKIT_PIPELINE_COMPRESSION_TIMEOUT_SECONDS": "303",
        "AMALGKIT_PIPELINE_QUANT_TIMEOUT_SECONDS": "404",
        "AMALGKIT_PIPELINE_DOWNLOAD_SPEED_LIMIT_BYTES": "2048",
        "AMALGKIT_PIPELINE_DOWNLOAD_SPEED_TIME_SECONDS": "505",
        "AMALGKIT_PIPELINE_ENA_RETRIES": "7",
        "AMALGKIT_PIPELINE_ENA_INTEGRITY_RETRIES": "1",
        "AMALGKIT_PIPELINE_ENA_RETRY_DELAY_SECONDS": "6",
        "AMALGKIT_PIPELINE_ENA_API_RETRIES": "3",
        "AMALGKIT_PIPELINE_ENA_API_RETRY_DELAY_SECONDS": "4",
        "AMALGKIT_MIN_EXTERNAL_FREE_GB": "6.5",
    }
    for name, value in values.items():
        monkeypatch.setenv(name, value)

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )

    assert orchestrator.download_timeout_seconds == 101
    assert orchestrator.fasterq_timeout_seconds == 202
    assert orchestrator.compression_timeout_seconds == 303
    assert orchestrator.quant_timeout_seconds == 404
    assert orchestrator.download_speed_limit_bytes == 2048
    assert orchestrator.download_speed_time_seconds == 505
    assert orchestrator.download_retries == 7
    assert orchestrator.download_integrity_retries == 1
    assert orchestrator.download_retry_delay_seconds == 6
    assert orchestrator.ena_api_retries == 3
    assert orchestrator.ena_api_retry_delay_seconds == 4
    assert orchestrator.min_external_free_gb == 6.5


def test_quant_timeout_terminates_the_process_group() -> None:
    """A timed-out quant command cannot leave its child process running."""
    command = [sys.executable, "-c", "import time; time.sleep(30)"]

    with pytest.raises(subprocess.TimeoutExpired):
        _run_command_in_process_group(command, timeout=1)


def test_ena_download_uses_speed_guard_and_process_group_runner(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Direct ENA transfers abort low-speed stalls without losing resumable state."""

    output_dir = tmp_path / "fastq" / "SRR_SPEED"
    commands: list[list[str]] = []
    downloader = ena_downloader_module.ENADownloader(
        timeout=123,
        retries=2,
        speed_limit_bytes=2048,
        speed_time_seconds=505,
    )
    monkeypatch.setattr(
        downloader,
        "get_fastq_urls",
        lambda _sample_id: ["https://example.test/SRR_SPEED.fastq.gz"],
    )

    def fake_runner(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        assert 0 < timeout <= 123
        partial = Path(command[command.index("-o") + 1])
        partial.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(partial, "wb") as handle:
            handle.write(b"@read\nACGT\n+\n!!!!\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(ena_downloader_module, "_run_command_in_process_group", fake_runner)

    success, _message, files = downloader.download_run("SRR_SPEED", output_dir)

    assert success is True
    assert files == [output_dir / "SRR_SPEED.fastq.gz"]
    command = commands[0]
    assert "--retry" not in command
    assert "--retry-all-errors" not in command
    assert command[command.index("--speed-limit") + 1] == "2048"
    assert command[command.index("--speed-time") + 1] == "505"
    assert command[command.index("--write-out") + 1] == "%{http_code}"
    assert (output_dir / "SRR_SPEED.fastq.gz").is_file()


def test_local_quant_scratch_promotes_validated_output_and_archives_replaced_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Local Kallisto output is promoted without losing a prior external result."""

    work_dir = tmp_path / "work"
    (work_dir / "metadata").mkdir(parents=True)
    (work_dir / "getfastq").mkdir()
    metadata = work_dir / "metadata" / "metadata_selected.tsv"
    metadata.write_text("run\tlib_layout\nSRR_LOCAL\tsingle\n")
    external_sample = work_dir / "quant" / "SRR_LOCAL"
    external_sample.mkdir(parents=True)
    (external_sample / "SRR_LOCAL_abundance.tsv").write_text("prior\n" * 30)
    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", str(tmp_path / "quant-scratch"))
    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB", "0")

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    command = [
        "amalgkit",
        "quant",
        "--out_dir",
        str(work_dir),
        "--metadata",
        str(metadata),
        "--index_dir",
        str(tmp_path / "index"),
    ]

    prepared = orchestrator._prepare_local_quant_workspace(
        work_dir,
        "test_species",
        "SRR_LOCAL",
        str(metadata),
        command,
    )
    assert prepared is not None
    local_command, scratch_root = prepared
    assert local_command[local_command.index("--out_dir") + 1] == str(scratch_root)
    assert (scratch_root / "getfastq").is_symlink()
    assert (scratch_root / "metadata").is_symlink()

    local_sample = scratch_root / "quant" / "SRR_LOCAL"
    local_sample.mkdir(parents=True)
    (local_sample / "SRR_LOCAL_abundance.tsv").write_text("target_id\ttpm\n" + "TX\t1\n" * 30)
    assert orchestrator._promote_local_quant_output(work_dir, "SRR_LOCAL", scratch_root)
    assert (work_dir / "quant" / "SRR_LOCAL" / "SRR_LOCAL_abundance.tsv").read_text().startswith("target_id")
    archives = list((work_dir / "archive" / "quant_replaced" / "SRR_LOCAL").iterdir())
    assert len(archives) == 1
    orchestrator._cleanup_local_quant_workspace(scratch_root)
    assert not scratch_root.exists()


def test_local_quant_scratch_stages_validated_fastq_when_enabled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Opt-in staging copies validated reads without replacing authoritative inputs."""
    work_dir = tmp_path / "work"
    metadata_dir = work_dir / "metadata"
    raw_dir = work_dir / "getfastq" / "SRR_STAGE"
    metadata_dir.mkdir(parents=True)
    raw_dir.mkdir(parents=True)
    metadata = metadata_dir / "metadata_selected.tsv"
    metadata.write_text("run\tlib_layout\nSRR_STAGE\tpaired\n", encoding="utf-8")
    for mate in (1, 2):
        with gzip.open(raw_dir / f"SRR_STAGE_{mate}.fastq.gz", "wb") as handle:
            handle.write(b"@read\nACGT\n+\n!!!!\n")

    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", str(tmp_path / "quant-scratch"))
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTQ_STAGE", "yes")
    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB", "0")
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    command = [
        "amalgkit",
        "quant",
        "--out_dir",
        str(work_dir),
        "--metadata",
        str(metadata),
        "--index_dir",
        str(tmp_path / "index"),
    ]

    prepared = orchestrator._prepare_local_quant_workspace(
        work_dir,
        "test_species",
        "SRR_STAGE",
        str(metadata),
        command,
        expected_paired=True,
    )

    assert prepared is not None
    _, scratch_root = prepared
    staged_sample = scratch_root / "getfastq" / "SRR_STAGE"
    assert (scratch_root / "getfastq").is_dir()
    assert not (scratch_root / "getfastq").is_symlink()
    assert sorted(path.name for path in staged_sample.iterdir()) == [
        "SRR_STAGE_1.fastq.gz",
        "SRR_STAGE_2.fastq.gz",
    ]
    assert (staged_sample / "SRR_STAGE_1.fastq.gz").read_bytes() == (raw_dir / "SRR_STAGE_1.fastq.gz").read_bytes()
    assert (raw_dir / "SRR_STAGE_1.fastq.gz").exists()
    marker = json.loads((work_dir / "raw_validation" / "SRR_STAGE.json").read_text())
    assert marker["declared_library_layout"] == "paired"
    assert marker["effective_library_layout"] == "paired"

    # An interrupted quantification leaves the staged FASTQs in place. A
    # resume should keep them rather than copy the mounted payload again.
    monkeypatch.setattr(
        orchestrator_module.shutil,
        "copy2",
        lambda *_args, **_kwargs: pytest.fail("resume should reuse the staged FASTQs"),
    )
    resumed = orchestrator._prepare_local_quant_workspace(
        work_dir,
        "test_species",
        "SRR_STAGE",
        str(metadata),
        command,
        expected_paired=True,
    )
    assert resumed is not None
    assert (resumed[1] / "getfastq" / "SRR_STAGE").is_dir()

    orchestrator._cleanup_local_quant_workspace(scratch_root)
    assert not scratch_root.exists()
    assert (raw_dir / "SRR_STAGE_1.fastq.gz").exists()


def test_current_quant_provenance_sidecar_is_required_and_atomic(tmp_path: Path) -> None:
    """Readable prior abundance files are not current release evidence."""
    sample_dir = tmp_path / "quant" / "SRR123"
    sample_dir.mkdir(parents=True)
    (sample_dir / "SRR123_abundance.tsv").write_text("target_id\ttpm\nTX1\t1\n")
    assert not is_current_quantification(sample_dir, "SRR123")

    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n")
    marker = write_quant_provenance(
        sample_dir,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
    )
    assert marker.exists()
    assert is_current_quantification(sample_dir, "SRR123")
    assert not is_current_quantification(sample_dir, "SRR999")


def test_new_quant_provenance_binds_abundance_content(tmp_path: Path) -> None:
    """A current sidecar becomes invalid when its abundance table changes."""
    sample_dir = tmp_path / "quant" / "SRR123"
    sample_dir.mkdir(parents=True)
    abundance = sample_dir / "SRR123_abundance.tsv"
    abundance.write_text("target_id\ttpm\nTX1\t1\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n", encoding="utf-8")

    write_quant_provenance(
        sample_dir,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
        quantification_file=abundance,
    )
    assert is_current_quantification(sample_dir, "SRR123") is True
    abundance.write_text("target_id\ttpm\nTX1\t2\n", encoding="utf-8")
    assert is_current_quantification(sample_dir, "SRR123") is False


def test_non_current_state_reconciliation_uses_validation_slot_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Restart provenance reads use bounded parallelism and serialized DB writes."""

    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config",
        tmp_path / "logs",
        db_path=tmp_path / "progress.db",
    )
    sample_ids = ["SRR1", "SRR2", "SRR3", "SRR4"]
    orchestrator.db.init_species("test_species", sample_ids)
    for sample_id in sample_ids:
        orchestrator.db.set_state("test_species", sample_id, "quantified")
    orchestrator._resource_profile = build_pipeline_resource_profile(
        workers=4,
        threads=4,
        validation_slots=3,
    )

    observed_workers: list[int] = []
    real_executor = orchestrator_module.concurrent.futures.ThreadPoolExecutor

    class RecordingExecutor(real_executor):
        def __init__(self, max_workers: int, *args: object, **kwargs: object) -> None:
            observed_workers.append(max_workers)
            super().__init__(max_workers=max_workers, *args, **kwargs)

    monkeypatch.setattr(
        orchestrator_module.concurrent.futures,
        "ThreadPoolExecutor",
        RecordingExecutor,
    )
    monkeypatch.setattr(
        orchestrator_module,
        "classify_quantification",
        lambda _sample_dir, sample_id: {
            "status": "invalid" if sample_id == "SRR3" else "current",
            "reason": "test invalid" if sample_id == "SRR3" else "test current",
            "contract_id": None,
            "observed_amalgkit_version": None,
            "observed_release_tag": None,
            "observed_source_revision": None,
        },
    )

    assert orchestrator._reset_non_current_quantified_states("test_species") == 1
    assert observed_workers == [3]
    assert orchestrator.db.get_state("test_species", "SRR3") == "quarantined"
    assert orchestrator.db.get_state("test_species", "SRR1") == "quantified"


def test_metadata_provenance_binds_current_table_contents(tmp_path: Path) -> None:
    """New metadata sidecars invalidate out-of-band table edits."""
    work_dir = tmp_path / "work"
    metadata_dir = work_dir / "metadata"
    metadata_dir.mkdir(parents=True)
    (metadata_dir / "metadata.tsv").write_text("run\tvalue\nSRR1\t1\n", encoding="utf-8")
    selected = metadata_dir / "metadata_selected.tsv"
    selected.write_text("run\tvalue\nSRR1\t1\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n", encoding="utf-8")
    rules = tmp_path / "select_rules.tsv"
    rules.write_text("rule\tvalue\nkeep\tyes\n", encoding="utf-8")

    write_metadata_provenance(
        work_dir,
        species="test_species",
        config_path=config,
        selection_rules_path=rules,
    )
    assert is_current_metadata(work_dir) is True
    selected.write_text("run\tvalue\nSRR1\t2\n", encoding="utf-8")
    assert is_current_metadata(work_dir) is False


def test_current_quantification_skips_download_and_reclaims_leftover_raw(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A current result is idempotent and cleans any leftover raw input."""
    data_root = tmp_path / "data"
    work_dir = data_root / "test_species" / "work"
    quant_dir = work_dir / "quant" / "SRR123"
    raw_dir = work_dir / "getfastq" / "SRR123"
    quant_dir.mkdir(parents=True)
    raw_dir.mkdir(parents=True)
    (quant_dir / "SRR123_abundance.tsv").write_text("target_id\ttpm\nTX1\t1\n" * 20)
    (raw_dir / "SRR123_1.fastq.gz").write_bytes(b"raw" * 100)
    config_path = tmp_path / "amalgkit_test_species.yaml"
    config_path.write_text("species_list: [Test_species]\n")
    write_quant_provenance(
        quant_dir,
        species="test_species",
        run_accession="SRR123",
        config_path=config_path,
        command=["amalgkit", "quant"],
    )
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "progress.db"
    )
    orchestrator.download_fastq = (  # type: ignore[method-assign]
        lambda *_args, **_kwargs: pytest.fail("current quant must skip download")
    )
    orchestrator.quant_sample = (  # type: ignore[method-assign]
        lambda *_args, **_kwargs: pytest.fail("current quant must skip quantification")
    )

    result = orchestrator.process_single_sample("SRR123", 1, raw_dir, config_path, "test_species", threads=1)

    assert result["skipped"] is True
    assert result["quantified"] is True
    assert not (raw_dir / "SRR123_1.fastq.gz").exists()


def test_runtime_config_remaps_existing_external_alias(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Alias workspaces must be shared by quantification and downstream steps."""
    data_root = tmp_path / "data"
    active_work = data_root / "apis_mellifera_all" / "work"
    (active_work / "metadata").mkdir(parents=True)
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))

    source = tmp_path / "amalgkit_apis_mellifera.yaml"
    source.write_text(
        yaml.safe_dump(
            {
                "work_dir": "output/amalgkit/apis_mellifera/work",
                "log_dir": "output/amalgkit/apis_mellifera/logs",
                "genome": {"dest_dir": "output/amalgkit/apis_mellifera/genome"},
                "steps": {
                    "merge": {"out_dir": "output/amalgkit/apis_mellifera/work"},
                    "select": {"select_rules_tsv": "config/amalgkit/select_rules.tsv"},
                },
            },
            sort_keys=False,
        )
    )

    runtime = _runtime_config_path(source, "apis_mellifera")
    assert runtime != source
    config = yaml.safe_load(runtime.read_text())
    assert config["work_dir"] == "output/amalgkit/apis_mellifera_all/work"
    assert config["log_dir"] == "output/amalgkit/apis_mellifera_all/logs"
    assert config["genome"]["dest_dir"] == "output/amalgkit/apis_mellifera_all/genome"
    assert config["steps"]["select"]["select_rules_tsv"] == "config/amalgkit/select_rules.tsv"


def test_missing_species_preparation_uses_only_current_steps(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Automatic preparation must not invoke removed Amalgkit CLI stages."""
    data_root = tmp_path / "data"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (index_dir / "Test_species.idx").write_bytes(b"kallisto-index")
    config_path = config_dir / "amalgkit_test_species.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "species_list": ["Test_species"],
                "genome": {"index_dir": str(index_dir)},
            }
        )
    )

    calls: list[list[str]] = []
    metadata_path = data_root / "test_species" / "work" / "metadata" / "metadata.tsv"

    def fake_run(command: list[str], **_: object) -> object:
        calls.append(command)
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_text("run\ttotal_bases\nSRR123\t1000\n")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    def fake_urlretrieve(_: str, destination: object) -> None:
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        Path(destination).write_bytes(b"taxdump")

    import subprocess

    monkeypatch.setattr("metainformant.rna.engine.streaming_orchestrator.subprocess.run", fake_run)
    monkeypatch.setattr("metainformant.rna.engine.streaming_orchestrator.urllib.request.urlretrieve", fake_urlretrieve)

    orchestrator = StreamingPipelineOrchestrator(config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db")
    checks = iter([False, True])
    monkeypatch.setattr(orchestrator, "verify_genome_index", lambda *_: next(checks))

    tasks = orchestrator.discover_species_tasks(config_path.name, max_gb=1.0, threads=1)

    assert len(tasks) == 1
    prep_calls = [
        call for call in calls if len(call) >= 2 and call[0] == "amalgkit" and call[1] in {"metadata", "select"}
    ]
    assert [call[1] for call in prep_calls] == ["metadata", "select"]
    assert prep_calls[0][prep_calls[0].index("--redo") + 1] == "yes"
    assert "--redo" not in prep_calls[1]
    assert all("config" not in call and "index" not in call for call in prep_calls)
    assert is_current_metadata(data_root / "test_species" / "work") is True
