"""Comprehensive tests for the RNA workflow cleanup module.

Tests for metainformant.rna.engine.workflow_cleanup functions using REAL file
operations only -- no mocking, no monkeypatching, no stubs.

Every test creates real files on disk via tmp_path and exercises the actual
function logic.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any, Dict

import pytest

from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig
from metainformant.rna.engine.workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_after_quant,
    cleanup_fastqs,
    cleanup_incorrectly_placed_sra_files,
    cleanup_temp_files,
    filter_metadata_for_unquantified,
    get_quantified_samples,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def workflow_config(tmp_path: Path) -> AmalgkitWorkflowConfig:
    """Create a minimal AmalgkitWorkflowConfig rooted in tmp_path."""
    work_dir = tmp_path / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    config = AmalgkitWorkflowConfig(
        work_dir=str(work_dir),
        threads=4,
        species_list=["test_species"],
    )
    return config


@pytest.fixture
def workflow_config_with_steps(tmp_path: Path) -> AmalgkitWorkflowConfig:
    """Config that sets custom step output directories via extra_config."""
    work_dir = tmp_path / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    config = AmalgkitWorkflowConfig(
        work_dir=str(work_dir),
        threads=4,
        species_list=["test_species"],
        steps={
            "getfastq": {"out_dir": str(work_dir / "custom_fastq")},
            "quant": {"out_dir": str(work_dir / "custom_quant")},
        },
    )
    return config


def _create_abundance(quant_dir: Path, sample_id: str, size: int = 200) -> Path:
    """Helper: create a realistic abundance.tsv file for a sample."""
    sample_dir = quant_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    abundance_file = sample_dir / "abundance.tsv"
    # Write enough content to exceed the 100-byte threshold used in get_quantified_samples
    header = "target_id\tlength\teff_length\test_counts\ttpm\n"
    row = "gene_0001\t1000\t800.5\t42.0\t3.14\n"
    content = header + (row * max(1, (size - len(header)) // len(row)))
    abundance_file.write_text(content)
    return abundance_file


def _create_fastq_files(fastq_dir: Path, sample_id: str) -> Path:
    """Helper: create a sample directory with FASTQ files."""
    sample_dir = fastq_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    (sample_dir / f"{sample_id}_1.fastq.gz").write_bytes(b"@SEQ\nACGT\n+\nIIII\n" * 100)
    (sample_dir / f"{sample_id}_2.fastq.gz").write_bytes(b"@SEQ\nTGCA\n+\nIIII\n" * 100)
    return sample_dir


def _create_metadata_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Helper: write a metadata TSV with a 'run' column."""
    fieldnames = list(rows[0].keys()) if rows else ["run", "organism", "layout"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ===========================================================================
# check_disk_space
# ===========================================================================


class TestCheckDiskSpace:
    """Tests for check_disk_space."""

    def test_returns_tuple(self, tmp_path: Path) -> None:
        """check_disk_space returns a (bool, float) tuple."""
        result = check_disk_space(tmp_path)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], bool)
        assert isinstance(result[1], float)

    def test_sufficient_space_with_low_threshold(self, tmp_path: Path) -> None:
        """With a near-zero threshold any real disk should pass."""
        ok, free_gb = check_disk_space(tmp_path, min_free_gb=0.0001)
        assert ok is True
        assert free_gb > 0.0

    def test_free_gb_positive(self, tmp_path: Path) -> None:
        """Reported free GB should be a positive number on a real filesystem."""
        _, free_gb = check_disk_space(tmp_path)
        assert free_gb > 0

    def test_nonexistent_path_returns_default(self) -> None:
        """Non-existent path triggers the except branch, returns (True, 100.0)."""
        result = check_disk_space(Path("/this/path/does/not/exist/at/all"))
        assert result == (True, 100.0)


# ===========================================================================
# check_disk_space_or_fail
# ===========================================================================


class TestCheckDiskSpaceOrFail:
    """Tests for check_disk_space_or_fail."""

    def test_returns_free_gb_when_ok(self, tmp_path: Path) -> None:
        """Should return free GB when space is sufficient."""
        free_gb = check_disk_space_or_fail(tmp_path, min_free_gb=0.0001, step_name="test_step")
        assert isinstance(free_gb, float)
        assert free_gb > 0

    def test_does_not_raise_with_low_threshold(self, tmp_path: Path) -> None:
        """No exception when threshold is trivially low."""
        # Should not raise
        check_disk_space_or_fail(tmp_path, min_free_gb=0.0, step_name="quant")

    def test_nonexistent_path_assumed_ok(self) -> None:
        """Non-existent path falls back to (True, 100.0) so no exception."""
        free_gb = check_disk_space_or_fail(
            Path("/nonexistent/path/xyz"),
            min_free_gb=50.0,
            step_name="getfastq",
        )
        assert free_gb == 100.0


# ===========================================================================
# cleanup_temp_files
# ===========================================================================


class TestCleanupTempFiles:
    """Tests for cleanup_temp_files."""

    def test_nonexistent_dir_noop(self, tmp_path: Path) -> None:
        """Should return immediately for a directory that doesn't exist."""
        missing = tmp_path / "no_such_dir"
        cleanup_temp_files(missing, max_size_gb=0.0)
        assert not missing.exists()

    def test_under_threshold_keeps_files(self, tmp_path: Path) -> None:
        """Files below the size threshold should be left alone."""
        temp_dir = tmp_path / "temp"
        temp_dir.mkdir()
        (temp_dir / "small_file.txt").write_text("hello")
        (temp_dir / "another.txt").write_text("world")

        # 50 GB threshold -- our tiny files will never hit it
        cleanup_temp_files(temp_dir, max_size_gb=50.0)

        assert (temp_dir / "small_file.txt").exists()
        assert (temp_dir / "another.txt").exists()

    def test_over_threshold_deletes_contents(self, tmp_path: Path) -> None:
        """When total size exceeds threshold, all contents are removed."""
        temp_dir = tmp_path / "temp"
        temp_dir.mkdir()

        # Create files and a subdirectory
        (temp_dir / "file_a.dat").write_bytes(b"x" * 1024)
        sub = temp_dir / "subdir"
        sub.mkdir()
        (sub / "file_b.dat").write_bytes(b"y" * 512)

        # Use a threshold of ~0 GB so even small files exceed it
        cleanup_temp_files(temp_dir, max_size_gb=0.000001)

        # Directory itself should still exist; contents should be gone
        assert temp_dir.exists()
        assert list(temp_dir.iterdir()) == []

    def test_empty_dir_under_threshold(self, tmp_path: Path) -> None:
        """Empty directory (0 bytes) should stay untouched."""
        temp_dir = tmp_path / "temp"
        temp_dir.mkdir()

        cleanup_temp_files(temp_dir, max_size_gb=0.0)
        assert temp_dir.exists()


# ===========================================================================
# cleanup_incorrectly_placed_sra_files
# ===========================================================================


class TestCleanupIncorrectlyPlacedSraFiles:
    """Tests for cleanup_incorrectly_placed_sra_files.

    The function scans ~/ncbi/public/sra and /tmp/ncbi/public/sra.
    We cannot easily create files in those system directories from a test,
    so we verify the function handles non-existent default locations gracefully.
    """

    def test_no_crash_when_default_locations_missing(self, tmp_path: Path) -> None:
        """Function should complete without error when default SRA dirs don't exist."""
        getfastq_dir = tmp_path / "getfastq"
        getfastq_dir.mkdir()
        # This should not raise -- the default locations almost certainly don't exist
        # or if they do, the function just scans them harmlessly
        cleanup_incorrectly_placed_sra_files(getfastq_dir)
        assert getfastq_dir.exists()

    def test_getfastq_dir_created_if_needed(self, tmp_path: Path) -> None:
        """The target directory should still exist after the call."""
        getfastq_dir = tmp_path / "getfastq"
        getfastq_dir.mkdir(parents=True, exist_ok=True)
        cleanup_incorrectly_placed_sra_files(getfastq_dir)
        assert getfastq_dir.is_dir()


# ===========================================================================
# cleanup_fastqs
# ===========================================================================


class TestCleanupFastqs:
    """Tests for cleanup_fastqs."""

    def test_cleans_primary_fastq_dir(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Should remove sample dirs under work_dir/fastq/getfastq/<sample>."""
        fastq_dir = workflow_config.work_dir / "fastq" / "getfastq"
        sample_dir = _create_fastq_files(fastq_dir, "SRR111111")

        assert sample_dir.exists()
        cleanup_fastqs(workflow_config, ["SRR111111"])
        assert not sample_dir.exists()

    def test_cleans_alt_fastq_dir(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Should also clean work_dir/fastq/<sample>."""
        alt_dir = workflow_config.work_dir / "fastq"
        sample_dir = _create_fastq_files(alt_dir, "SRR222222")

        cleanup_fastqs(workflow_config, ["SRR222222"])
        assert not sample_dir.exists()

    def test_cleans_getfastq_dir(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Should also clean work_dir/getfastq/<sample>."""
        gf_dir = workflow_config.work_dir / "getfastq"
        sample_dir = _create_fastq_files(gf_dir, "ERR333333")

        cleanup_fastqs(workflow_config, ["ERR333333"])
        assert not sample_dir.exists()

    def test_skips_empty_sample_ids(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Empty strings in sample_ids should be silently skipped."""
        cleanup_fastqs(workflow_config, ["", "", ""])
        # No error should have occurred

    def test_missing_sample_dir_no_error(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Non-existent sample dirs should be handled gracefully."""
        cleanup_fastqs(workflow_config, ["SRR_NOT_REAL"])

    def test_multiple_samples(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Should clean all listed samples."""
        fastq_base = workflow_config.work_dir / "fastq" / "getfastq"
        dirs = []
        for sid in ["SRR100001", "SRR100002", "SRR100003"]:
            dirs.append(_create_fastq_files(fastq_base, sid))

        cleanup_fastqs(workflow_config, ["SRR100001", "SRR100002", "SRR100003"])
        for d in dirs:
            assert not d.exists()

    def test_per_step_getfastq_dir(self, tmp_path: Path) -> None:
        """When per_step defines a getfastq out_dir, that path is also cleaned."""
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        custom_out = tmp_path / "custom_output"
        config = AmalgkitWorkflowConfig(
            work_dir=str(work_dir),
            threads=2,
            species_list=["sp"],
            per_step={"getfastq": {"out_dir": str(custom_out)}},
        )

        # Create files in the custom per_step path
        sample_dir = custom_out / "getfastq" / "DRR444444"
        sample_dir.mkdir(parents=True)
        (sample_dir / "DRR444444_1.fastq.gz").write_bytes(b"data")

        cleanup_fastqs(config, ["DRR444444"])
        assert not sample_dir.exists()


# ===========================================================================
# get_quantified_samples
# ===========================================================================


class TestGetQuantifiedSamples:
    """Tests for get_quantified_samples."""

    def test_empty_when_no_quant_dir(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """No quant directory => empty set."""
        result = get_quantified_samples(workflow_config)
        assert result == set()

    def test_detects_srr_abundance(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Finds SRR-prefixed samples with abundance.tsv."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR999001")
        _create_abundance(quant_dir, "SRR999002")

        result = get_quantified_samples(workflow_config)
        assert result == {"SRR999001", "SRR999002"}

    def test_detects_err_and_drr_prefixes(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Samples with ERR/DRR prefixes should also be detected."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "ERR000001")
        _create_abundance(quant_dir, "DRR000002")

        result = get_quantified_samples(workflow_config)
        assert "ERR000001" in result
        assert "DRR000002" in result

    def test_ignores_small_abundance_files(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Abundance files <= 100 bytes should be ignored (incomplete quant)."""
        quant_dir = workflow_config.work_dir / "quant"
        sample_dir = quant_dir / "SRR000099"
        sample_dir.mkdir(parents=True)
        tiny = sample_dir / "abundance.tsv"
        tiny.write_text("hdr\nrow\n")  # well under 100 bytes
        assert tiny.stat().st_size <= 100

        result = get_quantified_samples(workflow_config)
        assert "SRR000099" not in result

    def test_ignores_non_sra_prefixed_dirs(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Directories not starting with SRR/ERR/DRR should be skipped."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SAMPLE_001")  # wrong prefix

        result = get_quantified_samples(workflow_config)
        assert "SAMPLE_001" not in result

    def test_custom_quant_dir_via_steps(self, workflow_config_with_steps: AmalgkitWorkflowConfig) -> None:
        """When steps config specifies a custom quant out_dir, it is used."""
        config = workflow_config_with_steps
        custom_quant = Path(config.extra_config["steps"]["quant"]["out_dir"])
        _create_abundance(custom_quant, "SRR888001")

        result = get_quantified_samples(config)
        assert "SRR888001" in result

    def test_named_abundance_pattern(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Detects *_abundance.tsv pattern (alternative naming)."""
        quant_dir = workflow_config.work_dir / "quant"
        sample_dir = quant_dir / "SRR777001"
        sample_dir.mkdir(parents=True)
        ab = sample_dir / "SRR777001_abundance.tsv"
        ab.write_text("target_id\tlength\teff_length\test_counts\ttpm\n" + "g\t100\t90\t10\t1.0\n" * 10)

        result = get_quantified_samples(workflow_config)
        assert "SRR777001" in result


# ===========================================================================
# cleanup_after_quant
# ===========================================================================


class TestCleanupAfterQuant:
    """Tests for cleanup_after_quant."""

    def _setup_quantified_sample(
        self,
        config: AmalgkitWorkflowConfig,
        sample_id: str,
    ) -> tuple[Path, Path]:
        """Create both quant results and FASTQ files for a sample."""
        quant_dir = config.work_dir / "quant"
        _create_abundance(quant_dir, sample_id)

        fastq_dir = config.work_dir / "fastq" / "getfastq"
        sample_fastq_dir = _create_fastq_files(fastq_dir, sample_id)

        return quant_dir / sample_id, sample_fastq_dir

    def test_dry_run_does_not_delete(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Dry run should report files but leave them on disk."""
        _, fastq_dir = self._setup_quantified_sample(workflow_config, "SRR600001")

        result = cleanup_after_quant(workflow_config, dry_run=True)
        assert result["samples_cleaned"] >= 1
        assert result["fastq_files_deleted"] >= 1
        # Files must still exist
        assert fastq_dir.exists()
        fastq_files = list(fastq_dir.glob("*.fastq.gz"))
        assert len(fastq_files) == 2

    def test_actual_deletion(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Non-dry-run should actually delete the FASTQ files."""
        _, fastq_dir = self._setup_quantified_sample(workflow_config, "SRR600002")

        result = cleanup_after_quant(workflow_config, dry_run=False)
        assert result["samples_cleaned"] >= 1
        assert result["bytes_freed"] > 0
        # FASTQ files should be gone
        remaining = list(fastq_dir.glob("*.fastq.gz")) if fastq_dir.exists() else []
        assert remaining == []

    def test_no_quantified_samples(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """When no samples are quantified, nothing happens."""
        result = cleanup_after_quant(workflow_config)
        assert result["samples_cleaned"] == 0
        assert result["fastq_files_deleted"] == 0
        assert result["sra_files_deleted"] == 0

    def test_sra_files_counted_separately(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """SRA files (.sra, .sra.part) should be counted in sra_files_deleted."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR600003")

        fastq_base = workflow_config.work_dir / "fastq" / "getfastq" / "SRR600003"
        fastq_base.mkdir(parents=True)
        (fastq_base / "SRR600003.sra").write_bytes(b"sra_data" * 50)
        (fastq_base / "SRR600003.sra.part").write_bytes(b"partial" * 30)
        (fastq_base / "SRR600003_1.fastq.gz").write_bytes(b"fq_data" * 40)

        result = cleanup_after_quant(workflow_config, dry_run=False)
        assert result["sra_files_deleted"] >= 1
        assert result["fastq_files_deleted"] >= 1

    def test_cleanup_returns_proper_structure(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Return value should have all expected keys."""
        result = cleanup_after_quant(workflow_config)
        expected_keys = {"samples_cleaned", "fastq_files_deleted", "sra_files_deleted", "bytes_freed", "errors"}
        assert set(result.keys()) == expected_keys

    def test_missing_fastq_dir_no_error(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """If FASTQ dir doesn't exist but quant does, no crash."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR600004")
        # Deliberately do NOT create fastq directory

        result = cleanup_after_quant(workflow_config, dry_run=False)
        # Should return without error
        assert result["errors"] == []


# ===========================================================================
# filter_metadata_for_unquantified
# ===========================================================================


class TestFilterMetadataForUnquantified:
    """Tests for filter_metadata_for_unquantified."""

    def test_no_quantified_copies_source(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """When nothing is quantified, the source is copied as-is."""
        source = tmp_path / "source.tsv"
        output = tmp_path / "output.tsv"
        _create_metadata_tsv(
            source,
            [
                {"run": "SRR000001", "organism": "honeybee", "layout": "PAIRED"},
                {"run": "SRR000002", "organism": "honeybee", "layout": "SINGLE"},
            ],
        )

        count = filter_metadata_for_unquantified(workflow_config, source, output)
        assert count == 2
        assert output.exists()

    def test_filters_quantified_samples(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """Quantified samples should be excluded from the output."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR000001")

        source = tmp_path / "source.tsv"
        output = tmp_path / "output.tsv"
        _create_metadata_tsv(
            source,
            [
                {"run": "SRR000001", "organism": "bee", "layout": "PAIRED"},
                {"run": "SRR000002", "organism": "bee", "layout": "PAIRED"},
                {"run": "SRR000003", "organism": "bee", "layout": "SINGLE"},
            ],
        )

        count = filter_metadata_for_unquantified(workflow_config, source, output)
        assert count == 2

        # Verify contents
        with open(output, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            runs = [row["run"] for row in reader]
        assert "SRR000001" not in runs
        assert "SRR000002" in runs
        assert "SRR000003" in runs

    def test_all_quantified_returns_zero(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """When all samples are quantified, returns 0 and does not write output."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR000010")
        _create_abundance(quant_dir, "SRR000011")

        source = tmp_path / "source.tsv"
        output = tmp_path / "output.tsv"
        _create_metadata_tsv(
            source,
            [
                {"run": "SRR000010", "organism": "bee", "layout": "PAIRED"},
                {"run": "SRR000011", "organism": "bee", "layout": "PAIRED"},
            ],
        )

        count = filter_metadata_for_unquantified(workflow_config, source, output)
        assert count == 0

    def test_missing_source_returns_zero(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """Missing source metadata file returns 0 when nothing is quantified."""
        source = tmp_path / "nonexistent.tsv"
        output = tmp_path / "output.tsv"
        # With no quantified samples, the function tries to copy source and
        # then count lines -- if source doesn't exist it returns 0
        count = filter_metadata_for_unquantified(workflow_config, source, output)
        assert count == 0

    def test_output_preserves_all_columns(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """The filtered output should retain all original TSV columns."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR000020")

        source = tmp_path / "source.tsv"
        output = tmp_path / "output.tsv"
        _create_metadata_tsv(
            source,
            [
                {"run": "SRR000020", "organism": "ant", "layout": "PAIRED"},
                {"run": "SRR000021", "organism": "wasp", "layout": "SINGLE"},
            ],
        )

        filter_metadata_for_unquantified(workflow_config, source, output)

        with open(output, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["run"] == "SRR000021"
        assert rows[0]["organism"] == "wasp"
        assert rows[0]["layout"] == "SINGLE"

    def test_empty_metadata_no_data_rows(self, workflow_config: AmalgkitWorkflowConfig, tmp_path: Path) -> None:
        """Metadata with header only and no data rows."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR000030")

        source = tmp_path / "source.tsv"
        output = tmp_path / "output.tsv"
        # Write header only
        source.write_text("run\torganism\tlayout\n")

        count = filter_metadata_for_unquantified(workflow_config, source, output)
        assert count == 0


# ===========================================================================
# Edge-case / Integration tests
# ===========================================================================


class TestEdgeCases:
    """Edge-case and cross-function integration tests."""

    def test_check_disk_space_with_file_path(self, tmp_path: Path) -> None:
        """check_disk_space should handle a path that is a file (its parent volume)."""
        file_path = tmp_path / "test_file.txt"
        file_path.write_text("content")
        # The function calls os.statvfs which works on files too
        ok, free_gb = check_disk_space(file_path)
        assert isinstance(ok, bool)
        assert isinstance(free_gb, float)

    def test_cleanup_temp_files_nested_structure(self, tmp_path: Path) -> None:
        """Deeply nested directories should be fully removed when over threshold."""
        temp_dir = tmp_path / "deep_temp"
        deep = temp_dir / "a" / "b" / "c" / "d"
        deep.mkdir(parents=True)
        (deep / "file.dat").write_bytes(b"z" * 2048)
        (temp_dir / "top.txt").write_text("top level")

        cleanup_temp_files(temp_dir, max_size_gb=0.0000001)
        assert temp_dir.exists()
        assert list(temp_dir.iterdir()) == []

    def test_get_quantified_then_cleanup(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Integration: get_quantified_samples feeds into cleanup_fastqs."""
        quant_dir = workflow_config.work_dir / "quant"
        fastq_base = workflow_config.work_dir / "fastq" / "getfastq"

        for sid in ["SRR500001", "SRR500002"]:
            _create_abundance(quant_dir, sid)
            _create_fastq_files(fastq_base, sid)

        quantified = get_quantified_samples(workflow_config)
        assert len(quantified) == 2

        cleanup_fastqs(workflow_config, list(quantified))
        for sid in quantified:
            assert not (fastq_base / sid).exists()

    def test_cleanup_after_quant_with_multiple_fastq_extensions(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """cleanup_after_quant should handle .fq, .fq.gz, .amalgkit.fastq.gz extensions."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR700001")

        fastq_base = workflow_config.work_dir / "fastq" / "getfastq" / "SRR700001"
        fastq_base.mkdir(parents=True)
        (fastq_base / "SRR700001_1.fq").write_bytes(b"fq1" * 50)
        (fastq_base / "SRR700001_2.fq.gz").write_bytes(b"fq2" * 50)
        (fastq_base / "SRR700001.amalgkit.fastq.gz").write_bytes(b"amal" * 50)
        (fastq_base / "SRR700001.fastq").write_bytes(b"plain" * 50)

        result = cleanup_after_quant(workflow_config, dry_run=False)
        assert result["fastq_files_deleted"] == 4
        remaining = list(fastq_base.glob("*")) if fastq_base.exists() else []
        assert remaining == []

    def test_dry_run_bytes_freed_matches_actual(self, workflow_config: AmalgkitWorkflowConfig) -> None:
        """Dry-run should report the same bytes_freed as an actual run would."""
        quant_dir = workflow_config.work_dir / "quant"
        _create_abundance(quant_dir, "SRR800001")

        fastq_base = workflow_config.work_dir / "fastq" / "getfastq" / "SRR800001"
        fastq_base.mkdir(parents=True)
        data = b"x" * 4096
        (fastq_base / "SRR800001_1.fastq.gz").write_bytes(data)
        (fastq_base / "SRR800001_2.fastq.gz").write_bytes(data)

        dry_result = cleanup_after_quant(workflow_config, dry_run=True)
        # Files still exist after dry run, so we can run again for real
        actual_result = cleanup_after_quant(workflow_config, dry_run=False)

        assert dry_result["bytes_freed"] == actual_result["bytes_freed"]
        assert dry_result["fastq_files_deleted"] == actual_result["fastq_files_deleted"]
