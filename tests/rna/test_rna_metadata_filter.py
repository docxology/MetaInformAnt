"""Tests for RNA-seq metadata filtering utilities.

This module tests the metadata filter helpers that reduce a full ENA/SRA
metadata table down to the samples selected for downstream processing.

All tests follow real-implementation policy and use real TSV files on disk
via pytest's tmp_path fixture.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from metainformant.rna.amalgkit.metadata_filter import (
    filter_selected_metadata,
    get_sample_count,
    validate_filtered_metadata,
)


def _write_metadata(path: Path, rows: list[dict[str, str]]) -> Path:
    """Write a real metadata TSV to disk from row dicts."""

    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    return path


class TestFilterSelectedMetadata:
    """Test filter_selected_metadata() function."""

    def test_filters_to_selected_samples_and_preserves_original(self, tmp_path: Path):
        """Test default-flag filtering keeps only sampled/qualified/RNA rows."""
        metadata_path = _write_metadata(
            tmp_path / "metadata.tsv",
            [
                # Kept: sampled, qualified, no exclusion, RNA-Seq, TRANSCRIPTOMIC
                {
                    "run": "SRR000001",
                    "experiment": "SRX000001",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "",
                    "lib_strategy": "RNA-Seq",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                # Rejected: not sampled
                {
                    "run": "SRR000002",
                    "experiment": "SRX000002",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "no",
                    "is_qualified": "yes",
                    "exclusion": "",
                    "lib_strategy": "RNA-Seq",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                # Rejected: not qualified
                {
                    "run": "SRR000003",
                    "experiment": "SRX000003",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "yes",
                    "is_qualified": "no",
                    "exclusion": "",
                    "lib_strategy": "RNA-Seq",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                # Rejected: explicit exclusion reason
                {
                    "run": "SRR000004",
                    "experiment": "SRX000004",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "low quality reads",
                    "lib_strategy": "RNA-Seq",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                # Rejected: non-RNA library strategy
                {
                    "run": "SRR000005",
                    "experiment": "SRX000005",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "",
                    "lib_strategy": "WGS",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                # Rejected: non-TRANSCRIPTOMIC library source
                {
                    "run": "SRR000006",
                    "experiment": "SRX000006",
                    "scientific_name": "Homo sapiens",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "",
                    "lib_strategy": "RNA-Seq",
                    "lib_source": "GENOMIC",
                },
            ],
        )

        original_bytes = metadata_path.read_bytes()
        output_path = filter_selected_metadata(metadata_path)

        # Default output path is metadata_selected.tsv next to the input
        assert output_path == metadata_path.parent / "metadata_selected.tsv"
        assert output_path.exists()

        filtered = pd.read_csv(output_path, sep="\t")
        # Exactly the single passing row survives
        assert list(filtered["run"]) == ["SRR000001"]
        assert list(filtered["is_sampled"]) == ["yes"]
        assert list(filtered["is_qualified"]) == ["yes"]
        assert list(filtered["lib_strategy"]) == ["RNA-Seq"]
        assert list(filtered["lib_source"]) == ["TRANSCRIPTOMIC"]

        # Original file is untouched
        assert metadata_path.read_bytes() == original_bytes
        assert get_sample_count(metadata_path) == 6

        # Row-count round trip agrees with get_sample_count on both files
        assert get_sample_count(output_path) == 1

        # The result satisfies the downstream validator
        assert validate_filtered_metadata(metadata_path, output_path, expected_count=1) is True

    def test_custom_output_path_and_keeps_cDNA_and_exclusion_no(self, tmp_path: Path):
        """Test explicit output_path plus cDNA strategy and exclusion=no acceptance."""
        metadata_path = _write_metadata(
            tmp_path / "custom" / "metadata.tsv",
            [
                {
                    "run": "SRR100001",
                    "experiment": "SRX100001",
                    "scientific_name": "Mus musculus",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "no",
                    "lib_strategy": "cDNA",
                    "lib_source": "TRANSCRIPTOMIC",
                },
                {
                    "run": "SRR100002",
                    "experiment": "SRX100002",
                    "scientific_name": "Mus musculus",
                    "is_sampled": "no",
                    "is_qualified": "yes",
                    "exclusion": "no",
                    "lib_strategy": "cDNA",
                    "lib_source": "TRANSCRIPTOMIC",
                },
            ],
        )
        custom_output = tmp_path / "elsewhere" / "filtered.tsv"
        custom_output.parent.mkdir(parents=True, exist_ok=True)

        output_path = filter_selected_metadata(metadata_path, output_path=custom_output)

        assert output_path == custom_output
        assert output_path.exists()
        assert output_path.parent == (tmp_path / "elsewhere")
        assert get_sample_count(output_path) == 1
        filtered = pd.read_csv(output_path, sep="\t")
        assert list(filtered["run"]) == ["SRR100001"]

    def test_missing_required_columns_raise_value_error(self, tmp_path: Path):
        """Test that absent is_sampled/is_qualified/exclusion columns are rejected."""
        metadata_path = _write_metadata(
            tmp_path / "metadata.tsv",
            [
                {
                    "run": "SRR200001",
                    "experiment": "SRX200001",
                    "scientific_name": "Homo sapiens",
                }
            ],
        )

        with pytest.raises(ValueError, match="Required columns missing"):
            filter_selected_metadata(metadata_path)

    def test_missing_input_file_raises_file_not_found_error(self, tmp_path: Path):
        """Test FileNotFoundError for a nonexistent metadata file."""
        with pytest.raises(FileNotFoundError, match="Metadata file not found"):
            filter_selected_metadata(tmp_path / "does_not_exist.tsv")


class TestGetSampleCount:
    """Test get_sample_count() function."""

    def test_counts_data_rows_excluding_header(self, tmp_path: Path):
        """Test that the returned count excludes the header row."""
        metadata_path = _write_metadata(
            tmp_path / "metadata.tsv",
            [
                {"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"},
                {"run": "SRR2", "experiment": "SRX2", "scientific_name": "B"},
                {"run": "SRR3", "experiment": "SRX3", "scientific_name": "C"},
            ],
        )

        assert get_sample_count(metadata_path) == 3

    def test_missing_file_raises_file_not_found_error(self, tmp_path: Path):
        """Test FileNotFoundError for a nonexistent metadata file."""
        with pytest.raises(FileNotFoundError, match="Metadata file not found"):
            get_sample_count(tmp_path / "does_not_exist.tsv")

    def test_corrupt_file_raises_value_error(self, tmp_path: Path):
        """Test ValueError when the file cannot be parsed as a TSV."""
        corrupt_path = tmp_path / "corrupt.tsv"
        corrupt_path.write_bytes(b"\x00\xff\xfe\x01binary garbage \x9c not utf-8 or tsv")

        with pytest.raises(ValueError, match="Failed to read metadata file"):
            get_sample_count(corrupt_path)


class TestValidateFilteredMetadata:
    """Test validate_filtered_metadata() function."""

    def test_same_row_count_original_and_filtered_returns_true(self, tmp_path: Path):
        """Test validation passes when the filter removed no rows."""
        rows = [
            {"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"},
            {"run": "SRR2", "experiment": "SRX2", "scientific_name": "B"},
        ]
        original_path = _write_metadata(tmp_path / "original.tsv", rows)
        filtered_path = _write_metadata(tmp_path / "filtered.tsv", rows)

        assert validate_filtered_metadata(original_path, filtered_path) is True

    def test_filtered_larger_than_original_raises_value_error(self, tmp_path: Path):
        """Test ValueError when the filtered file has more samples than the original."""
        original_path = _write_metadata(
            tmp_path / "original.tsv",
            [{"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"}],
        )
        filtered_path = _write_metadata(
            tmp_path / "filtered.tsv",
            [
                {"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"},
                {"run": "SRR2", "experiment": "SRX2", "scientific_name": "B"},
            ],
        )

        with pytest.raises(ValueError, match="more samples than original"):
            validate_filtered_metadata(original_path, filtered_path)

    def test_expected_count_mismatch_warns_but_returns_true(self, tmp_path: Path):
        """Test expected_count mismatch only warns (source behavior), still True."""
        rows = [
            {"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"},
            {"run": "SRR2", "experiment": "SRX2", "scientific_name": "B"},
        ]
        original_path = _write_metadata(tmp_path / "original.tsv", rows)
        filtered_path = _write_metadata(tmp_path / "filtered.tsv", rows)

        assert validate_filtered_metadata(original_path, filtered_path, expected_count=5) is True

    def test_missing_required_columns_raise_value_error(self, tmp_path: Path):
        """Test ValueError when the filtered file lacks run/experiment/scientific_name."""
        original_path = _write_metadata(
            tmp_path / "original.tsv",
            [{"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"}],
        )
        filtered_path = _write_metadata(tmp_path / "filtered.tsv", [{"run": "SRR1"}])

        with pytest.raises(ValueError, match="missing required columns"):
            validate_filtered_metadata(original_path, filtered_path)

    def test_empty_filtered_file_raises_value_error(self, tmp_path: Path):
        """Test ValueError when the filtered file contains no samples."""
        original_path = _write_metadata(
            tmp_path / "original.tsv",
            [{"run": "SRR1", "experiment": "SRX1", "scientific_name": "A"}],
        )
        empty_path = tmp_path / "empty.tsv"
        empty_path.write_text("run\texperiment\tscientific_name\n")

        with pytest.raises(ValueError, match="empty"):
            validate_filtered_metadata(original_path, empty_path)
