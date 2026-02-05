"""Tests for RNA metadata utilities.

Tests the metadata manipulation functions used in amalgkit workflows,
including deduplication and scientific name handling.
"""

from __future__ import annotations

import pytest
from pathlib import Path
import pandas as pd

from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata


class TestDeduplicateMetadata:
    """Test metadata deduplication functionality."""

    def test_deduplicate_removes_duplicate_runs(self, tmp_path: Path) -> None:
        """Test that duplicate runs are removed, keeping first occurrence."""
        # Create test metadata with duplicate runs
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR001", "SRR002", "SRR003", "SRR003"],
                "scientific_name": [
                    "Homo sapiens",
                    "Homo sapiens",
                    "Mus musculus",
                    "Rattus norvegicus",
                    "Rattus norvegicus",
                ],
                "sample_group": ["group1", "group1", "group2", "group3", "group3"],
            }
        )

        input_path = tmp_path / "metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        # Deduplicate
        result = deduplicate_metadata(input_path)

        assert result is True

        # Verify deduplication
        df_result = pd.read_csv(input_path, sep="\t")
        assert len(df_result) == 3  # Only 3 unique runs
        assert list(df_result["run"]) == ["SRR001", "SRR002", "SRR003"]

    def test_deduplicate_prefers_real_names_over_placeholders(self, tmp_path: Path) -> None:
        """Test that real scientific names are preferred over placeholders."""
        # Create metadata with placeholder names
        placeholder = "Please add in format: Genus species"
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR001", "SRR002", "SRR002"],
                "scientific_name": [placeholder, "Homo sapiens", "Mus musculus", placeholder],
                "sample_group": ["group1", "group1", "group2", "group2"],
            }
        )

        input_path = tmp_path / "metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        # Deduplicate
        result = deduplicate_metadata(input_path)

        assert result is True

        # Verify real names are kept
        df_result = pd.read_csv(input_path, sep="\t")
        assert len(df_result) == 2

        # SRR001 should have "Homo sapiens" (real name), not placeholder
        srr001_name = df_result[df_result["run"] == "SRR001"]["scientific_name"].values[0]
        assert srr001_name == "Homo sapiens"

        # SRR002 should have "Mus musculus" (real name)
        srr002_name = df_result[df_result["run"] == "SRR002"]["scientific_name"].values[0]
        assert srr002_name == "Mus musculus"

    def test_deduplicate_with_custom_output_path(self, tmp_path: Path) -> None:
        """Test deduplication with separate output path."""
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR001", "SRR002"],
                "scientific_name": ["Species A", "Species A", "Species B"],
            }
        )

        input_path = tmp_path / "input.tsv"
        output_path = tmp_path / "output.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path, output_path)

        assert result is True
        assert output_path.exists()

        # Input should be unchanged
        df_input = pd.read_csv(input_path, sep="\t")
        assert len(df_input) == 3

        # Output should be deduplicated
        df_output = pd.read_csv(output_path, sep="\t")
        assert len(df_output) == 2

    def test_deduplicate_handles_missing_file(self, tmp_path: Path) -> None:
        """Test that missing files return False."""
        missing_path = tmp_path / "nonexistent.tsv"

        result = deduplicate_metadata(missing_path)

        assert result is False

    def test_deduplicate_handles_empty_file(self, tmp_path: Path) -> None:
        """Test handling of empty metadata files."""
        # Create empty TSV with just headers
        input_path = tmp_path / "empty.tsv"
        with open(input_path, "w") as f:
            f.write("run\tscientific_name\n")

        result = deduplicate_metadata(input_path)

        assert result is True

        # Should still be valid but empty
        df = pd.read_csv(input_path, sep="\t")
        assert len(df) == 0

    def test_deduplicate_without_scientific_name_column(self, tmp_path: Path) -> None:
        """Test deduplication when scientific_name column is missing."""
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR001", "SRR002"],
                "sample_group": ["group1", "group1", "group2"],
            }
        )

        input_path = tmp_path / "metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path)

        assert result is True

        # Should still deduplicate by run
        df_result = pd.read_csv(input_path, sep="\t")
        assert len(df_result) == 2

    def test_deduplicate_preserves_all_columns(self, tmp_path: Path) -> None:
        """Test that all metadata columns are preserved."""
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR001", "SRR002"],
                "scientific_name": ["Homo sapiens", "Homo sapiens", "Mus musculus"],
                "sample_group": ["brain", "brain", "liver"],
                "library_strategy": ["RNA-Seq", "RNA-Seq", "RNA-Seq"],
                "platform": ["ILLUMINA", "ILLUMINA", "ILLUMINA"],
                "extra_col": ["a", "b", "c"],
            }
        )

        input_path = tmp_path / "metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path)

        assert result is True

        df_result = pd.read_csv(input_path, sep="\t")
        assert len(df_result) == 2
        assert list(df_result.columns) == [
            "run",
            "scientific_name",
            "sample_group",
            "library_strategy",
            "platform",
            "extra_col",
        ]

    def test_deduplicate_with_many_placeholders(self, tmp_path: Path) -> None:
        """Test handling when many entries have placeholder names."""
        placeholder = "Please add in format: Genus species"
        metadata = pd.DataFrame(
            {
                "run": ["SRR001", "SRR002", "SRR003", "SRR004", "SRR005"],
                "scientific_name": [placeholder, "Real species", placeholder, placeholder, "Another species"],
                "sample_group": ["g1", "g2", "g3", "g4", "g5"],
            }
        )

        input_path = tmp_path / "metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path)

        assert result is True

        df_result = pd.read_csv(input_path, sep="\t")
        # All runs should be unique, no deduplication needed
        assert len(df_result) == 5

    def test_deduplicate_large_dataset(self, tmp_path: Path) -> None:
        """Test performance with larger dataset."""
        # Create 1000 samples with some duplicates
        runs = [f"SRR{i:06d}" for i in range(500)] * 2
        names = ["Species " + chr(65 + (i % 26)) for i in range(1000)]

        metadata = pd.DataFrame(
            {
                "run": runs,
                "scientific_name": names,
            }
        )

        input_path = tmp_path / "large_metadata.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path)

        assert result is True

        df_result = pd.read_csv(input_path, sep="\t")
        # Should have exactly 500 unique runs
        assert len(df_result) == 500


class TestMetadataUtilsIntegration:
    """Integration tests for metadata utilities."""

    def test_deduplicate_real_world_pattern(self, tmp_path: Path) -> None:
        """Test with realistic amalgkit metadata patterns."""
        # Simulate amalgkit integrate output with duplicates
        placeholder = "Please add in format: Genus species"
        metadata = pd.DataFrame(
            {
                "run": ["SRR14740487", "SRR14740487", "SRR14740488", "SRR14740488", "SRR14740489"],
                "scientific_name": [
                    "Pogonomyrmex barbatus",
                    placeholder,
                    placeholder,
                    "Pogonomyrmex barbatus",
                    "Pogonomyrmex barbatus",
                ],
                "sample_group": ["worker", "worker", "queen", "queen", "male"],
                "library_strategy": ["RNA-Seq"] * 5,
                "library_layout": ["PAIRED"] * 5,
                "instrument_model": ["Illumina HiSeq 2500"] * 5,
                "bio_material": ["caste:worker", "caste:worker", "caste:queen", "caste:queen", "caste:male"],
            }
        )

        input_path = tmp_path / "metadata_integrate.tsv"
        metadata.to_csv(input_path, sep="\t", index=False)

        result = deduplicate_metadata(input_path)

        assert result is True

        df_result = pd.read_csv(input_path, sep="\t")

        # Should have 3 unique runs
        assert len(df_result) == 3
        assert set(df_result["run"]) == {"SRR14740487", "SRR14740488", "SRR14740489"}

        # Real names should be preserved over placeholders
        for _, row in df_result.iterrows():
            assert row["scientific_name"] != placeholder or row["run"] not in ["SRR14740487", "SRR14740488"]
