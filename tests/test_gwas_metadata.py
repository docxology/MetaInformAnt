"""Tests for GWAS sample metadata loading, validation, and merging."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.data.metadata import (
    get_geographic_coordinates,
    get_population_labels,
    load_sample_metadata,
    merge_metadata_with_phenotypes,
    validate_metadata,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_tsv(path: Path, lines: list[str]) -> Path:
    """Write tab-separated lines to *path* and return it."""
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
# load_sample_metadata
# ---------------------------------------------------------------------------


class TestLoadSampleMetadata:
    """Tests for loading TSV metadata files."""

    def test_load_basic(self, tmp_path: Path) -> None:
        """Load a well-formed TSV with several columns."""
        tsv = _write_tsv(
            tmp_path / "meta.tsv",
            [
                "sample_id\tsubspecies\tpopulation\tlocation",
                "S001\tligustica\tIT_North\tMilan",
                "S002\tcarnica\tAT_East\tVienna",
                "S003\tmellifera\tFR_West\tParis",
            ],
        )
        result = load_sample_metadata(tsv)
        assert result["status"] == "success"
        assert result["n_samples"] == 3
        assert set(result["columns"]) == {"subspecies", "population", "location"}
        assert result["metadata"]["S001"]["subspecies"] == "ligustica"
        assert result["metadata"]["S002"]["location"] == "Vienna"

    def test_load_missing_file(self, tmp_path: Path) -> None:
        """Return error status when file does not exist."""
        result = load_sample_metadata(tmp_path / "nonexistent.tsv")
        assert result["status"] == "error"
        assert result["n_samples"] == 0

    def test_load_empty_file(self, tmp_path: Path) -> None:
        """Return error status for a completely empty file."""
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        result = load_sample_metadata(empty)
        assert result["status"] == "error"
        assert result["n_samples"] == 0

    def test_load_missing_sample_id_column(self, tmp_path: Path) -> None:
        """Return error when header lacks sample_id."""
        tsv = _write_tsv(
            tmp_path / "bad_header.tsv",
            ["name\tpopulation", "S001\tIT_North"],
        )
        result = load_sample_metadata(tsv)
        assert result["status"] == "error"
        assert "sample_id" in result["message"]

    def test_load_skips_blank_sample_id(self, tmp_path: Path) -> None:
        """Rows with empty sample_id should be silently skipped."""
        tsv = _write_tsv(
            tmp_path / "blanks.tsv",
            [
                "sample_id\tpopulation",
                "S001\tIT_North",
                "\tGhost",
                "S002\tAT_East",
            ],
        )
        result = load_sample_metadata(tsv)
        assert result["status"] == "success"
        assert result["n_samples"] == 2
        assert "S001" in result["metadata"]
        assert "S002" in result["metadata"]

    def test_load_header_only(self, tmp_path: Path) -> None:
        """A file with only a header row should yield zero samples."""
        tsv = _write_tsv(tmp_path / "header_only.tsv", ["sample_id\tpopulation"])
        result = load_sample_metadata(tsv)
        assert result["status"] == "success"
        assert result["n_samples"] == 0


# ---------------------------------------------------------------------------
# merge_metadata_with_phenotypes
# ---------------------------------------------------------------------------


class TestMergeMetadataWithPhenotypes:
    """Tests for merging metadata with phenotype data."""

    def test_full_merge(self) -> None:
        """All samples have phenotype data."""
        metadata = {
            "S001": {"population": "IT_North"},
            "S002": {"population": "AT_East"},
        }
        phenotypes = {"S001": 1.5, "S002": 2.3}
        result = merge_metadata_with_phenotypes(metadata, phenotypes)
        assert result["status"] == "success"
        assert result["n_matched"] == 2
        assert result["n_unmatched"] == 0
        assert result["merged"]["S001"]["phenotype"] == 1.5

    def test_partial_merge(self) -> None:
        """Some samples lack phenotype data."""
        metadata = {
            "S001": {"population": "IT_North"},
            "S002": {"population": "AT_East"},
        }
        phenotypes = {"S001": 1.5}
        result = merge_metadata_with_phenotypes(metadata, phenotypes)
        assert result["n_matched"] == 1
        assert result["n_unmatched"] == 1
        assert result["merged"]["S002"]["phenotype"] is None

    def test_empty_inputs(self) -> None:
        """Both inputs empty should succeed with zero counts."""
        result = merge_metadata_with_phenotypes({}, {})
        assert result["status"] == "success"
        assert result["n_matched"] == 0


# ---------------------------------------------------------------------------
# validate_metadata
# ---------------------------------------------------------------------------


class TestValidateMetadata:
    """Tests for metadata completeness validation."""

    def test_complete(self) -> None:
        """All expected samples present."""
        metadata = {"S001": {}, "S002": {}, "S003": {}}
        result = validate_metadata(metadata, ["S001", "S002", "S003"])
        assert result["valid"] is True
        assert result["completeness"] == pytest.approx(1.0)
        assert result["missing_samples"] == []
        assert result["extra_samples"] == []

    def test_missing_samples(self) -> None:
        """Some expected samples not in metadata."""
        metadata = {"S001": {}}
        result = validate_metadata(metadata, ["S001", "S002", "S003"])
        assert result["valid"] is False
        assert set(result["missing_samples"]) == {"S002", "S003"}
        assert result["completeness"] == pytest.approx(1 / 3)

    def test_extra_samples(self) -> None:
        """Metadata has samples not in expected list."""
        metadata = {"S001": {}, "S002": {}, "S099": {}}
        result = validate_metadata(metadata, ["S001", "S002"])
        assert result["valid"] is True
        assert result["extra_samples"] == ["S099"]
        assert result["completeness"] == pytest.approx(1.0)

    def test_empty_expected(self) -> None:
        """Empty expected list: completeness is 0.0 (avoid division by zero)."""
        result = validate_metadata({"S001": {}}, [])
        assert result["valid"] is True
        assert result["completeness"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# get_population_labels
# ---------------------------------------------------------------------------


class TestGetPopulationLabels:
    """Tests for population label extraction."""

    def test_basic(self) -> None:
        metadata = {
            "S001": {"population": "IT_North", "location": "Milan"},
            "S002": {"population": "AT_East", "location": "Vienna"},
        }
        labels = get_population_labels(metadata)
        assert labels == {"S001": "IT_North", "S002": "AT_East"}

    def test_custom_column(self) -> None:
        metadata = {"S001": {"subspecies": "ligustica"}}
        labels = get_population_labels(metadata, column="subspecies")
        assert labels == {"S001": "ligustica"}

    def test_missing_column(self) -> None:
        """Samples without the requested column are omitted."""
        metadata = {"S001": {"population": "IT_North"}, "S002": {"location": "Vienna"}}
        labels = get_population_labels(metadata)
        assert "S001" in labels
        assert "S002" not in labels

    def test_empty_value_skipped(self) -> None:
        metadata = {"S001": {"population": ""}, "S002": {"population": "AT_East"}}
        labels = get_population_labels(metadata)
        assert labels == {"S002": "AT_East"}


# ---------------------------------------------------------------------------
# get_geographic_coordinates
# ---------------------------------------------------------------------------


class TestGetGeographicCoordinates:
    """Tests for geographic coordinate extraction."""

    def test_basic(self) -> None:
        metadata = {
            "S001": {"latitude": "45.46", "longitude": "9.19"},
            "S002": {"latitude": "48.21", "longitude": "16.37"},
        }
        coords = get_geographic_coordinates(metadata)
        assert len(coords) == 2
        assert coords[0]["sample_id"] == "S001"
        assert coords[0]["latitude"] == pytest.approx(45.46)
        assert coords[0]["longitude"] == pytest.approx(9.19)

    def test_non_numeric_skipped(self) -> None:
        """Non-numeric lat/lon are silently skipped."""
        metadata = {
            "S001": {"latitude": "not_a_number", "longitude": "9.19"},
            "S002": {"latitude": "48.21", "longitude": "16.37"},
        }
        coords = get_geographic_coordinates(metadata)
        assert len(coords) == 1
        assert coords[0]["sample_id"] == "S002"

    def test_missing_columns(self) -> None:
        """Samples without lat/lon columns produce no coordinates."""
        metadata = {"S001": {"population": "IT_North"}}
        coords = get_geographic_coordinates(metadata)
        assert coords == []

    def test_partial_coordinates_skipped(self) -> None:
        """Samples with only latitude or only longitude are skipped."""
        metadata = {
            "S001": {"latitude": "45.0", "longitude": ""},
            "S002": {"latitude": "", "longitude": "9.0"},
        }
        coords = get_geographic_coordinates(metadata)
        assert coords == []

    def test_roundtrip_through_load(self, tmp_path: Path) -> None:
        """Coordinates loaded from TSV can be extracted correctly."""
        tsv = _write_tsv(
            tmp_path / "geo.tsv",
            [
                "sample_id\tlatitude\tlongitude",
                "S001\t45.46\t9.19",
                "S002\t48.21\t16.37",
            ],
        )
        loaded = load_sample_metadata(tsv)
        assert loaded["status"] == "success"
        coords = get_geographic_coordinates(loaded["metadata"])
        assert len(coords) == 2
        assert coords[0]["latitude"] == pytest.approx(45.46)
