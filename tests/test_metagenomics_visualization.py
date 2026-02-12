"""Tests for metagenomics visualization submodule.

Tests all plot functions output valid files.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import pytest

from metainformant.metagenomics.visualization.plots import (
    plot_alpha_diversity,
    plot_heatmap,
    plot_krona_chart,
    plot_ordination,
    plot_rarefaction_curves,
    plot_stacked_bar,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def taxonomy_data() -> dict[str, float]:
    """Taxonomy lineage data for Krona chart."""
    return {
        "Bacteria;Firmicutes;Bacilli": 0.35,
        "Bacteria;Proteobacteria;Gammaproteobacteria": 0.25,
        "Bacteria;Bacteroidetes;Bacteroidia": 0.20,
        "Archaea;Euryarchaeota;Methanobacteria": 0.10,
        "Bacteria;Actinobacteria;Actinomycetia": 0.10,
    }


@pytest.fixture
def abundance_data() -> dict[str, dict[str, float]]:
    """Per-sample taxonomic abundance for stacked bar."""
    return {
        "sample_A": {"Firmicutes": 40, "Proteobacteria": 30, "Bacteroidetes": 20, "Actinobacteria": 10},
        "sample_B": {"Firmicutes": 10, "Proteobacteria": 50, "Bacteroidetes": 25, "Actinobacteria": 15},
        "sample_C": {"Firmicutes": 25, "Proteobacteria": 25, "Bacteroidetes": 25, "Actinobacteria": 25},
    }


@pytest.fixture
def otu_table() -> dict[str, dict[str, int]]:
    """OTU table for rarefaction curves and alpha diversity."""
    return {
        "sample_A": {"OTU1": 100, "OTU2": 50, "OTU3": 30, "OTU4": 10, "OTU5": 5},
        "sample_B": {"OTU1": 20, "OTU2": 80, "OTU3": 60, "OTU4": 30, "OTU5": 25},
    }


@pytest.fixture
def distance_matrix_dict() -> dict[str, dict[str, float]]:
    """Pairwise distance matrix as nested dict for ordination."""
    return {
        "s1": {"s1": 0.0, "s2": 0.3, "s3": 0.7, "s4": 0.9},
        "s2": {"s1": 0.3, "s2": 0.0, "s3": 0.5, "s4": 0.8},
        "s3": {"s1": 0.7, "s2": 0.5, "s3": 0.0, "s4": 0.4},
        "s4": {"s1": 0.9, "s2": 0.8, "s3": 0.4, "s4": 0.0},
    }


# ---------------------------------------------------------------------------
# Tests: plot_krona_chart
# ---------------------------------------------------------------------------


class TestPlotKronaChart:
    """Tests for Krona-style sunburst chart."""

    def test_creates_output_file(self, tmp_path: Path, taxonomy_data: dict[str, float]) -> None:
        out = tmp_path / "krona.png"
        result = plot_krona_chart(taxonomy_data, output_path=out)
        assert out.exists()
        assert out.stat().st_size > 0
        assert str(result) == str(out)

    def test_empty_data(self, tmp_path: Path) -> None:
        out = tmp_path / "krona_empty.png"
        plot_krona_chart({}, output_path=out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Tests: plot_stacked_bar
# ---------------------------------------------------------------------------


class TestPlotStackedBar:
    """Tests for stacked bar chart."""

    def test_creates_png(self, tmp_path: Path, abundance_data: dict[str, dict[str, float]]) -> None:
        out = tmp_path / "stacked.png"
        result = plot_stacked_bar(abundance_data, output_path=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_empty_data(self, tmp_path: Path) -> None:
        out = tmp_path / "stacked_empty.png"
        plot_stacked_bar({}, output_path=out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Tests: plot_rarefaction_curves
# ---------------------------------------------------------------------------


class TestPlotRarefactionCurves:
    """Tests for rarefaction curve plots."""

    def test_creates_png(self, tmp_path: Path, otu_table: dict[str, dict[str, int]]) -> None:
        out = tmp_path / "rarefaction.png"
        result = plot_rarefaction_curves(otu_table, step=20, output_path=out, seed=42)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_empty_data(self, tmp_path: Path) -> None:
        out = tmp_path / "rarefaction_empty.png"
        plot_rarefaction_curves({}, output_path=out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Tests: plot_ordination
# ---------------------------------------------------------------------------


class TestPlotOrdination:
    """Tests for ordination plot."""

    def test_creates_png(self, tmp_path: Path, distance_matrix_dict: dict[str, dict[str, float]]) -> None:
        out = tmp_path / "ordination.png"
        result = plot_ordination(distance_matrix_dict, output_path=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_with_groups(self, tmp_path: Path, distance_matrix_dict: dict[str, dict[str, float]]) -> None:
        out = tmp_path / "ordination_groups.png"
        groups = {"s1": "A", "s2": "A", "s3": "B", "s4": "B"}
        plot_ordination(distance_matrix_dict, output_path=out, groups=groups)
        assert out.exists()

    def test_too_few_samples(self, tmp_path: Path) -> None:
        out = tmp_path / "ordination_small.png"
        dm = {"s1": {"s1": 0.0, "s2": 0.5}, "s2": {"s1": 0.5, "s2": 0.0}}
        plot_ordination(dm, output_path=out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Tests: plot_alpha_diversity
# ---------------------------------------------------------------------------


class TestPlotAlphaDiversity:
    """Tests for alpha diversity boxplots."""

    def test_creates_png(self, tmp_path: Path, otu_table: dict[str, dict[str, int]]) -> None:
        out = tmp_path / "alpha.png"
        result = plot_alpha_diversity(otu_table, output_path=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_with_groups(self, tmp_path: Path, otu_table: dict[str, dict[str, int]]) -> None:
        out = tmp_path / "alpha_groups.png"
        groups = {"sample_A": "control", "sample_B": "treatment"}
        plot_alpha_diversity(otu_table, output_path=out, groups=groups)
        assert out.exists()

    def test_single_metric(self, tmp_path: Path, otu_table: dict[str, dict[str, int]]) -> None:
        out = tmp_path / "alpha_shannon.png"
        plot_alpha_diversity(otu_table, metrics=["shannon"], output_path=out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Tests: plot_heatmap
# ---------------------------------------------------------------------------


class TestPlotHeatmap:
    """Tests for abundance heatmap."""

    def test_creates_png(self, tmp_path: Path, abundance_data: dict[str, dict[str, float]]) -> None:
        out = tmp_path / "heatmap.png"
        result = plot_heatmap(abundance_data, output_path=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_no_cluster(self, tmp_path: Path, abundance_data: dict[str, dict[str, float]]) -> None:
        out = tmp_path / "heatmap_nocluster.png"
        plot_heatmap(abundance_data, output_path=out, cluster=False)
        assert out.exists()

    def test_empty_data(self, tmp_path: Path) -> None:
        out = tmp_path / "heatmap_empty.png"
        plot_heatmap({}, output_path=out)
        assert out.exists()
