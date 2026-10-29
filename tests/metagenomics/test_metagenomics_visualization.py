"""Tests for metagenomics visualization submodule.

Tests all plot functions output valid files.
Uses real implementations -- real-implementation policy.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.metagenomics.diversity.metrics import ordination
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
        plot_stacked_bar(abundance_data, output_path=out)
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
        plot_rarefaction_curves(otu_table, step=20, output_path=out, seed=42)
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
        plot_ordination(distance_matrix_dict, output_path=out)
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
        plot_alpha_diversity(otu_table, output_path=out)
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
        plot_heatmap(abundance_data, output_path=out)
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


# ---------------------------------------------------------------------------
# Regression tests: plot_ordination(method=) selects the algorithm
# ---------------------------------------------------------------------------


def _pairwise(points: np.ndarray) -> np.ndarray:
    return np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)


def _plotted_points(fig: object, n_samples: int) -> np.ndarray:
    ax = fig.axes[0]  # type: ignore[attr-defined]
    offsets = [np.asarray(c.get_offsets()).reshape(-1, 2) for c in ax.collections]
    points = np.vstack(offsets)
    assert points.shape == (n_samples, 2)
    return points


class TestPlotOrdinationMethod:
    """`method` used to change only the title; the plot was always PCoA."""

    DM = {
        "s1": {"s1": 0.0, "s2": 0.2, "s3": 0.9, "s4": 0.4, "s5": 0.7},
        "s2": {"s1": 0.2, "s2": 0.0, "s3": 0.6, "s4": 0.8, "s5": 0.3},
        "s3": {"s1": 0.9, "s2": 0.6, "s3": 0.0, "s4": 0.5, "s5": 0.95},
        "s4": {"s1": 0.4, "s2": 0.8, "s3": 0.5, "s4": 0.0, "s5": 0.35},
        "s5": {"s1": 0.7, "s2": 0.3, "s3": 0.95, "s4": 0.35, "s5": 0.0},
    }

    @classmethod
    def _matrix(cls) -> list[list[float]]:
        samples = sorted(cls.DM)
        return [[cls.DM[a][b] for b in samples] for a in samples]

    def test_nmds_plots_the_nmds_configuration(self) -> None:
        fig = plot_ordination(self.DM, method="NMDS", seed=7)
        try:
            plotted = _plotted_points(fig, 5)
            expected = np.array(ordination(self._matrix(), method="nmds", n_components=2, seed=7)["coordinates"])
            # Ordinations are defined up to rotation/reflection, so compare pairwise distances
            np.testing.assert_allclose(_pairwise(plotted), _pairwise(expected), atol=1e-9)
            assert fig.axes[0].get_xlabel() == "NMDS1"
            assert fig.axes[0].get_ylabel() == "NMDS2"
            assert "NMDS" in fig.axes[0].get_title()
        finally:
            plt.close(fig)

    def test_pcoa_plots_the_pcoa_configuration(self) -> None:
        fig = plot_ordination(self.DM, method="PCoA")
        try:
            plotted = _plotted_points(fig, 5)
            expected = np.array(ordination(self._matrix(), method="pcoa", n_components=2)["coordinates"])
            np.testing.assert_allclose(_pairwise(plotted), _pairwise(expected), atol=1e-5)
            assert fig.axes[0].get_xlabel().startswith("PC1 (")  # "PC1 (35.2%)": share of positive-eigenvalue variance
            assert fig.axes[0].get_xlabel().endswith("%)")
        finally:
            plt.close(fig)

    def test_methods_give_different_configurations(self) -> None:
        pcoa_fig = plot_ordination(self.DM, method="PCoA")
        nmds_fig = plot_ordination(self.DM, method="NMDS", seed=7)
        try:
            gap = np.abs(_pairwise(_plotted_points(pcoa_fig, 5)) - _pairwise(_plotted_points(nmds_fig, 5))).max()
            assert gap > 0.05
        finally:
            plt.close(pcoa_fig)
            plt.close(nmds_fig)

    def test_method_is_case_insensitive(self) -> None:
        lower = plot_ordination(self.DM, method="nmds", seed=3)
        upper = plot_ordination(self.DM, method="NMDS", seed=3)
        try:
            np.testing.assert_allclose(_plotted_points(lower, 5), _plotted_points(upper, 5))
        finally:
            plt.close(lower)
            plt.close(upper)

    def test_nmds_seed_is_reproducible_and_the_best_start_is_the_same_solution(self) -> None:
        """The seed draws the random starts; keeping the best of several makes the solution seed-independent."""
        figs = [plot_ordination(self.DM, method="NMDS", seed=s) for s in (1, 1, 2)]
        try:
            first, repeat, other = (_plotted_points(f, 5) for f in figs)
            np.testing.assert_allclose(first, repeat)
            np.testing.assert_allclose(_pairwise(first), _pairwise(other), atol=1e-6)
        finally:
            for f in figs:
                plt.close(f)
        starts = [ordination(self._matrix(), method="nmds", n_components=2, seed=s)["stress_by_start"] for s in (1, 2)]
        assert starts[0][1:] != starts[1][1:]  # the seeds really gave different random starts

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="method"):
            plot_ordination(self.DM, method="tsne")

    def test_unknown_method_raises_even_with_too_few_samples(self) -> None:
        with pytest.raises(ValueError, match="method"):
            plot_ordination({"s1": {"s1": 0.0}}, method="umap")

    def test_nmds_writes_a_file(self, tmp_path: Path) -> None:
        out = tmp_path / "nmds.png"
        assert plot_ordination(self.DM, method="NMDS", output_path=out) == str(out)
        assert out.stat().st_size > 0
