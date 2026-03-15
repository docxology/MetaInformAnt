"""Tests for metainformant.spatial.visualization -- all plot functions.

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from scipy import sparse as sp_sparse

from metainformant.spatial.visualization.plots import (
    plot_cell_type_map,
    plot_deconvolution_pie,
    plot_gene_expression_map,
    plot_interaction_heatmap,
    plot_neighborhood_graph,
    plot_spatial_autocorrelation,
    plot_spatial_scatter,
    plot_tissue_overlay,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def coordinates() -> np.ndarray:
    rng = np.random.RandomState(42)
    return rng.uniform(0, 100, (50, 2))


@pytest.fixture()
def numeric_values() -> np.ndarray:
    rng = np.random.RandomState(42)
    return rng.standard_normal(50)


@pytest.fixture()
def categorical_values() -> np.ndarray:
    rng = np.random.RandomState(42)
    return rng.choice(["A", "B", "C"], size=50)


@pytest.fixture()
def tissue_image() -> np.ndarray:
    """Small fake tissue image (100x100 RGB)."""
    rng = np.random.RandomState(42)
    return rng.randint(0, 255, (100, 100, 3), dtype=np.uint8)


@pytest.fixture()
def spatial_graph(coordinates: np.ndarray) -> sp_sparse.csr_matrix:
    from metainformant.spatial.analysis.clustering import build_spatial_graph

    return build_spatial_graph(coordinates, method="knn", n_neighbors=5)


@pytest.fixture()
def spatial_dataset_obj() -> SimpleNamespace:
    """Minimal spatial dataset-like object for plot_gene_expression_map / plot_cell_type_map."""
    rng = np.random.RandomState(42)
    n_spots, n_genes = 30, 5
    return SimpleNamespace(
        expression=np.abs(rng.standard_normal((n_spots, n_genes))),
        coordinates=rng.uniform(0, 100, (n_spots, 2)),
        gene_names=["CD3E", "CD8A", "MS4A1", "COL1A1", "EPCAM"],
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestPlotSpatialScatter:
    def test_numeric_scatter_creates_png(
        self, coordinates: np.ndarray, numeric_values: np.ndarray, tmp_path: Path
    ) -> None:
        out = tmp_path / "scatter_num.png"
        fig = plot_spatial_scatter(coordinates, numeric_values, out)
        assert out.exists()
        assert out.stat().st_size > 0
        assert fig is not None

    def test_categorical_scatter(self, coordinates: np.ndarray, categorical_values: np.ndarray, tmp_path: Path) -> None:
        out = tmp_path / "scatter_cat.png"
        fig = plot_spatial_scatter(coordinates, categorical_values, out)
        assert out.exists()
        assert fig is not None


class TestPlotTissueOverlay:
    def test_creates_output(
        self,
        coordinates: np.ndarray,
        numeric_values: np.ndarray,
        tissue_image: np.ndarray,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "overlay.png"
        fig = plot_tissue_overlay(coordinates, numeric_values, tissue_image, out)
        assert out.exists()
        assert out.stat().st_size > 0
        assert fig is not None


class TestPlotGeneExpressionMap:
    def test_creates_png(self, spatial_dataset_obj: SimpleNamespace, tmp_path: Path) -> None:
        out = tmp_path / "gene_expr.png"
        fig = plot_gene_expression_map(spatial_dataset_obj, "CD3E", out)
        assert out.exists()
        assert fig is not None

    def test_missing_gene_raises(self, spatial_dataset_obj: SimpleNamespace, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="not found"):
            plot_gene_expression_map(spatial_dataset_obj, "NONEXISTENT", tmp_path / "x.png")


class TestPlotCellTypeMap:
    def test_creates_png(self, spatial_dataset_obj: SimpleNamespace, tmp_path: Path) -> None:
        rng = np.random.RandomState(42)
        labels = rng.choice(["T-cell", "B-cell", "Fibroblast"], size=30)
        out = tmp_path / "cell_type.png"
        fig = plot_cell_type_map(spatial_dataset_obj, labels, out)
        assert out.exists()
        assert fig is not None


class TestPlotDeconvolutionPie:
    def test_creates_png(self, coordinates: np.ndarray, tmp_path: Path) -> None:
        rng = np.random.RandomState(42)
        fractions = rng.dirichlet([1, 1, 1], size=50)
        out = tmp_path / "deconv_pie.png"
        fig = plot_deconvolution_pie(
            coordinates,
            fractions,
            out,
            cell_type_names=["A", "B", "C"],
            pie_radius=3.0,
        )
        assert out.exists()
        assert fig is not None


class TestPlotSpatialAutocorrelation:
    def test_lisa_cluster_map(self, coordinates: np.ndarray, tmp_path: Path) -> None:
        rng = np.random.RandomState(42)
        scores = rng.standard_normal(50)
        labels = rng.choice(["HH", "LL", "HL", "LH", "NS"], size=50)
        out = tmp_path / "lisa.png"
        fig = plot_spatial_autocorrelation(coordinates, scores, out, cluster_labels=list(labels))
        assert out.exists()
        assert fig is not None

    def test_continuous_map(self, coordinates: np.ndarray, numeric_values: np.ndarray, tmp_path: Path) -> None:
        out = tmp_path / "autocorr_continuous.png"
        fig = plot_spatial_autocorrelation(coordinates, numeric_values, out)
        assert out.exists()
        assert fig is not None


class TestPlotInteractionHeatmap:
    def test_creates_png(self, tmp_path: Path) -> None:
        rng = np.random.RandomState(42)
        matrix = rng.standard_normal((3, 3))
        out = tmp_path / "interaction.png"
        fig = plot_interaction_heatmap(matrix, out, cell_type_names=["A", "B", "C"])
        assert out.exists()
        assert fig is not None


class TestPlotNeighborhoodGraph:
    def test_creates_png(
        self,
        coordinates: np.ndarray,
        spatial_graph: sp_sparse.csr_matrix,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "neighborhood.png"
        fig = plot_neighborhood_graph(coordinates, spatial_graph, out)
        assert out.exists()
        assert fig is not None

    def test_with_node_colors(
        self,
        coordinates: np.ndarray,
        spatial_graph: sp_sparse.csr_matrix,
        numeric_values: np.ndarray,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "neighborhood_colored.png"
        fig = plot_neighborhood_graph(coordinates, spatial_graph, out, node_colors=numeric_values)
        assert out.exists()
        assert fig is not None
