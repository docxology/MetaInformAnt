"""Tests for single-cell visualization functions.

Real implementation testing for UMAP, tSNE, PCA, trajectory, marker
expression, QC, and cluster comparison plots.
No mocking used - all tests use real data and matplotlib Agg backend.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for tests

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from metainformant.singlecell.data.preprocessing import SingleCellData

# Import visualization functions
try:
    from metainformant.singlecell.visualization.visualization import (
        plot_cluster_comparison,
        plot_marker_expression,
        plot_pca,
        plot_qc_metrics,
        plot_trajectory,
        plot_tsne,
        plot_umap,
    )

    VIZ_AVAILABLE = True
except ImportError:
    VIZ_AVAILABLE = False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_singlecell_data(
    n_cells: int = 60,
    n_genes: int = 40,
    seed: int = 42,
    add_umap: bool = False,
    add_tsne: bool = False,
    add_pca: bool = False,
    add_clusters: bool = False,
    add_qc: bool = False,
) -> SingleCellData:
    """Build a SingleCellData object with optional embedding columns.

    Uses ``pd.option_context('future.infer_string', False)`` so that string
    columns are stored with the legacy ``object`` dtype.  The visualization
    source code calls ``np.issubdtype(colors.dtype, np.number)`` which is
    incompatible with pandas 3.x ``StringDtype``.
    """
    rng = np.random.RandomState(seed)
    X = rng.exponential(1.0, size=(n_cells, n_genes)).astype(np.float64)

    gene_names = [f"gene_{i}" for i in range(n_genes)]
    cell_names = [f"cell_{i}" for i in range(n_cells)]

    with pd.option_context("future.infer_string", False):
        obs = pd.DataFrame(index=cell_names)
        var = pd.DataFrame(index=gene_names)

        if add_umap:
            obs["UMAP_1"] = rng.randn(n_cells)
            obs["UMAP_2"] = rng.randn(n_cells)

        if add_tsne:
            obs["tSNE_1"] = rng.randn(n_cells)
            obs["tSNE_2"] = rng.randn(n_cells)

        if add_pca:
            for i in range(3):
                obs[f"PC{i+1}"] = rng.randn(n_cells)

        if add_clusters:
            obs["cluster"] = [f"c{i % 3}" for i in range(n_cells)]
            obs["leiden"] = [f"L{i % 4}" for i in range(n_cells)]

        if add_qc:
            obs["n_counts"] = rng.randint(500, 5000, size=n_cells).astype(float)
            obs["n_genes"] = rng.randint(100, 1000, size=n_cells).astype(float)
            obs["pct_mito"] = rng.uniform(0, 15, size=n_cells)

    return SingleCellData(X=X, obs=obs, var=var, uns={})


@pytest.fixture(autouse=True)
def _close_plots():
    """Close all matplotlib figures after each test to free memory."""
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_umap
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotUMAP:
    """Tests for plot_umap function."""

    def test_basic_umap(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        fig = plot_umap(data)
        assert fig is not None
        assert isinstance(fig, plt.Figure)

    def test_umap_with_categorical_color(self) -> None:
        data = _make_singlecell_data(add_umap=True, add_clusters=True)
        fig = plot_umap(data, color="cluster")
        assert fig is not None

    def test_umap_with_numeric_color(self) -> None:
        data = _make_singlecell_data(add_umap=True, add_qc=True)
        fig = plot_umap(data, color="n_counts")
        assert fig is not None

    def test_umap_missing_coordinates_raises(self) -> None:
        data = _make_singlecell_data()  # no UMAP
        with pytest.raises(Exception):
            plot_umap(data)

    def test_umap_with_kwargs(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        fig = plot_umap(data, alpha=0.3, s=5, figsize=(10, 8))
        assert fig is not None

    def test_umap_invalid_input_type(self) -> None:
        with pytest.raises(Exception):
            plot_umap("not_a_singlecelldata")


# ---------------------------------------------------------------------------
# plot_tsne
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotTSNE:
    """Tests for plot_tsne function."""

    def test_basic_tsne(self) -> None:
        data = _make_singlecell_data(add_tsne=True)
        fig = plot_tsne(data)
        assert fig is not None
        assert isinstance(fig, plt.Figure)

    def test_tsne_with_categorical_color(self) -> None:
        data = _make_singlecell_data(add_tsne=True, add_clusters=True)
        fig = plot_tsne(data, color="cluster")
        assert fig is not None

    def test_tsne_with_numeric_color(self) -> None:
        data = _make_singlecell_data(add_tsne=True, add_qc=True)
        fig = plot_tsne(data, color="n_counts")
        assert fig is not None

    def test_tsne_missing_coordinates_raises(self) -> None:
        data = _make_singlecell_data()
        with pytest.raises(Exception):
            plot_tsne(data)


# ---------------------------------------------------------------------------
# plot_pca
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotPCA:
    """Tests for plot_pca function."""

    def test_basic_pca_2d(self) -> None:
        data = _make_singlecell_data(add_pca=True)
        fig = plot_pca(data, n_components=2)
        assert fig is not None
        assert isinstance(fig, plt.Figure)

    def test_pca_3d(self) -> None:
        data = _make_singlecell_data(add_pca=True)
        fig = plot_pca(data, n_components=3)
        assert fig is not None

    def test_pca_with_color(self) -> None:
        data = _make_singlecell_data(add_pca=True, add_clusters=True)
        fig = plot_pca(data, color="cluster", n_components=2)
        assert fig is not None

    def test_pca_with_numeric_color_3d(self) -> None:
        data = _make_singlecell_data(add_pca=True, add_qc=True)
        fig = plot_pca(data, color="n_counts", n_components=3)
        assert fig is not None

    def test_pca_missing_coordinates_raises(self) -> None:
        data = _make_singlecell_data()
        with pytest.raises(Exception):
            plot_pca(data, n_components=2)

    def test_pca_invalid_n_components_raises(self) -> None:
        data = _make_singlecell_data(add_pca=True)
        with pytest.raises(Exception):
            plot_pca(data, n_components=4)


# ---------------------------------------------------------------------------
# plot_trajectory
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotTrajectory:
    """Tests for plot_trajectory function."""

    def test_dpt_trajectory(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        # Add pseudotime and DPT metadata
        data.obs["dpt_pseudotime"] = np.linspace(0, 1, data.n_obs)
        data.uns["dpt"] = {"method": "diffusion_pseudotime", "root_cell": 0}
        fig = plot_trajectory(data, trajectory_key="dpt")
        assert fig is not None

    def test_trajectory_with_root_cell_marker(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        data.obs["dpt_pseudotime"] = np.linspace(0, 1, data.n_obs)
        data.uns["dpt"] = {
            "method": "diffusion_pseudotime",
            "root_cell": 5,
        }
        fig = plot_trajectory(data, trajectory_key="dpt")
        assert fig is not None

    def test_paga_trajectory(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        data.uns["paga"] = {"groups": "cluster", "connectivity_matrix": []}
        fig = plot_trajectory(data, trajectory_key="paga", color_by_pseudotime=False)
        assert fig is not None

    def test_trajectory_missing_key_raises(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        with pytest.raises(Exception):
            plot_trajectory(data, trajectory_key="nonexistent")

    def test_trajectory_no_embedding_raises(self) -> None:
        data = _make_singlecell_data()
        data.uns["dpt"] = {"root_cell": 0}
        with pytest.raises(Exception):
            plot_trajectory(data, trajectory_key="dpt")


# ---------------------------------------------------------------------------
# plot_marker_expression
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotMarkerExpression:
    """Tests for plot_marker_expression function."""

    def test_dotplot(self) -> None:
        data = _make_singlecell_data()
        fig = plot_marker_expression(data, marker_genes=["gene_0", "gene_1", "gene_2"], method="dotplot")
        assert fig is not None

    def test_heatmap(self) -> None:
        data = _make_singlecell_data()
        fig = plot_marker_expression(
            data,
            marker_genes=["gene_0", "gene_1"],
            method="heatmap",
        )
        assert fig is not None

    def test_violin(self) -> None:
        data = _make_singlecell_data()
        fig = plot_marker_expression(
            data,
            marker_genes=["gene_0", "gene_1", "gene_2"],
            method="violin",
        )
        assert fig is not None

    def test_missing_genes_raises(self) -> None:
        data = _make_singlecell_data()
        with pytest.raises(Exception):
            plot_marker_expression(data, marker_genes=["NOT_A_GENE"], method="dotplot")

    def test_invalid_method_raises(self) -> None:
        data = _make_singlecell_data()
        with pytest.raises(Exception):
            plot_marker_expression(data, marker_genes=["gene_0"], method="invalid_method")


# ---------------------------------------------------------------------------
# plot_qc_metrics
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotQCMetrics:
    """Tests for plot_qc_metrics function."""

    def test_basic_qc_plot(self) -> None:
        data = _make_singlecell_data(add_qc=True)
        fig = plot_qc_metrics(data)
        assert fig is not None
        assert isinstance(fig, plt.Figure)

    def test_qc_without_metrics(self) -> None:
        """Should still produce a figure (possibly empty subplots)."""
        data = _make_singlecell_data()
        fig = plot_qc_metrics(data)
        assert fig is not None

    def test_qc_with_kwargs(self) -> None:
        data = _make_singlecell_data(add_qc=True)
        fig = plot_qc_metrics(data, figsize=(14, 12))
        assert fig is not None


# ---------------------------------------------------------------------------
# plot_cluster_comparison
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not VIZ_AVAILABLE, reason="visualization not available")
class TestPlotClusterComparison:
    """Tests for plot_cluster_comparison function."""

    def test_single_clustering(self) -> None:
        data = _make_singlecell_data(add_umap=True, add_clusters=True)
        fig = plot_cluster_comparison(data, cluster_cols=["cluster"])
        assert fig is not None

    def test_multiple_clusterings(self) -> None:
        data = _make_singlecell_data(add_umap=True, add_clusters=True)
        fig = plot_cluster_comparison(data, cluster_cols=["cluster", "leiden"])
        assert fig is not None

    def test_explicit_embedding_cols(self) -> None:
        data = _make_singlecell_data(add_umap=True, add_clusters=True)
        fig = plot_cluster_comparison(
            data,
            cluster_cols=["cluster"],
            embedding_cols=["UMAP_1", "UMAP_2"],
        )
        assert fig is not None

    def test_missing_cluster_col_raises(self) -> None:
        data = _make_singlecell_data(add_umap=True)
        with pytest.raises(Exception):
            plot_cluster_comparison(data, cluster_cols=["nonexistent"])

    def test_no_embedding_raises(self) -> None:
        data = _make_singlecell_data(add_clusters=True)
        with pytest.raises(Exception):
            plot_cluster_comparison(data, cluster_cols=["cluster"])
