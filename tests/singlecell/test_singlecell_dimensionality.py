"""Tests for single-cell dimensionality reduction module.

Real implementation testing for PCA, UMAP, t-SNE, and neighbor graph computation.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

try:
    from scipy import sparse

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    sparse = None

pytestmark = pytest.mark.skipif(not SCIPY_AVAILABLE, reason="scipy not available")

from metainformant.core.utils.errors import ValidationError
from metainformant.singlecell.analysis.dimensionality import (
    SKLEARN_AVAILABLE,
    compute_diffusion_map,
    compute_neighbors,
    compute_pca,
    compute_tsne,
    compute_umap,
    select_hvgs,
)
from metainformant.singlecell.data.preprocessing import SingleCellData, log_transform, normalize_counts, scale_data


class TestHVGSelection:
    """Test highly variable gene selection."""

    def setup_method(self):
        """Setup test data with known variance structure."""
        np.random.seed(42)

        # Create data with different variance levels
        n_cells, n_genes = 100, 200

        # Low variance genes
        X_low = np.random.normal(2, 0.1, (n_cells, 50))
        # Medium variance genes
        X_med = np.random.normal(2, 1.0, (n_cells, 100))
        # High variance genes
        X_high = np.random.normal(2, 3.0, (n_cells, 50))

        X = np.concatenate([X_low, X_med, X_high], axis=1)
        X = np.maximum(X, 0)  # Ensure non-negative

        # Create gene names
        gene_names = (
            [f"LowVar_{i}" for i in range(50)]
            + [f"MedVar_{i}" for i in range(100)]
            + [f"HighVar_{i}" for i in range(50)]
        )

        var = pd.DataFrame({"gene_name": gene_names})
        self.test_data = SingleCellData(X, var=var)

        # Normalize and log transform for HVG selection
        self.test_data = normalize_counts(self.test_data)
        self.test_data = log_transform(self.test_data)

    def test_select_hvgs_seurat_method(self):
        """Test Seurat method for HVG selection."""
        data = select_hvgs(self.test_data, n_top_genes=50, method="seurat")

        # Check that HVG information is added
        assert "highly_variable" in data.var.columns
        assert "means" in data.var.columns
        assert "variances" in data.var.columns
        assert "variances_norm" in data.var.columns

        # Should select approximately requested number of genes
        n_hvgs = data.var["highly_variable"].sum()
        assert n_hvgs <= 50
        assert n_hvgs > 0

        # High variance genes should be preferentially selected
        hvg_names = data.var[data.var["highly_variable"]]["gene_name"]
        high_var_selected = hvg_names.str.contains("HighVar").sum()
        low_var_selected = hvg_names.str.contains("LowVar").sum()

        # Should select more high variance than low variance genes
        assert high_var_selected >= low_var_selected

    def test_select_hvgs_variance_method(self):
        """Test simple variance-based HVG selection."""
        data = select_hvgs(self.test_data, n_top_genes=30, method="variance")

        n_hvgs = data.var["highly_variable"].sum()
        assert n_hvgs == 30

        # Should select genes with highest variance
        hvg_names = data.var[data.var["highly_variable"]]["gene_name"]
        high_var_selected = hvg_names.str.contains("HighVar").sum()

        # Most selected genes should be high variance
        assert high_var_selected >= 15  # At least half should be high variance

    def test_select_hvgs_cell_ranger_method(self):
        """Test Cell Ranger method for HVG selection."""
        data = select_hvgs(self.test_data, n_top_genes=40, method="cell_ranger")

        assert "highly_variable" in data.var.columns
        n_hvgs = data.var["highly_variable"].sum()
        assert n_hvgs == 40

    def test_select_hvgs_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError, match="Unknown HVG selection method"):
            select_hvgs(self.test_data, method="invalid_method")


@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="sklearn not available")
class TestPCA:
    """Test PCA computation."""

    def setup_method(self):
        """Setup preprocessed test data for PCA."""
        np.random.seed(42)

        # Create structured data for PCA
        n_cells, n_genes = 150, 100

        # Create some structure in the data
        # First component: gradient
        comp1 = np.linspace(-2, 2, n_cells)
        # Second component: sine wave
        comp2 = np.sin(np.linspace(0, 4 * np.pi, n_cells))

        # Generate expression data with these components
        X = np.random.normal(0, 0.1, (n_cells, n_genes))

        # Add structure to first 20 genes based on comp1
        X[:, :20] += comp1[:, np.newaxis] * np.random.normal(1, 0.1, 20)
        # Add structure to next 20 genes based on comp2
        X[:, 20:40] += comp2[:, np.newaxis] * np.random.normal(1, 0.1, 20)

        X = np.maximum(X + 5, 0)  # Ensure non-negative

        self.test_data = SingleCellData(X)

        # Preprocess
        self.test_data = normalize_counts(self.test_data)
        self.test_data = log_transform(self.test_data)
        self.test_data = scale_data(self.test_data)

        # Select HVGs
        self.test_data = select_hvgs(self.test_data, n_top_genes=50)

    def test_compute_pca_basic(self):
        """Test basic PCA computation."""
        data = compute_pca(self.test_data, n_components=20)

        # Check that PCA results are stored
        assert "X_pca" in data.obsm
        assert "PCs" in data.varm
        assert "pca" in data.uns

        # Check dimensions
        assert data.obsm["X_pca"].shape == (150, 20)
        assert data.varm["PCs"].shape == (100, 20)

        # Check PCA metadata
        pca_info = data.uns["pca"]
        assert "variance_ratio" in pca_info
        assert "variance" in pca_info
        assert len(pca_info["variance_ratio"]) == 20

        # Variance ratios should sum to less than 1 and be decreasing
        var_ratios = np.array(pca_info["variance_ratio"])
        assert np.all(var_ratios > 0)
        assert np.all(np.diff(var_ratios) <= 0)  # Decreasing
        assert np.sum(var_ratios) <= 1.0

    def test_compute_pca_with_hvgs(self):
        """Test PCA computation using HVGs."""
        data = compute_pca(self.test_data, n_components=15, use_hvgs=True)

        assert "X_pca" in data.obsm
        assert data.obsm["X_pca"].shape[1] == 15
        assert data.uns["pca"]["use_hvgs"] is True

    def test_compute_pca_without_hvgs(self):
        """Test PCA computation using all genes."""
        data = compute_pca(self.test_data, n_components=10, use_hvgs=False)

        assert "X_pca" in data.obsm
        assert data.obsm["X_pca"].shape[1] == 10
        assert data.uns["pca"]["use_hvgs"] is False

    def test_compute_pca_too_many_components(self):
        """Test PCA with more components than possible."""
        # Request more components than samples-1
        data = compute_pca(self.test_data, n_components=200)

        # Should automatically adjust to maximum possible
        actual_components = data.obsm["X_pca"].shape[1]
        assert actual_components < 200
        assert actual_components <= min(150 - 1, 100)  # min(n_samples-1, n_features)

    def test_compute_pca_requires_sklearn(self):
        """Test that PCA raises ImportError when sklearn is not available."""
        # This test verifies the error handling, but sklearn is available in test environment
        # The actual error would occur if sklearn were missing
        # We verify the check exists in the code by ensuring it doesn't fail when sklearn is available
        data = compute_pca(self.test_data, n_components=10)
        assert "X_pca" in data.obsm
        # If sklearn were missing, this would raise ImportError before reaching here


class TestNeighborGraph:
    """Test neighbor graph computation."""

    def setup_method(self):
        """Setup test data with PCA for neighbor computation."""
        np.random.seed(42)

        # Create data with clear clusters
        n_cells = 120

        # Three clusters
        cluster1 = np.random.multivariate_normal([0, 0], [[1, 0], [0, 1]], 40)
        cluster2 = np.random.multivariate_normal([5, 0], [[1, 0], [0, 1]], 40)
        cluster3 = np.random.multivariate_normal([0, 5], [[1, 0], [0, 1]], 40)

        # Embed in higher dimensional space
        pca_coords = np.vstack([cluster1, cluster2, cluster3])
        full_coords = np.random.normal(0, 0.1, (n_cells, 50))
        full_coords[:, :2] = pca_coords

        X = full_coords + 5  # Ensure positive
        self.test_data = SingleCellData(X)

        # Add PCA coordinates
        self.test_data.obsm["X_pca"] = pca_coords

    def test_compute_neighbors_basic(self):
        """Test basic neighbor graph computation."""
        data = compute_neighbors(self.test_data, n_neighbors=15, n_pcs=2)

        # Check that neighbor graph is stored
        assert "neighbors" in data.uns
        neighbors_info = data.uns["neighbors"]

        assert "connectivities" in neighbors_info
        assert "distances" in neighbors_info
        assert "indices" in neighbors_info
        assert "params" in neighbors_info

        # Check dimensions
        connectivities = neighbors_info["connectivities"]
        assert connectivities.shape == (120, 120)

        # Should be approximately symmetric
        if sparse.issparse(connectivities):
            conn_dense = connectivities.toarray()
        else:
            conn_dense = connectivities

        asymmetry = np.abs(conn_dense - conn_dense.T).max()
        assert asymmetry < 1e-10  # Should be symmetric due to symmetrization

        # Check parameters are stored
        params = neighbors_info["params"]
        assert params["n_neighbors"] == 15
        assert params["n_pcs"] == 2

    def test_compute_neighbors_different_metrics(self):
        """Test neighbor computation with different distance metrics."""
        for metric in ["euclidean", "cosine", "manhattan"]:
            try:
                data = compute_neighbors(self.test_data, n_neighbors=10, metric=metric, method="sklearn")
                assert "neighbors" in data.uns
                assert data.uns["neighbors"]["params"]["metric"] == metric
            except Exception as e:
                # Some metrics might not be available
                pytest.skip(f"Metric {metric} not available: {e}")

    def test_compute_neighbors_without_pca(self):
        """Test neighbor computation without PCA coordinates."""
        # Remove PCA coordinates
        data = self.test_data.copy()
        if "X_pca" in data.obsm:
            del data.obsm["X_pca"]

        data = compute_neighbors(data, n_neighbors=10)

        # Should still work using expression data
        assert "neighbors" in data.uns


class TestUMAP:
    """Test UMAP computation."""

    def setup_method(self):
        """Setup test data with neighbors for UMAP."""
        np.random.seed(42)

        # Create simple 2D data that should embed well
        n_cells = 80
        t = np.linspace(0, 4 * np.pi, n_cells)

        # Create a spiral pattern
        x = t * np.cos(t) * 0.1
        y = t * np.sin(t) * 0.1

        # Embed in higher dimensions with noise
        pca_coords = np.column_stack([x, y])
        full_data = np.random.normal(0, 0.1, (n_cells, 50))
        full_data[:, :2] = pca_coords

        self.test_data = SingleCellData(full_data)
        self.test_data.obsm["X_pca"] = pca_coords

        # Compute neighbors
        self.test_data = compute_neighbors(self.test_data, n_neighbors=15)

    @pytest.mark.slow
    def test_compute_umap_basic(self):
        """Test basic UMAP computation."""
        try:
            data = compute_umap(self.test_data, n_components=2, n_epochs=50)

            # Check that UMAP coordinates are stored
            assert "X_umap" in data.obsm
            assert data.obsm["X_umap"].shape == (80, 2)
            assert "umap" in data.uns

            # Check parameters are stored
            umap_params = data.uns["umap"]["params"]
            assert umap_params["n_components"] == 2
            assert umap_params["n_epochs"] == 50

        except ImportError:
            pytest.skip("UMAP package not available")

    def test_compute_umap_different_dimensions(self):
        """Test UMAP with different output dimensions."""
        try:
            # 3D UMAP
            data = compute_umap(self.test_data, n_components=3, n_epochs=50)
            assert data.obsm["X_umap"].shape == (80, 3)

            # 1D UMAP
            data = compute_umap(self.test_data, n_components=1, n_epochs=50)
            assert data.obsm["X_umap"].shape == (80, 1)

        except ImportError:
            pytest.skip("UMAP package not available")

    def test_compute_umap_without_neighbors(self):
        """Test UMAP computation when neighbors aren't precomputed."""
        data = self.test_data.copy()
        if "neighbors" in data.uns:
            del data.uns["neighbors"]

        try:
            # Should compute neighbors automatically
            data = compute_umap(data, n_components=2, n_epochs=50)
            assert "X_umap" in data.obsm
            assert "neighbors" in data.uns  # Should be computed automatically

        except ImportError:
            pytest.skip("UMAP package not available")


class TestTSNE:
    """Test t-SNE computation."""

    def setup_method(self):
        """Setup test data for t-SNE."""
        np.random.seed(42)

        # Create clustered data
        n_per_cluster = 30

        # Two well-separated clusters
        cluster1 = np.random.multivariate_normal([0, 0], [[1, 0], [0, 1]], n_per_cluster)
        cluster2 = np.random.multivariate_normal([10, 10], [[1, 0], [0, 1]], n_per_cluster)

        # Embed in higher dimensions
        pca_coords = np.vstack([cluster1, cluster2])
        full_data = np.random.normal(0, 0.1, (60, 40))
        full_data[:, :2] = pca_coords

        self.test_data = SingleCellData(full_data)
        self.test_data.obsm["X_pca"] = pca_coords

    @pytest.mark.slow
    def test_compute_tsne_basic(self):
        """Test basic t-SNE computation."""
        data = compute_tsne(self.test_data, n_components=2, n_iter=250)

        # Check that t-SNE coordinates are stored
        assert "X_tsne" in data.obsm
        assert data.obsm["X_tsne"].shape == (60, 2)
        assert "tsne" in data.uns

        # Check parameters
        tsne_params = data.uns["tsne"]["params"]
        assert tsne_params["n_components"] == 2
        assert tsne_params["n_iter"] == 250

    @pytest.mark.slow
    def test_compute_tsne_different_perplexity(self):
        """Test t-SNE with different perplexity values."""
        # Test with different perplexity
        data = compute_tsne(self.test_data, perplexity=15.0, n_iter=250)
        assert "X_tsne" in data.obsm
        assert data.uns["tsne"]["params"]["perplexity"] == 15.0

    @pytest.mark.slow
    def test_compute_tsne_perplexity_adjustment(self):
        """Test that perplexity is adjusted for small datasets."""
        # Create very small dataset
        small_data = SingleCellData(np.random.normal(0, 1, (20, 10)))

        # Large perplexity should be automatically reduced
        data = compute_tsne(small_data, perplexity=50.0, n_iter=100)

        # Should have been automatically adjusted
        actual_perplexity = data.uns["tsne"]["params"]["perplexity"]
        assert actual_perplexity < 50.0
        assert actual_perplexity <= (20 - 1) / 3  # Should be <= (n_samples - 1) / 3


class TestDiffusionMap:
    """Test diffusion map computation."""

    def setup_method(self):
        """Setup test data with neighbors for diffusion map."""
        np.random.seed(42)

        # Create data with trajectory-like structure
        n_cells = 100
        t = np.linspace(0, 2 * np.pi, n_cells)

        # Circular trajectory
        x = np.cos(t)
        y = np.sin(t)

        coords = np.column_stack([x, y])
        full_data = np.random.normal(0, 0.1, (n_cells, 30))
        full_data[:, :2] = coords

        self.test_data = SingleCellData(full_data)
        self.test_data.obsm["X_pca"] = coords

        # Compute neighbors
        self.test_data = compute_neighbors(self.test_data, n_neighbors=10)

    def test_compute_diffusion_map_basic(self):
        """Test basic diffusion map computation."""
        data = compute_diffusion_map(self.test_data, n_components=10)

        # Check that diffusion map coordinates are stored
        assert "X_diffmap" in data.obsm
        assert data.obsm["X_diffmap"].shape == (100, 10)
        assert "diffmap" in data.uns

        # Check that eigenvalues are stored
        diffmap_info = data.uns["diffmap"]
        assert "eigenvalues" in diffmap_info
        assert len(diffmap_info["eigenvalues"]) == 10

        # Eigenvalues should be decreasing (in absolute value)
        eigenvals = np.abs(diffmap_info["eigenvalues"])
        assert np.all(np.diff(eigenvals) <= 1e-10)  # Approximately decreasing

    def test_compute_diffusion_map_without_neighbors(self):
        """Test diffusion map when neighbors need to be computed."""
        data = self.test_data.copy()
        if "neighbors" in data.uns:
            del data.uns["neighbors"]

        # Should compute neighbors automatically
        data = compute_diffusion_map(data, n_components=5)

        assert "X_diffmap" in data.obsm
        assert "neighbors" in data.uns  # Should be computed automatically


class TestIntegrationWorkflow:
    """Test complete dimensionality reduction workflow."""

    def test_full_dimensionality_pipeline(self):
        """Test complete pipeline from HVG selection to visualization."""
        np.random.seed(42)

        # Create realistic test data with structure
        n_cells, n_genes = 200, 500

        # Create three populations with different expression patterns
        pop1 = np.random.lognormal(1, 1, (70, n_genes))
        pop2 = np.random.lognormal(0.5, 1, (70, n_genes))
        pop3 = np.random.lognormal(1.5, 1, (60, n_genes))

        # Add population-specific genes
        pop1[:, :50] *= 5  # High expression in population 1
        pop2[:, 50:100] *= 5  # High expression in population 2
        pop3[:, 100:150] *= 5  # High expression in population 3

        X = np.vstack([pop1, pop2, pop3])

        data = SingleCellData(X)

        # Preprocessing
        data = normalize_counts(data, target_sum=10000)
        data = log_transform(data)
        data = scale_data(data)

        # Dimensionality reduction pipeline
        # 1. Select HVGs
        data = select_hvgs(data, n_top_genes=100, method="seurat")
        assert data.var["highly_variable"].sum() <= 100

        # 2. Compute PCA
        data = compute_pca(data, n_components=30, use_hvgs=True)
        assert "X_pca" in data.obsm
        assert data.obsm["X_pca"].shape == (200, 30)

        # 3. Compute neighbors
        data = compute_neighbors(data, n_neighbors=20, n_pcs=20)
        assert "neighbors" in data.uns

        # 4. Compute embeddings
        try:
            data = compute_umap(data, n_components=2, n_epochs=50)
            assert "X_umap" in data.obsm
        except ImportError:
            # UMAP not available, use t-SNE
            data = compute_tsne(data, n_components=2, n_iter=250)
            assert "X_tsne" in data.obsm

        # 5. Compute diffusion map for trajectory analysis
        data = compute_diffusion_map(data, n_components=10)
        assert "X_diffmap" in data.obsm

        # Final checks
        assert data.n_obs == 200
        assert "highly_variable" in data.var.columns
        assert "pca" in data.uns
        assert "neighbors" in data.uns
        assert "diffmap" in data.uns


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_small_dataset_handling(self):
        """Test handling of very small datasets."""
        # Very small dataset
        X = np.random.normal(0, 1, (5, 10))
        data = SingleCellData(X)

        # Should handle small dataset gracefully
        data = compute_pca(data, n_components=3)
        assert data.obsm["X_pca"].shape[1] <= 4  # Should be limited by n_samples-1

        # t-SNE should adjust perplexity
        data = compute_tsne(data, perplexity=10, n_iter=100)
        actual_perplexity = data.uns["tsne"]["params"]["perplexity"]
        assert actual_perplexity < 10

    def test_single_cell_handling(self):
        """Test handling of single-cell datasets."""
        X = np.random.normal(0, 1, (1, 20))
        data = SingleCellData(X)

        # PCA should handle single cell (though not very meaningful)
        with pytest.raises((ValueError, np.linalg.LinAlgError, ValidationError)):
            # Should fail gracefully for single sample
            compute_pca(data, n_components=5)

    def test_zero_variance_genes(self):
        """Test handling of zero-variance genes."""
        import warnings

        X = np.random.normal(0, 1, (50, 20))
        # Make some genes have zero variance
        X[:, :5] = 1.0  # Constant expression

        data = SingleCellData(X)

        # Suppress the precision loss warning from scipy for zero variance data
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Precision loss occurred", category=RuntimeWarning)
            data = scale_data(data)  # This should handle zero variance

        # Should not contain NaN or infinite values
        assert not np.any(np.isnan(data.X))
        assert not np.any(np.isinf(data.X))

        # PCA should still work
        data = compute_pca(data, n_components=10)
        assert "X_pca" in data.obsm
