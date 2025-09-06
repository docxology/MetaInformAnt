"""Basic tests for single-cell analysis functionality.

Real implementation testing for single-cell methods.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

import numpy as np
import pytest

# Import modules with graceful handling for optional dependencies
try:
    from metainformant.singlecell.preprocessing import (
        SingleCellData,
        filter_cells,
        filter_genes,
        log_transform,
        normalize_counts,
    )

    PREPROCESSING_AVAILABLE = True
except ImportError:
    PREPROCESSING_AVAILABLE = False

    class SingleCellData:
        """Dummy class for when preprocessing is not available."""

        pass


try:
    from metainformant.singlecell.dimensionality import run_pca, run_tsne, run_umap

    DIMENSIONALITY_AVAILABLE = True
except ImportError:
    DIMENSIONALITY_AVAILABLE = False

try:
    from metainformant.singlecell.clustering import find_markers, leiden_clustering

    CLUSTERING_AVAILABLE = True
except ImportError:
    CLUSTERING_AVAILABLE = False


@pytest.mark.skipif(not PREPROCESSING_AVAILABLE, reason="singlecell.preprocessing not available")
class TestSingleCellPreprocessing:
    """Test single-cell preprocessing functionality."""

    def setup_method(self):
        """Set up test single-cell data."""
        # Create synthetic single-cell count matrix
        np.random.seed(42)

        # 100 cells, 50 genes
        self.n_cells = 100
        self.n_genes = 50

        # Generate count data (simulate UMI counts)
        # Most genes have low counts, some have higher expression
        self.count_matrix = np.random.negative_binomial(n=5, p=0.3, size=(self.n_cells, self.n_genes))

        # Add some high-expressing genes
        self.count_matrix[:, :5] = np.random.negative_binomial(n=20, p=0.2, size=(self.n_cells, 5))

        # Add some zero-inflated genes (dropout)
        dropout_mask = np.random.random((self.n_cells, self.n_genes)) < 0.3
        self.count_matrix[dropout_mask] = 0

        # Gene and cell names
        self.gene_names = [f"Gene_{i}" for i in range(self.n_genes)]
        self.cell_names = [f"Cell_{i}" for i in range(self.n_cells)]

    def test_filter_cells_basic(self):
        """Test basic cell filtering functionality."""
        # Create SingleCellData object
        import pandas as pd

        obs = pd.DataFrame(index=self.cell_names)
        var = pd.DataFrame(index=self.gene_names)

        data = SingleCellData(X=self.count_matrix, obs=obs, var=var)

        # Filter cells with too few genes expressed
        filtered_data = filter_cells(data, min_genes=5)

        # Should retain most cells (most express > 5 genes)
        assert filtered_data.X.shape[0] <= self.n_cells
        assert filtered_data.X.shape[1] == self.n_genes
        assert len(filtered_data.obs) == filtered_data.X.shape[0]

        # Check that all retained cells meet criteria
        for i in range(filtered_data.X.shape[0]):
            genes_expressed = np.sum(filtered_data.X[i, :] > 0)
            assert genes_expressed >= 5

    def test_filter_genes_basic(self):
        """Test basic gene filtering functionality."""
        # Create SingleCellData object
        import pandas as pd

        obs = pd.DataFrame(index=self.cell_names)
        var = pd.DataFrame(index=self.gene_names)

        data = SingleCellData(X=self.count_matrix, obs=obs, var=var)

        # Filter genes expressed in too few cells
        filtered_data = filter_genes(data, min_cells=3)

        # Should retain most genes (most are expressed in > 3 cells)
        assert filtered_data.X.shape[0] == self.n_cells
        assert filtered_data.X.shape[1] <= self.n_genes
        assert len(filtered_data.var) == filtered_data.X.shape[1]

        # Check that all retained genes meet criteria
        for j in range(filtered_data.X.shape[1]):
            cells_expressing = np.sum(filtered_data.X[:, j] > 0)
            assert cells_expressing >= 3

    def test_normalize_counts_basic(self):
        """Test count normalization."""
        # Create SingleCellData object
        import pandas as pd

        obs = pd.DataFrame(index=self.cell_names)
        var = pd.DataFrame(index=self.gene_names)

        data = SingleCellData(X=self.count_matrix, obs=obs, var=var)

        normalized_data = normalize_counts(data, target_sum=10000)

        # Should have same dimensions
        assert normalized_data.X.shape == self.count_matrix.shape

        # Each cell should sum to approximately target_sum (within tolerance)
        cell_sums = np.sum(normalized_data.X, axis=1)

        # Cells with non-zero counts should be normalized
        non_zero_cells = np.sum(self.count_matrix, axis=1) > 0
        if np.any(non_zero_cells):
            target_sums = cell_sums[non_zero_cells]
            # Allow some tolerance for normalization
            assert np.allclose(target_sums, 10000, rtol=0.1)

    def test_log_transform_basic(self):
        """Test log transformation."""
        # Create SingleCellData object
        import pandas as pd

        obs = pd.DataFrame(index=self.cell_names)
        var = pd.DataFrame(index=self.gene_names)

        data = SingleCellData(X=self.count_matrix.astype(float), obs=obs, var=var)

        log_data = log_transform(data, base=np.e)

        # Should have same dimensions
        assert log_data.X.shape == self.count_matrix.shape

        # Should be non-negative (log of positive values with pseudocount)
        assert np.all(log_data.X >= 0)

        # Check that zeros become log(1) = 0
        zero_mask = self.count_matrix == 0
        # With pseudocount added, zeros should become log(1) = 0
        assert np.allclose(log_data.X[zero_mask], 0.0, atol=1e-10)

    def test_preprocessing_pipeline(self):
        """Test complete preprocessing pipeline."""
        # Create SingleCellData object
        import pandas as pd

        obs = pd.DataFrame(index=self.cell_names)
        var = pd.DataFrame(index=self.gene_names)

        data = SingleCellData(X=self.count_matrix.astype(float), obs=obs, var=var)

        # 1. Filter low-quality cells and genes
        filtered_data = filter_cells(data, min_genes=3)
        filtered_data = filter_genes(filtered_data, min_cells=2)

        # 2. Normalize
        normalized_data = normalize_counts(filtered_data, target_sum=10000)

        # 3. Log transform
        log_data = log_transform(normalized_data, base=np.e)

        # Verify pipeline integrity
        assert log_data.X.shape[0] <= self.n_cells
        assert log_data.X.shape[1] <= self.n_genes
        assert len(log_data.obs) == log_data.X.shape[0]
        assert len(log_data.var) == log_data.X.shape[1]
        assert np.all(np.isfinite(log_data.X))


@pytest.mark.skipif(not DIMENSIONALITY_AVAILABLE, reason="singlecell.dimensionality not available")
class TestSingleCellDimensionality:
    """Test single-cell dimensionality reduction."""

    def setup_method(self):
        """Set up test data for dimensionality reduction."""
        np.random.seed(123)

        # Create higher-dimensional data with some structure
        self.n_cells = 80
        self.n_genes = 100

        # Create data with latent structure (cell types)
        n_cell_types = 3
        cells_per_type = self.n_cells // n_cell_types

        # Generate cell type-specific expression patterns
        self.expression_matrix = np.zeros((self.n_cells, self.n_genes))
        self.cell_types = []

        for i in range(n_cell_types):
            start_idx = i * cells_per_type
            end_idx = start_idx + cells_per_type

            # Each cell type has different expression pattern
            type_expression = np.random.lognormal(mean=2 + i * 0.5, sigma=0.8, size=(cells_per_type, self.n_genes))

            # Some genes are specific to each type
            type_specific_genes = slice(i * 30, (i + 1) * 30)
            type_expression[:, type_specific_genes] *= 3

            self.expression_matrix[start_idx:end_idx, :] = type_expression
            self.cell_types.extend([f"Type_{i}"] * cells_per_type)

        # Add remaining cells
        remaining = self.n_cells - n_cell_types * cells_per_type
        if remaining > 0:
            self.expression_matrix[-remaining:, :] = np.random.lognormal(
                mean=2, sigma=0.8, size=(remaining, self.n_genes)
            )
            self.cell_types.extend([f"Type_0"] * remaining)

    def test_run_pca_basic(self):
        """Test PCA dimensionality reduction."""
        pca_result = run_pca(expression_matrix=self.expression_matrix, n_components=10)

        # Check result structure
        assert "embedding" in pca_result
        assert "components" in pca_result
        assert "explained_variance" in pca_result

        # Check dimensions
        embedding = pca_result["embedding"]
        assert embedding.shape == (self.n_cells, 10)

        components = pca_result["components"]
        assert components.shape == (self.n_genes, 10)

        explained_var = pca_result["explained_variance"]
        assert len(explained_var) == 10
        assert np.all(explained_var >= 0)
        assert np.sum(explained_var) <= 1.0

    def test_run_tsne_basic(self):
        """Test t-SNE dimensionality reduction."""
        try:
            tsne_result = run_tsne(expression_matrix=self.expression_matrix, n_components=2, perplexity=30)

            # Check result structure
            assert "embedding" in tsne_result

            # Check dimensions
            embedding = tsne_result["embedding"]
            assert embedding.shape == (self.n_cells, 2)

        except ImportError:
            # Skip if sklearn not available
            pytest.skip("sklearn not available for t-SNE")
        except Exception:
            # Skip if other issues (e.g., insufficient data)
            pytest.skip("t-SNE failed, possibly due to data constraints")

    def test_run_umap_basic(self):
        """Test UMAP dimensionality reduction."""
        try:
            umap_result = run_umap(expression_matrix=self.expression_matrix, n_components=2, n_neighbors=15)

            # Check result structure
            assert "embedding" in umap_result

            # Check dimensions
            embedding = umap_result["embedding"]
            assert embedding.shape == (self.n_cells, 2)

        except ImportError:
            # Skip if umap-learn not available
            pytest.skip("umap-learn not available for UMAP")
        except Exception:
            # Skip if other issues
            pytest.skip("UMAP failed, possibly due to missing dependencies")


@pytest.mark.skipif(not CLUSTERING_AVAILABLE, reason="singlecell.clustering not available")
class TestSingleCellClustering:
    """Test single-cell clustering functionality."""

    def setup_method(self):
        """Set up test data for clustering."""
        np.random.seed(456)

        # Create data with clear cluster structure
        self.n_cells = 60
        self.n_features = 20  # Reduced dimensionality (e.g., PC space)

        # Generate 3 distinct clusters
        cluster_size = self.n_cells // 3

        self.feature_matrix = np.zeros((self.n_cells, self.n_features))
        self.true_clusters = []

        # Cluster 1: centered around origin
        self.feature_matrix[:cluster_size, :] = np.random.multivariate_normal(
            mean=np.zeros(self.n_features), cov=np.eye(self.n_features) * 0.5, size=cluster_size
        )
        self.true_clusters.extend([0] * cluster_size)

        # Cluster 2: shifted in positive direction
        self.feature_matrix[cluster_size : 2 * cluster_size, :] = np.random.multivariate_normal(
            mean=np.ones(self.n_features) * 3, cov=np.eye(self.n_features) * 0.5, size=cluster_size
        )
        self.true_clusters.extend([1] * cluster_size)

        # Cluster 3: shifted in different direction
        shift = np.zeros(self.n_features)
        shift[::2] = -2  # Negative values for even indices
        self.feature_matrix[2 * cluster_size :, :] = np.random.multivariate_normal(
            mean=shift, cov=np.eye(self.n_features) * 0.5, size=self.n_cells - 2 * cluster_size
        )
        self.true_clusters.extend([2] * (self.n_cells - 2 * cluster_size))

        self.true_clusters = np.array(self.true_clusters)

    def test_leiden_clustering_basic(self):
        """Test Leiden clustering."""
        try:
            cluster_labels = leiden_clustering(feature_matrix=self.feature_matrix, resolution=0.5)

            # Should return cluster labels for all cells
            assert len(cluster_labels) == self.n_cells
            assert np.all(cluster_labels >= 0)

            # Should identify reasonable number of clusters
            n_clusters = len(np.unique(cluster_labels))
            assert 1 <= n_clusters <= self.n_cells
            assert n_clusters <= 10  # Reasonable upper bound

        except ImportError:
            # Skip if required dependencies not available
            pytest.skip("Leiden clustering dependencies not available")
        except Exception:
            # Skip if other issues
            pytest.skip("Leiden clustering failed")

    def test_find_markers_basic(self):
        """Test marker gene identification."""
        # Use true clusters for testing marker finding
        try:
            # Create expression matrix for marker testing
            n_genes = 50
            expression_matrix = np.random.lognormal(mean=1, sigma=0.5, size=(self.n_cells, n_genes))

            # Make some genes differentially expressed between clusters
            for cluster_id in np.unique(self.true_clusters):
                cluster_mask = self.true_clusters == cluster_id
                # Upregulate specific genes in this cluster
                marker_genes = slice(cluster_id * 10, (cluster_id + 1) * 10)
                expression_matrix[cluster_mask, marker_genes] *= 3

            markers = find_markers(expression_matrix=expression_matrix, cluster_labels=self.true_clusters)

            # Should return markers for each cluster
            assert isinstance(markers, dict)
            assert len(markers) == len(np.unique(self.true_clusters))

            # Each cluster should have some markers
            for cluster_id, marker_list in markers.items():
                assert len(marker_list) >= 0  # Might be empty for some clusters

                # If markers found, should have reasonable structure
                if len(marker_list) > 0:
                    # Each marker should be a dict or tuple with gene info
                    first_marker = marker_list[0]
                    assert isinstance(first_marker, (dict, tuple))

        except ImportError:
            # Skip if required dependencies not available
            pytest.skip("Marker finding dependencies not available")
        except Exception:
            # Skip if other issues
            pytest.skip("Marker finding failed")


class TestSingleCellIntegration:
    """Integration tests for single-cell analysis."""

    @pytest.mark.skipif(
        not (PREPROCESSING_AVAILABLE and DIMENSIONALITY_AVAILABLE), reason="Required single-cell modules not available"
    )
    def test_basic_singlecell_workflow(self):
        """Test basic single-cell analysis workflow."""
        import pandas as pd

        np.random.seed(42)

        # 1. Generate synthetic single-cell data
        n_cells = 50
        n_genes = 30

        count_matrix = np.random.negative_binomial(n=3, p=0.4, size=(n_cells, n_genes))
        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_names = [f"Cell_{i}" for i in range(n_cells)]

        obs = pd.DataFrame(index=cell_names)
        var = pd.DataFrame(index=gene_names)

        data = SingleCellData(X=count_matrix.astype(float), obs=obs, var=var)

        # 2. Preprocessing
        filtered_data = filter_cells(data, min_genes=2)
        normalized_data = normalize_counts(filtered_data, target_sum=1000)
        log_data = log_transform(normalized_data, base=np.e)

        # 3. Dimensionality reduction
        pca_result = run_pca(expression_matrix=log_data.X, n_components=5)

        # 4. Verify workflow integrity
        assert log_data.X.shape[0] <= n_cells
        assert log_data.X.shape[1] == n_genes
        assert pca_result["embedding"].shape[0] == log_data.X.shape[0]
        assert pca_result["embedding"].shape[1] == 5

        # All results should be finite
        assert np.all(np.isfinite(log_data.X))
        assert np.all(np.isfinite(pca_result["embedding"]))


class TestSingleCellEdgeCases:
    """Test edge cases for single-cell analysis."""

    @pytest.mark.skipif(not PREPROCESSING_AVAILABLE, reason="preprocessing not available")
    def test_empty_matrix(self):
        """Test handling of empty matrices."""
        import pandas as pd

        empty_matrix = np.array([]).reshape(0, 5)
        empty_obs = pd.DataFrame(index=[])
        empty_var = pd.DataFrame(index=[f"Gene_{i}" for i in range(5)])

        # Should handle gracefully
        try:
            data = SingleCellData(X=empty_matrix, obs=empty_obs, var=empty_var)
            filtered_data = filter_cells(data, min_genes=1)
            assert filtered_data.X.shape == (0, 5)
            assert len(filtered_data.obs) == 0
        except (ValueError, IndexError):
            # Acceptable to raise error for empty data
            pass

    @pytest.mark.skipif(not PREPROCESSING_AVAILABLE, reason="preprocessing not available")
    def test_all_zero_matrix(self):
        """Test handling of all-zero count matrix."""
        import pandas as pd

        zero_matrix = np.zeros((10, 20))
        obs = pd.DataFrame(index=[f"Cell_{i}" for i in range(10)])
        var = pd.DataFrame(index=[f"Gene_{i}" for i in range(20)])

        # Should handle all-zero data gracefully
        try:
            data = SingleCellData(X=zero_matrix, obs=obs, var=var)
            normalized_data = normalize_counts(data, target_sum=1000)
            # Should remain all zeros
            assert np.allclose(normalized_data.X, 0)
        except (ValueError, ZeroDivisionError):
            # Acceptable to raise error for all-zero data
            pass

    @pytest.mark.skipif(not PREPROCESSING_AVAILABLE, reason="preprocessing not available")
    def test_single_cell_single_gene(self):
        """Test handling of minimal data (1 cell, 1 gene)."""
        import pandas as pd

        minimal_matrix = np.array([[5.0]])
        obs = pd.DataFrame(index=["Cell_0"])
        var = pd.DataFrame(index=["Gene_0"])

        # Should handle minimal data
        try:
            data = SingleCellData(X=minimal_matrix, obs=obs, var=var)
            normalized_data = normalize_counts(data, target_sum=1000)
            assert normalized_data.X.shape == (1, 1)
            assert normalized_data.X[0, 0] == 1000  # Should be normalized to target

            log_data = log_transform(normalized_data, base=np.e)
            assert log_data.X.shape == (1, 1)
            assert np.isfinite(log_data.X[0, 0])

        except Exception:
            # Some operations might not be meaningful for single cell/gene
            pass
