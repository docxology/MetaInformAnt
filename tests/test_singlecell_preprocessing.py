"""Tests for single-cell preprocessing module.

Real implementation testing for single-cell data loading, quality control,
and preprocessing functions. No mocking used - all tests use real data
and computational methods.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

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

from metainformant.singlecell.preprocessing import (
    SingleCellData,
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    load_count_matrix,
    log_transform,
    normalize_counts,
    scale_data,
)


class TestSingleCellData:
    """Test the SingleCellData container class."""

    def test_singlecelldata_initialization(self):
        """Test basic initialization of SingleCellData."""
        # Create test data
        X = np.random.poisson(2, (100, 50)).astype(float)

        data = SingleCellData(X)

        assert data.n_obs == 100
        assert data.n_vars == 50
        assert data.X.shape == (100, 50)
        assert len(data.obs) == 100
        assert len(data.var) == 50

    def test_singlecelldata_with_metadata(self):
        """Test initialization with metadata."""
        X = np.random.poisson(2, (50, 20)).astype(float)
        obs = pd.DataFrame(
            {"sample": [f"sample_{i//10}" for i in range(50)], "batch": ["A"] * 25 + ["B"] * 25},
            index=[f"cell_{i}" for i in range(50)],
        )
        var = pd.DataFrame({"gene_name": [f"Gene_{i}" for i in range(20)]}, index=[f"ENSG_{i:05d}" for i in range(20)])

        data = SingleCellData(X, obs=obs, var=var)

        assert data.n_obs == 50
        assert data.n_vars == 20
        assert "sample" in data.obs.columns
        assert "gene_name" in data.var.columns
        assert data.obs.index[0] == "cell_0"
        assert data.var.index[0] == "ENSG_00000"

    def test_singlecelldata_sparse_matrix(self):
        """Test with sparse matrix input."""
        X_dense = np.random.poisson(1, (30, 15)).astype(float)
        X_sparse = sparse.csr_matrix(X_dense)

        data = SingleCellData(X_sparse)

        assert data.n_obs == 30
        assert data.n_vars == 15
        assert sparse.issparse(data.X)

    def test_singlecelldata_copy(self):
        """Test copying SingleCellData objects."""
        X = np.random.poisson(2, (20, 10)).astype(float)
        data = SingleCellData(X)
        data.obs["test"] = range(20)

        data_copy = data.copy()

        # Modify original
        data.obs.loc[0, "test"] = 999

        # Copy should be unchanged
        assert data_copy.obs.loc[0, "test"] == 0
        assert data.obs.loc[0, "test"] == 999


class TestDataLoading:
    """Test data loading functions."""

    def test_load_count_matrix_csv(self, tmp_path):
        """Test loading CSV count matrix."""
        # Create test CSV file
        data = pd.DataFrame(
            np.random.poisson(3, (20, 10)),
            index=[f"cell_{i}" for i in range(20)],
            columns=[f"gene_{i}" for i in range(10)],
        )

        csv_path = tmp_path / "test_counts.csv"
        data.to_csv(csv_path)

        # Load data
        sc_data = load_count_matrix(csv_path, format="csv")

        assert sc_data.n_obs == 20
        assert sc_data.n_vars == 10
        assert sparse.issparse(sc_data.X)
        assert len(sc_data.obs) == 20
        assert len(sc_data.var) == 10

    def test_load_count_matrix_tsv(self, tmp_path):
        """Test loading TSV count matrix."""
        # Create test TSV file
        data = pd.DataFrame(
            np.random.poisson(2, (15, 8)),
            index=[f"CELL_{i}" for i in range(15)],
            columns=[f"GENE_{i}" for i in range(8)],
        )

        tsv_path = tmp_path / "test_counts.tsv"
        data.to_csv(tsv_path, sep="\t")

        # Load data
        sc_data = load_count_matrix(tsv_path, format="tsv")

        assert sc_data.n_obs == 15
        assert sc_data.n_vars == 8
        assert "CELL_0" in sc_data.obs.index.astype(str)
        assert "GENE_0" in sc_data.var.index.astype(str)

    def test_load_count_matrix_transpose(self, tmp_path):
        """Test transpose functionality."""
        # Create test data (genes x cells format)
        data = pd.DataFrame(
            np.random.poisson(2, (10, 20)),  # 10 genes, 20 cells
            index=[f"gene_{i}" for i in range(10)],
            columns=[f"cell_{i}" for i in range(20)],
        )

        csv_path = tmp_path / "test_counts_transposed.csv"
        data.to_csv(csv_path)

        # Load with transpose=True
        sc_data = load_count_matrix(csv_path, format="csv", transpose=True)

        # Should now be 20 cells x 10 genes
        assert sc_data.n_obs == 20
        assert sc_data.n_vars == 10


class TestQualityControl:
    """Test quality control functions."""

    def setup_method(self):
        """Setup test data with known characteristics."""
        # Create test data with mitochondrial and ribosomal genes
        np.random.seed(42)

        # 100 cells, 50 genes
        n_cells, n_genes = 100, 50

        # Create expression matrix
        X = np.random.poisson(2, (n_cells, n_genes)).astype(float)

        # Create gene names with some MT and ribosomal genes
        gene_names = []
        for i in range(n_genes):
            if i < 5:
                gene_names.append(f"MT-CO{i+1}")  # Mitochondrial genes
            elif i < 10:
                gene_names.append(f"RPS{i-4}")  # Ribosomal genes
            else:
                gene_names.append(f"Gene_{i}")

        var = pd.DataFrame({"gene_name": gene_names})
        self.test_data = SingleCellData(X, var=var)

    def test_calculate_qc_metrics(self):
        """Test QC metrics calculation."""
        data = calculate_qc_metrics(self.test_data)

        # Check that QC metrics are calculated
        assert "total_counts" in data.obs.columns
        assert "n_genes" in data.obs.columns
        assert "pct_mt" in data.obs.columns
        assert "pct_ribo" in data.obs.columns

        # Check gene-level metrics
        assert "n_cells" in data.var.columns
        assert "total_counts" in data.var.columns
        assert "mean_expression" in data.var.columns

        # Check values are reasonable
        assert all(data.obs["total_counts"] >= 0)
        assert all(data.obs["n_genes"] >= 0)
        assert all((data.obs["pct_mt"] >= 0) & (data.obs["pct_mt"] <= 100))
        assert all((data.obs["pct_ribo"] >= 0) & (data.obs["pct_ribo"] <= 100))

    def test_filter_cells(self):
        """Test cell filtering."""
        data = calculate_qc_metrics(self.test_data)

        # Apply filtering
        filtered_data = filter_cells(data, min_genes=5, min_counts=10, max_pct_mt=50.0)

        # Should have fewer or equal cells
        assert filtered_data.n_obs <= data.n_obs
        assert filtered_data.n_vars == data.n_vars

        # Check that filtered cells meet criteria
        assert all(filtered_data.obs["n_genes"] >= 5)
        assert all(filtered_data.obs["total_counts"] >= 10)
        assert all(filtered_data.obs["pct_mt"] <= 50.0)

    def test_filter_genes(self):
        """Test gene filtering."""
        data = calculate_qc_metrics(self.test_data)

        # Apply filtering
        filtered_data = filter_genes(data, min_cells=2, min_counts=5)

        # Should have fewer or equal genes
        assert filtered_data.n_vars <= data.n_vars
        assert filtered_data.n_obs == data.n_obs

        # Check that filtered genes meet criteria
        assert all(filtered_data.var["n_cells"] >= 2)
        assert all(filtered_data.var["total_counts"] >= 5)


class TestNormalization:
    """Test normalization functions."""

    def setup_method(self):
        """Setup test data for normalization."""
        np.random.seed(42)

        # Create test data with varying library sizes
        X = np.random.poisson(2, (50, 30)).astype(float)

        # Make some cells have much higher counts
        X[:10, :] *= 5  # High count cells
        X[10:20, :] *= 0.2  # Low count cells

        self.test_data = SingleCellData(X)

    def test_normalize_counts_total_count(self):
        """Test total count normalization."""
        normalized_data = normalize_counts(self.test_data, target_sum=10000, method="total_count")

        # Check that library sizes are normalized
        if sparse.issparse(normalized_data.X):
            library_sizes = np.array(normalized_data.X.sum(axis=1)).flatten()
        else:
            library_sizes = normalized_data.X.sum(axis=1)

        # Should be close to target sum (allowing for cells with zero counts)
        non_zero_cells = library_sizes > 0
        if non_zero_cells.any():
            assert np.allclose(library_sizes[non_zero_cells], 10000, rtol=1e-10)

        # Check that normalization info is stored
        assert "normalization" in normalized_data.uns
        assert normalized_data.uns["normalization"]["method"] == "total_count"

    def test_normalize_counts_median(self):
        """Test median normalization."""
        normalized_data = normalize_counts(self.test_data, method="median")

        # Check that data is normalized
        assert normalized_data.X.shape == self.test_data.X.shape
        assert "normalization" in normalized_data.uns
        assert normalized_data.uns["normalization"]["method"] == "median"

    def test_log_transform(self):
        """Test log transformation."""
        # First normalize
        normalized_data = normalize_counts(self.test_data)

        # Then log transform
        log_data = log_transform(normalized_data, base=np.e)

        # Check that transformation is applied
        assert log_data.X.shape == normalized_data.X.shape
        assert "log_transformed" in log_data.uns
        assert log_data.uns["log_transformed"] is True
        assert log_data.uns["log_base"] == np.e

        # Values should be non-negative (log of positive values + 1)
        if sparse.issparse(log_data.X):
            assert np.all(log_data.X.data >= 0)
        else:
            assert np.all(log_data.X >= 0)

    def test_scale_data(self):
        """Test data scaling."""
        # Prepare data
        normalized_data = normalize_counts(self.test_data)
        log_data = log_transform(normalized_data)

        # Scale data
        scaled_data = scale_data(log_data, zero_center=True, max_value=10)

        # Check scaling properties
        assert scaled_data.X.shape == log_data.X.shape
        assert "scaled" in scaled_data.uns

        # Values should be clipped to max_value
        assert np.all(np.abs(scaled_data.X) <= 10)

        # Should be approximately zero-centered (per gene)
        gene_means = np.mean(scaled_data.X, axis=0)
        assert np.allclose(gene_means, 0, atol=1e-10)


class TestIntegrationWorkflow:
    """Test complete preprocessing workflow."""

    def test_full_preprocessing_pipeline(self, tmp_path):
        """Test complete preprocessing pipeline with real data."""
        # Create realistic test data
        np.random.seed(42)
        n_cells, n_genes = 200, 100

        # Simulate count data with some structure
        X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes)).astype(float)

        # Add some dropout (zeros)
        dropout_mask = np.random.random((n_cells, n_genes)) < 0.7
        X[dropout_mask] = 0

        # Create gene names including MT and ribosomal genes
        gene_names = []
        for i in range(n_genes):
            if i < 10:
                gene_names.append(f"MT-{i}")
            elif i < 20:
                gene_names.append(f"RPS{i-9}")
            else:
                gene_names.append(f"Gene_{i}")

        var = pd.DataFrame({"gene_name": gene_names})
        obs = pd.DataFrame({"sample": [f"sample_{i//40}" for i in range(n_cells)]})

        # Save as CSV and load
        data_df = pd.DataFrame(X, columns=gene_names)
        csv_path = tmp_path / "test_pipeline.csv"
        data_df.to_csv(csv_path)

        # Load data
        data = load_count_matrix(csv_path, format="csv")
        data.var = var  # Add gene metadata
        data.obs = obs  # Add cell metadata

        # Run full pipeline
        # 1. QC metrics
        data = calculate_qc_metrics(data)
        assert "pct_mt" in data.obs.columns

        # 2. Filter cells and genes
        data = filter_cells(data, min_genes=5, min_counts=10)
        data = filter_genes(data, min_cells=3)

        # 3. Normalize and transform
        data = normalize_counts(data, target_sum=10000)
        data = log_transform(data)
        data = scale_data(data, max_value=10)

        # Check final state
        assert data.n_obs > 0  # Should have some cells remaining
        assert data.n_vars > 0  # Should have some genes remaining
        assert "scaled" in data.uns
        assert "log_transformed" in data.uns
        assert "normalization" in data.uns

        # Data should be properly scaled
        gene_means = np.mean(data.X, axis=0)
        assert np.allclose(gene_means, 0, atol=1e-10)
        assert np.all(np.abs(data.X) <= 10)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_data_handling(self):
        """Test handling of empty or minimal data."""
        # Single cell, single gene
        X = np.array([[5.0]])
        data = SingleCellData(X)

        # Should handle QC calculation
        data = calculate_qc_metrics(data)
        assert data.obs["total_counts"].iloc[0] == 5.0
        assert data.obs["n_genes"].iloc[0] == 1.0

    def test_all_zero_data(self):
        """Test handling of all-zero data."""
        X = np.zeros((10, 5))
        data = SingleCellData(X)

        # Should handle normalization of zero data
        normalized = normalize_counts(data)

        # All values should remain zero
        if sparse.issparse(normalized.X):
            assert normalized.X.nnz == 0
        else:
            assert np.all(normalized.X == 0)

    def test_single_nonzero_gene(self):
        """Test handling of data with only one expressing gene."""
        X = np.zeros((20, 10))
        X[:, 0] = np.random.poisson(3, 20)  # Only first gene expressed

        data = SingleCellData(X)
        data = calculate_qc_metrics(data)

        # Should calculate metrics correctly
        assert all(data.obs["n_genes"] <= 1)  # At most 1 gene per cell
        assert data.var["n_cells"].iloc[0] > 0  # First gene should have expressing cells
        assert all(data.var["n_cells"].iloc[1:] == 0)  # Other genes should have zero cells

    def test_dimension_mismatch_error(self):
        """Test that dimension mismatches raise appropriate errors."""
        X = np.random.poisson(2, (10, 5))
        obs = pd.DataFrame(index=range(15))  # Wrong number of cells

        with pytest.raises(ValueError, match="obs length"):
            SingleCellData(X, obs=obs)

    def test_invalid_normalization_method(self):
        """Test invalid normalization method raises error."""
        X = np.random.poisson(2, (10, 5))
        data = SingleCellData(X)

        with pytest.raises(ValueError, match="Unknown normalization method"):
            normalize_counts(data, method="invalid_method")
