"""Comprehensive tests for multiomics module.

Tests cover data integration, sample alignment, joint analysis methods,
and integration workflows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.multiomics.analysis.integration import (
    MultiOmicsData,
    canonical_correlation,
    integrate_omics_data,
    joint_nmf,
    joint_pca,
)


class TestMultiOmicsData:
    """Tests for MultiOmicsData class."""

    def test_basic_creation(self):
        """Test creating MultiOmicsData with two layers."""
        np.random.seed(42)
        genomics = pd.DataFrame(np.random.randn(10, 50), index=[f"S{i}" for i in range(10)])
        transcriptomics = pd.DataFrame(np.random.randn(10, 100), index=[f"S{i}" for i in range(10)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        assert data.n_samples == 10
        assert len(data.samples) == 10
        assert set(data.layer_names) == {"genomics", "transcriptomics"}

    def test_sample_alignment(self):
        """Test automatic sample alignment across layers."""
        # Create data with overlapping but not identical samples
        genomics = pd.DataFrame(np.random.randn(5, 10), index=["S0", "S1", "S2", "S3", "S4"])
        transcriptomics = pd.DataFrame(np.random.randn(5, 10), index=["S1", "S2", "S3", "S4", "S5"])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        # Only common samples should be retained
        assert data.n_samples == 4
        assert set(data.samples) == {"S1", "S2", "S3", "S4"}

    def test_get_layer(self):
        """Test retrieving specific omics layer."""
        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        genomics_retrieved = data.get_layer("genomics")
        assert genomics_retrieved.shape == (5, 10)
        assert list(genomics_retrieved.index) == data.samples

    def test_get_layer_error(self):
        """Test error when requesting non-existent layer."""
        data = MultiOmicsData(genomics=pd.DataFrame(np.random.randn(5, 10)))

        with pytest.raises(KeyError, match="Layer 'transcriptomics' not available"):
            data.get_layer("transcriptomics")

    def test_subset_samples(self):
        """Test subsetting samples."""
        genomics = pd.DataFrame(np.random.randn(10, 50), index=[f"S{i}" for i in range(10)])
        transcriptomics = pd.DataFrame(np.random.randn(10, 100), index=[f"S{i}" for i in range(10)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        subset = data.subset_samples(["S0", "S1", "S2"])

        assert subset.n_samples == 3
        assert set(subset.samples) == {"S0", "S1", "S2"}
        assert subset.get_layer("genomics").shape == (3, 50)

    def test_subset_features(self):
        """Test subsetting features."""
        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)], columns=[f"SNP{i}" for i in range(10)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)], columns=[f"GENE{i}" for i in range(20)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        subset = data.subset_features({"transcriptomics": ["GENE0", "GENE1", "GENE2"]})

        assert subset.get_layer("transcriptomics").shape == (5, 3)
        assert subset.get_layer("genomics").shape == (5, 10)  # Unchanged

    def test_metadata_integration(self):
        """Test metadata integration."""
        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)])
        metadata = pd.DataFrame({"age": [20, 30, 40, 50, 60], "sex": ["M", "F", "M", "F", "M"]}, index=[f"S{i}" for i in range(5)])

        data = MultiOmicsData(genomics=genomics, metadata=metadata)

        assert data.metadata is not None
        assert len(data.metadata) == 5


class TestIntegrationFunction:
    """Tests for integrate_omics_data function."""

    def test_integrate_from_dataframes(self):
        """Test integration from pandas DataFrames."""
        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)])

        data_dict = {"genomics": genomics, "transcriptomics": transcriptomics}

        integrated = integrate_omics_data(data_dict)

        assert integrated.n_samples == 5
        assert set(integrated.layer_names) == {"genomics", "transcriptomics"}

    def test_integrate_from_files(self, tmp_path: Path):
        """Test integration from file paths."""
        # Create test CSV files
        genomics_file = tmp_path / "genomics.csv"
        transcriptomics_file = tmp_path / "transcriptomics.csv"

        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)])

        genomics.to_csv(genomics_file)
        transcriptomics.to_csv(transcriptomics_file)

        data_dict = {"genomics": str(genomics_file), "transcriptomics": str(transcriptomics_file)}

        integrated = integrate_omics_data(data_dict)

        assert integrated.n_samples == 5
        assert set(integrated.layer_names) == {"genomics", "transcriptomics"}


class TestJointPCA:
    """Tests for joint PCA functionality."""

    def test_joint_pca_basic(self):
        """Test basic joint PCA."""
        np.random.seed(42)
        genomics = pd.DataFrame(np.random.randn(10, 50), index=[f"S{i}" for i in range(10)])
        transcriptomics = pd.DataFrame(np.random.randn(10, 100), index=[f"S{i}" for i in range(10)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        embeddings, loadings, variance = joint_pca(data, n_components=5)

        assert embeddings.shape == (10, 5)
        assert "genomics" in loadings
        assert "transcriptomics" in loadings
        assert len(variance) == 5
        assert all(v >= 0 for v in variance)  # Variance should be non-negative

    def test_joint_pca_with_weights(self):
        """Test joint PCA with layer weights."""
        genomics = pd.DataFrame(np.random.randn(5, 10), index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        embeddings, loadings, variance = joint_pca(data, n_components=3, layer_weights={"genomics": 2.0, "transcriptomics": 1.0})

        assert embeddings.shape == (5, 3)
        assert len(variance) == 3

    def test_joint_pca_standardization(self):
        """Test joint PCA with and without standardization."""
        # Create data with different scales
        genomics = pd.DataFrame(np.random.randn(5, 10) * 10, index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.random.randn(5, 20), index=[f"S{i}" for i in range(5)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        # With standardization
        emb1, _, _ = joint_pca(data, n_components=3, standardize=True)

        # Without standardization
        emb2, _, _ = joint_pca(data, n_components=3, standardize=False)

        # Results should differ
        assert not np.allclose(emb1, emb2)


class TestJointNMF:
    """Tests for joint NMF functionality."""

    def test_joint_nmf_basic(self):
        """Test basic joint NMF."""
        np.random.seed(42)
        # NMF requires non-negative data
        genomics = pd.DataFrame(np.abs(np.random.randn(10, 50)), index=[f"S{i}" for i in range(10)])
        transcriptomics = pd.DataFrame(np.abs(np.random.randn(10, 100)), index=[f"S{i}" for i in range(10)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        W, H = joint_nmf(data, n_components=5, max_iter=50, random_state=42)

        assert W.shape == (10, 5)
        assert "genomics" in H
        assert "transcriptomics" in H
        assert H["genomics"].shape == (5, 50)
        assert H["transcriptomics"].shape == (5, 100)

        # Factors should be non-negative
        assert np.all(W >= 0)
        assert np.all(H["genomics"] >= 0)

    def test_joint_nmf_with_regularization(self):
        """Test joint NMF with different regularization."""
        genomics = pd.DataFrame(np.abs(np.random.randn(5, 10)), index=[f"S{i}" for i in range(5)])
        transcriptomics = pd.DataFrame(np.abs(np.random.randn(5, 20)), index=[f"S{i}" for i in range(5)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        W, H = joint_nmf(data, n_components=3, regularization=0.1, max_iter=30, random_state=42)

        assert W.shape == (5, 3)


class TestCanonicalCorrelation:
    """Tests for canonical correlation analysis."""

    def test_canonical_correlation_basic(self):
        """Test basic canonical correlation analysis."""
        np.random.seed(42)
        # Create correlated data
        genomics = pd.DataFrame(np.random.randn(20, 30), index=[f"S{i}" for i in range(20)])
        # Make transcriptomics partially correlated with genomics
        transcriptomics = pd.DataFrame(genomics.values[:, :20] * 0.7 + np.random.randn(20, 20) * 0.3, index=[f"S{i}" for i in range(20)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        X_c, Y_c, X_w, Y_w, correlations = canonical_correlation(data, ("genomics", "transcriptomics"), n_components=5)

        assert X_c.shape == (20, 5)
        assert Y_c.shape == (20, 5)
        assert X_w.shape == (30, 5)
        assert Y_w.shape == (20, 5)
        assert len(correlations) == 5
        assert all(0 <= c <= 1 for c in correlations)  # Correlations should be in [0, 1]

    def test_canonical_correlation_invalid_layer(self):
        """Test error handling for invalid layer names."""
        genomics = pd.DataFrame(np.random.randn(10, 20), index=[f"S{i}" for i in range(10)])
        data = MultiOmicsData(genomics=genomics)

        with pytest.raises(ValueError, match="Layer transcriptomics not found"):
            canonical_correlation(data, ("genomics", "transcriptomics"))

    def test_canonical_correlation_high_correlation(self):
        """Test CCA with highly correlated layers."""
        np.random.seed(42)
        base_data = np.random.randn(15, 25)

        genomics = pd.DataFrame(base_data, index=[f"S{i}" for i in range(15)])
        # Transcriptomics is very similar to genomics
        transcriptomics = pd.DataFrame(base_data[:, :20] * 0.95 + np.random.randn(15, 20) * 0.05, index=[f"S{i}" for i in range(15)])

        data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)

        _, _, _, _, correlations = canonical_correlation(data, ("genomics", "transcriptomics"), n_components=3)

        # First canonical correlation should be high
        assert correlations[0] > 0.8
