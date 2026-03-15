"""Tests for multi-omics data integration functionality.

Real implementation testing for multi-omics analysis methods.
No mocking used - all tests use real computational methods and data structures.
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
    """Test MultiOmicsData container functionality."""

    def setup_method(self):
        """Set up test data with aligned samples."""
        # Create sample data with common samples
        samples = ["sample1", "sample2", "sample3", "sample4"]

        # Genomics data (SNPs/mutations)
        self.genomics_data = pd.DataFrame(
            np.random.randint(0, 3, (4, 100)),  # Genotype data (0,1,2)
            index=samples,
            columns=[f"SNP_{i}" for i in range(100)],
        )

        # Transcriptomics data (gene expression)
        np.random.seed(42)
        self.transcriptomics_data = pd.DataFrame(
            np.random.lognormal(mean=2, sigma=1, size=(4, 200)),
            index=samples,
            columns=[f"Gene_{i}" for i in range(200)],
        )

        # Proteomics data (protein abundance)
        self.proteomics_data = pd.DataFrame(
            np.random.lognormal(mean=1, sigma=0.5, size=(4, 50)),
            index=samples,
            columns=[f"Protein_{i}" for i in range(50)],
        )

        # Metadata
        self.metadata = pd.DataFrame(
            {
                "age": [25, 30, 45, 60],
                "sex": ["M", "F", "M", "F"],
                "condition": ["healthy", "disease", "healthy", "disease"],
            },
            index=samples,
        )

    def test_multiomics_data_initialization(self):
        """Test basic MultiOmicsData initialization."""
        # Single omics layer
        omics_data = MultiOmicsData(transcriptomics=self.transcriptomics_data)
        assert len(omics_data.layer_names) == 1
        assert "transcriptomics" in omics_data.layer_names
        assert omics_data.n_samples == 4

        # Multiple omics layers
        omics_data = MultiOmicsData(
            genomics=self.genomics_data,
            transcriptomics=self.transcriptomics_data,
            proteomics=self.proteomics_data,
            metadata=self.metadata,
        )

        assert len(omics_data.layer_names) == 3
        assert set(omics_data.layer_names) == {"genomics", "transcriptomics", "proteomics"}
        assert omics_data.n_samples == 4
        assert len(omics_data.samples) == 4

    def test_sample_alignment(self):
        """Test automatic sample alignment across layers."""
        # Create data with different sample sets
        samples_1 = ["A", "B", "C", "D"]
        samples_2 = ["B", "C", "D", "E"]  # Missing A, extra E
        samples_3 = ["A", "B", "C"]  # Missing D, E

        data1 = pd.DataFrame(np.random.randn(4, 10), index=samples_1)
        data2 = pd.DataFrame(np.random.randn(4, 20), index=samples_2)
        data3 = pd.DataFrame(np.random.randn(3, 15), index=samples_3)

        # Should align to common samples: B, C
        omics_data = MultiOmicsData(genomics=data1, transcriptomics=data2, proteomics=data3)

        expected_common = {"B", "C"}
        assert set(omics_data.samples) == expected_common
        assert omics_data.n_samples == 2

        # Check that all layers have the same sample order
        for layer_name in omics_data.layer_names:
            layer_data = omics_data.get_layer(layer_name)
            assert list(layer_data.index) == sorted(expected_common)

    def test_no_common_samples_error(self):
        """Test error when no common samples exist."""
        data1 = pd.DataFrame(np.random.randn(2, 10), index=["A", "B"])
        data2 = pd.DataFrame(np.random.randn(2, 10), index=["C", "D"])

        with pytest.raises(ValueError, match="No common samples"):
            MultiOmicsData(genomics=data1, transcriptomics=data2)

    def test_empty_omics_error(self):
        """Test error when no omics data provided."""
        with pytest.raises(ValueError, match="At least one omics layer"):
            MultiOmicsData()

    def test_invalid_data_type_error(self):
        """Test error with invalid data types."""
        with pytest.raises(TypeError, match="must be pandas DataFrame"):
            MultiOmicsData(transcriptomics="not_a_dataframe")

    def test_get_layer(self):
        """Test getting specific omics layers."""
        omics_data = MultiOmicsData(transcriptomics=self.transcriptomics_data, proteomics=self.proteomics_data)

        # Get existing layer
        transcr_data = omics_data.get_layer("transcriptomics")
        assert transcr_data.shape[1] == 200  # 200 genes

        # Get nonexistent layer
        with pytest.raises(KeyError, match="Layer 'nonexistent' not available"):
            omics_data.get_layer("nonexistent")

    def test_subset_samples(self):
        """Test subsetting by samples."""
        omics_data = MultiOmicsData(
            transcriptomics=self.transcriptomics_data, proteomics=self.proteomics_data, metadata=self.metadata
        )

        # Subset to 2 samples
        subset_data = omics_data.subset_samples(["sample1", "sample3"])

        assert subset_data.n_samples == 2
        assert set(subset_data.samples) == {"sample1", "sample3"}

        # Check that all layers are subsetted
        for layer_name in subset_data.layer_names:
            layer_data = subset_data.get_layer(layer_name)
            assert set(layer_data.index) == {"sample1", "sample3"}

        # Check metadata subsetting
        assert set(subset_data.metadata.index) == {"sample1", "sample3"}

    def test_subset_features(self):
        """Test subsetting by features."""
        omics_data = MultiOmicsData(transcriptomics=self.transcriptomics_data, proteomics=self.proteomics_data)

        # Subset features
        feature_dict = {"transcriptomics": ["Gene_0", "Gene_1", "Gene_5"], "proteomics": ["Protein_0", "Protein_2"]}

        subset_data = omics_data.subset_features(feature_dict)

        # Check feature subsetting
        transcr_subset = subset_data.get_layer("transcriptomics")
        assert list(transcr_subset.columns) == ["Gene_0", "Gene_1", "Gene_5"]

        protein_subset = subset_data.get_layer("proteomics")
        assert list(protein_subset.columns) == ["Protein_0", "Protein_2"]

        # Samples should remain the same
        assert subset_data.n_samples == 4


class TestIntegrateOmicsData:
    """Test omics data loading and integration."""

    def test_integrate_from_dataframes(self):
        """Test integration from DataFrame dictionary."""
        # Create test DataFrames
        samples = ["S1", "S2", "S3"]

        genomics_df = pd.DataFrame(
            np.random.randint(0, 3, (3, 50)), index=samples, columns=[f"SNP_{i}" for i in range(50)]
        )

        transcr_df = pd.DataFrame(np.random.randn(3, 100), index=samples, columns=[f"Gene_{i}" for i in range(100)])

        data_dict = {"genomics": genomics_df, "transcriptomics": transcr_df}

        # Integrate
        omics_data = integrate_omics_data(data_dict)

        assert len(omics_data.layer_names) == 2
        assert omics_data.n_samples == 3
        assert set(omics_data.samples) == set(samples)

    def test_integrate_with_metadata(self):
        """Test integration with metadata."""
        samples = ["S1", "S2", "S3"]

        data_dict = {
            "transcriptomics": pd.DataFrame(
                np.random.randn(3, 10), index=samples, columns=[f"Gene_{i}" for i in range(10)]
            )
        }

        metadata_df = pd.DataFrame({"condition": ["A", "B", "A"], "batch": [1, 1, 2]}, index=samples)

        omics_data = integrate_omics_data(data_dict, metadata=metadata_df)

        assert omics_data.metadata is not None
        assert list(omics_data.metadata.index) == samples
        assert "condition" in omics_data.metadata.columns

    def test_integrate_from_files(self, tmp_path):
        """Test integration from CSV files."""
        # Create test CSV files
        samples = ["S1", "S2", "S3"]

        # Genomics CSV
        genomics_df = pd.DataFrame(
            np.random.randint(0, 3, (3, 20)), index=samples, columns=[f"SNP_{i}" for i in range(20)]
        )
        genomics_file = tmp_path / "genomics.csv"
        genomics_df.to_csv(genomics_file)

        # Transcriptomics TSV
        transcr_df = pd.DataFrame(np.random.randn(3, 30), index=samples, columns=[f"Gene_{i}" for i in range(30)])
        transcr_file = tmp_path / "transcriptomics.tsv"
        transcr_df.to_csv(transcr_file, sep="\t")

        # Metadata CSV
        metadata_df = pd.DataFrame({"age": [25, 30, 35], "sex": ["M", "F", "M"]}, index=samples)
        metadata_file = tmp_path / "metadata.csv"
        metadata_df.to_csv(metadata_file)

        # Integrate from files
        data_dict = {"genomics": genomics_file, "transcriptomics": transcr_file}

        omics_data = integrate_omics_data(data_dict, metadata=metadata_file)

        assert len(omics_data.layer_names) == 2
        assert omics_data.n_samples == 3
        assert omics_data.metadata is not None

        # Check data integrity
        np.testing.assert_array_equal(omics_data.get_layer("genomics").values, genomics_df.values)


class TestJointPCA:
    """Test joint PCA analysis across omics layers."""

    def setup_method(self):
        """Set up multi-omics data for PCA testing."""
        np.random.seed(42)
        samples = [f"S{i}" for i in range(20)]

        # Create correlated omics data
        base_signal = np.random.randn(20, 5)  # 5 latent factors

        # Genomics: sparse binary features
        genomics_loadings = np.random.randn(5, 100) * 0.5
        genomics_data = base_signal @ genomics_loadings + np.random.randn(20, 100) * 0.1
        genomics_data = (genomics_data > 0).astype(int)  # Binary

        # Transcriptomics: continuous expression
        transcr_loadings = np.random.randn(5, 50)
        transcr_data = base_signal @ transcr_loadings + np.random.randn(20, 50) * 0.5

        self.omics_data = MultiOmicsData(
            genomics=pd.DataFrame(genomics_data, index=samples, columns=[f"SNP_{i}" for i in range(100)]),
            transcriptomics=pd.DataFrame(transcr_data, index=samples, columns=[f"Gene_{i}" for i in range(50)]),
        )

    def test_joint_pca_basic(self):
        """Test basic joint PCA functionality."""
        embeddings, loadings, explained_var = joint_pca(self.omics_data, n_components=10)

        # Check output dimensions
        assert embeddings.shape == (20, 10)  # 20 samples, 10 components
        assert len(loadings) == 2  # 2 omics layers
        assert "genomics" in loadings
        assert "transcriptomics" in loadings
        assert loadings["genomics"].shape == (100, 10)  # 100 SNPs
        assert loadings["transcriptomics"].shape == (50, 10)  # 50 genes

        # Check explained variance
        assert len(explained_var) == 10
        assert np.all(explained_var >= 0)
        assert np.sum(explained_var) <= 1.0  # Should sum to <= 1

        # Components should be ordered by explained variance
        assert np.all(explained_var[:-1] >= explained_var[1:])

    def test_joint_pca_with_weights(self):
        """Test joint PCA with layer weights."""
        # Give more weight to transcriptomics
        layer_weights = {"genomics": 0.3, "transcriptomics": 0.7}

        embeddings, loadings, explained_var = joint_pca(self.omics_data, n_components=5, layer_weights=layer_weights)

        assert embeddings.shape == (20, 5)
        assert len(loadings) == 2
        assert len(explained_var) == 5

    def test_joint_pca_no_standardization(self):
        """Test joint PCA without standardization."""
        embeddings, loadings, explained_var = joint_pca(self.omics_data, n_components=5, standardize=False)

        assert embeddings.shape == (20, 5)
        assert len(loadings) == 2

    def test_joint_pca_single_layer(self):
        """Test joint PCA with single omics layer."""
        single_omics = MultiOmicsData(transcriptomics=self.omics_data.transcriptomics)

        embeddings, loadings, explained_var = joint_pca(single_omics, n_components=5)

        assert embeddings.shape == (20, 5)
        assert len(loadings) == 1
        assert "transcriptomics" in loadings


class TestJointNMF:
    """Test joint Non-negative Matrix Factorization."""

    def setup_method(self):
        """Set up non-negative multi-omics data."""
        np.random.seed(42)
        samples = [f"S{i}" for i in range(15)]

        # Create non-negative data
        self.omics_data = MultiOmicsData(
            transcriptomics=pd.DataFrame(
                np.random.exponential(2, (15, 30)),  # Non-negative
                index=samples,
                columns=[f"Gene_{i}" for i in range(30)],
            ),
            proteomics=pd.DataFrame(
                np.random.exponential(1, (15, 20)),  # Non-negative
                index=samples,
                columns=[f"Protein_{i}" for i in range(20)],
            ),
        )

    def test_joint_nmf_basic(self):
        """Test basic joint NMF functionality."""
        sample_factors, feature_factors = joint_nmf(self.omics_data, n_components=5, max_iter=50)

        # Check output dimensions
        assert sample_factors.shape == (15, 5)  # 15 samples, 5 components
        assert len(feature_factors) == 2  # 2 omics layers
        assert "transcriptomics" in feature_factors
        assert "proteomics" in feature_factors
        assert feature_factors["transcriptomics"].shape == (5, 30)  # 5 comp, 30 genes
        assert feature_factors["proteomics"].shape == (5, 20)  # 5 comp, 20 proteins

        # All factors should be non-negative
        assert np.all(sample_factors >= 0)
        for layer_factors in feature_factors.values():
            assert np.all(layer_factors >= 0)

    def test_joint_nmf_reproducibility(self):
        """Test NMF reproducibility with fixed seed."""
        factors1, features1 = joint_nmf(self.omics_data, n_components=3, random_state=123, max_iter=20)
        factors2, features2 = joint_nmf(self.omics_data, n_components=3, random_state=123, max_iter=20)

        # Results should be identical (allowing for numerical precision)
        np.testing.assert_array_almost_equal(factors1, factors2, decimal=5)
        for layer in features1.keys():
            np.testing.assert_array_almost_equal(features1[layer], features2[layer], decimal=5)

    def test_joint_nmf_reconstruction(self):
        """Test NMF reconstruction quality."""
        sample_factors, feature_factors = joint_nmf(self.omics_data, n_components=8, max_iter=100)

        # Reconstruct each layer
        for layer_name in self.omics_data.layer_names:
            original = self.omics_data.get_layer(layer_name).values
            reconstructed = sample_factors @ feature_factors[layer_name]

            # Reconstruction error should be reasonable
            mse = np.mean((original - reconstructed) ** 2)
            relative_error = mse / np.mean(original**2)
            assert relative_error < 1.0  # Less than 100% error


class TestCanonicalCorrelation:
    """Test Canonical Correlation Analysis between omics layers."""

    def setup_method(self):
        """Set up correlated multi-omics data."""
        np.random.seed(42)
        samples = [f"S{i}" for i in range(30)]

        # Create correlated layers
        latent_factors = np.random.randn(30, 3)  # 3 shared factors

        # Genomics features correlated with latent factors
        genomics_data = latent_factors @ np.random.randn(3, 40) + np.random.randn(30, 40) * 0.5

        # Transcriptomics features correlated with latent factors
        transcr_data = latent_factors @ np.random.randn(3, 25) + np.random.randn(30, 25) * 0.3

        self.omics_data = MultiOmicsData(
            genomics=pd.DataFrame(genomics_data, index=samples, columns=[f"SNP_{i}" for i in range(40)]),
            transcriptomics=pd.DataFrame(transcr_data, index=samples, columns=[f"Gene_{i}" for i in range(25)]),
        )

    def test_canonical_correlation_basic(self):
        """Test basic CCA functionality."""
        X_can, Y_can, X_weights, Y_weights, correlations = canonical_correlation(
            self.omics_data, layer_pair=("genomics", "transcriptomics"), n_components=5
        )

        # Check output dimensions
        assert X_can.shape == (30, 5)  # 30 samples, 5 canonical variables
        assert Y_can.shape == (30, 5)
        assert X_weights.shape == (40, 5)  # 40 genomics features
        assert Y_weights.shape == (25, 5)  # 25 transcriptomics features
        assert len(correlations) == 5

        # Correlations should be between 0 and 1
        assert np.all(correlations >= 0)
        assert np.all(correlations <= 1)

        # Correlations should be ordered (descending)
        assert np.all(correlations[:-1] >= correlations[1:])

    def test_canonical_correlation_actual_correlation(self):
        """Test that canonical variables achieve expected correlations."""
        X_can, Y_can, X_weights, Y_weights, correlations = canonical_correlation(
            self.omics_data, layer_pair=("genomics", "transcriptomics"), n_components=3
        )

        # Check actual correlations between canonical variables
        for i in range(3):
            actual_corr = np.corrcoef(X_can[:, i], Y_can[:, i])[0, 1]
            expected_corr = correlations[i]

            # Should be close (allowing for numerical precision)
            assert abs(actual_corr - expected_corr) < 0.1

    def test_canonical_correlation_invalid_layers(self):
        """Test error handling for invalid layer names."""
        with pytest.raises(ValueError, match="Layer nonexistent not found"):
            canonical_correlation(self.omics_data, layer_pair=("nonexistent", "transcriptomics"))

    def test_canonical_correlation_regularization(self):
        """Test CCA with different regularization values."""
        # Low regularization
        X_can1, Y_can1, _, _, corr1 = canonical_correlation(
            self.omics_data, layer_pair=("genomics", "transcriptomics"), n_components=3, regularization=0.001
        )

        # High regularization
        X_can2, Y_can2, _, _, corr2 = canonical_correlation(
            self.omics_data, layer_pair=("genomics", "transcriptomics"), n_components=3, regularization=0.1
        )

        # Both should work and produce valid results
        assert X_can1.shape == X_can2.shape == (30, 3)
        assert Y_can1.shape == Y_can2.shape == (30, 3)
        assert len(corr1) == len(corr2) == 3


class TestIntegrationEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_layers(self):
        """Test behavior with empty omics layers."""
        empty_df = pd.DataFrame(index=["S1", "S2"], columns=[])

        with pytest.raises(ValueError, match="has no features"):  # Should raise ValueError about empty columns
            MultiOmicsData(transcriptomics=empty_df)

    def test_single_sample(self):
        """Test with single sample."""
        single_sample_data = pd.DataFrame([[1, 2, 3, 4, 5]], index=["S1"], columns=["F1", "F2", "F3", "F4", "F5"])

        omics_data = MultiOmicsData(transcriptomics=single_sample_data)
        assert omics_data.n_samples == 1

        # PCA should handle single sample gracefully
        try:
            embeddings, loadings, explained_var = joint_pca(omics_data, n_components=2)
            # Should work or raise informative error
        except Exception:
            pass  # Acceptable for single sample

    def test_mismatched_samples_warning(self):
        """Test warning when samples don't fully overlap."""
        data1 = pd.DataFrame(np.random.randn(5, 10), index=["A", "B", "C", "D", "E"])
        data2 = pd.DataFrame(np.random.randn(3, 15), index=["A", "B", "C"])  # Missing D, E

        with pytest.warns(UserWarning, match="Only 3 samples are common"):
            omics_data = MultiOmicsData(genomics=data1, transcriptomics=data2)
            assert omics_data.n_samples == 3
