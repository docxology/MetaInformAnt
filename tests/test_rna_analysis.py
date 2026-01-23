"""Tests for RNA analysis module functions.

Tests protein integration, translation efficiency, and ribosome profiling
functions following the NO_MOCKING policy.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.protein_integration import (
    calculate_translation_efficiency,
    predict_protein_abundance_from_rna,
    ribosome_profiling_integration,
)


class TestCalculateTranslationEfficiency:
    """Tests for translation efficiency calculation."""

    def test_ratio_method_basic(self) -> None:
        """Test ratio method with simple data."""
        rna_df = pd.DataFrame(
            {"gene1": [10.0, 20.0, 30.0], "gene2": [5.0, 10.0, 15.0]},
            index=["sample1", "sample2", "sample3"]
        )
        protein_df = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0], "gene2": [0.5, 1.0, 1.5]},
            index=["sample1", "sample2", "sample3"]
        )

        result = calculate_translation_efficiency(rna_df, protein_df, method="ratio")

        assert len(result) == 2
        assert "gene_id" in result.columns
        assert "efficiency" in result.columns
        assert "method" in result.columns
        assert all(result["method"] == "ratio")

        # Check efficiency values are reasonable (protein/RNA ratio)
        gene1_eff = result[result["gene_id"] == "gene1"]["efficiency"].values[0]
        assert 0.09 < gene1_eff < 0.11  # Should be ~0.1 (1/10)

    def test_correlation_method(self) -> None:
        """Test correlation method with correlated data."""
        # Create perfectly correlated data
        rna_df = pd.DataFrame(
            {"gene1": [10.0, 20.0, 30.0, 40.0], "gene2": [5.0, 10.0, 15.0, 20.0]},
            index=["s1", "s2", "s3", "s4"]
        )
        protein_df = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0, 4.0], "gene2": [2.5, 5.0, 7.5, 10.0]},
            index=["s1", "s2", "s3", "s4"]
        )

        result = calculate_translation_efficiency(rna_df, protein_df, method="correlation")

        assert len(result) == 2
        # Perfect correlation should give efficiency close to 1.0
        for _, row in result.iterrows():
            assert row["efficiency"] > 0.99

    def test_empty_dataframes(self) -> None:
        """Test handling of empty DataFrames."""
        result = calculate_translation_efficiency(
            pd.DataFrame(),
            pd.DataFrame()
        )
        assert len(result) == 0
        assert list(result.columns) == ["gene_id", "efficiency", "method"]

    def test_no_common_samples(self) -> None:
        """Test when no samples overlap between datasets."""
        rna_df = pd.DataFrame({"gene1": [1, 2]}, index=["s1", "s2"])
        protein_df = pd.DataFrame({"gene1": [1, 2]}, index=["s3", "s4"])

        result = calculate_translation_efficiency(rna_df, protein_df)
        assert len(result) == 0

    def test_no_common_genes(self) -> None:
        """Test when no genes overlap between datasets."""
        rna_df = pd.DataFrame({"gene1": [1, 2]}, index=["s1", "s2"])
        protein_df = pd.DataFrame({"gene2": [1, 2]}, index=["s1", "s2"])

        result = calculate_translation_efficiency(rna_df, protein_df)
        assert len(result) == 0

    def test_handles_zero_values(self) -> None:
        """Test that zero RNA values are handled correctly."""
        rna_df = pd.DataFrame(
            {"gene1": [10.0, 0.0, 30.0], "gene2": [0.0, 10.0, 15.0]},
            index=["s1", "s2", "s3"]
        )
        protein_df = pd.DataFrame(
            {"gene1": [1.0, 0.0, 3.0], "gene2": [0.0, 1.0, 1.5]},
            index=["s1", "s2", "s3"]
        )

        result = calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        # Should still compute efficiency for valid pairs
        assert len(result) >= 1


class TestPredictProteinAbundance:
    """Tests for protein abundance prediction."""

    def test_basic_prediction(self) -> None:
        """Test basic protein abundance prediction."""
        rna_df = pd.DataFrame(
            {"gene1": [10.0, 20.0, 30.0], "gene2": [5.0, 10.0, 15.0]},
            index=["s1", "s2", "s3"]
        )

        predictions = predict_protein_abundance_from_rna(rna_df, method="linear")

        assert predictions.shape == rna_df.shape
        assert list(predictions.index) == list(rna_df.index)
        assert list(predictions.columns) == list(rna_df.columns)

    def test_with_training_data(self) -> None:
        """Test prediction with training data."""
        # Training data
        train_rna = pd.DataFrame(
            {"gene1": [10.0, 20.0], "gene2": [5.0, 10.0]},
            index=["train1", "train2"]
        )
        train_protein = pd.DataFrame(
            {"gene1": [1.0, 2.0], "gene2": [2.5, 5.0]},
            index=["train1", "train2"]
        )

        # Test data
        test_rna = pd.DataFrame(
            {"gene1": [30.0, 40.0], "gene2": [15.0, 20.0]},
            index=["test1", "test2"]
        )

        predictions = predict_protein_abundance_from_rna(
            test_rna,
            training_rna=train_rna,
            training_protein=train_protein
        )

        assert predictions.shape == test_rna.shape
        # gene1 should scale by 0.1 (protein/RNA ratio from training)
        # gene2 should scale by 0.5 (protein/RNA ratio from training)
        assert predictions.loc["test1", "gene1"] < test_rna.loc["test1", "gene1"]

    def test_empty_input(self) -> None:
        """Test with empty DataFrame."""
        result = predict_protein_abundance_from_rna(pd.DataFrame())
        assert result.empty


class TestRibosomeProfilingIntegration:
    """Tests for ribosome profiling integration."""

    def test_basic_integration(self) -> None:
        """Test basic ribosome profiling integration."""
        rna_df = pd.DataFrame(
            {"gene1": [10.0, 20.0, 30.0], "gene2": [5.0, 10.0, 15.0], "gene3": [100.0, 200.0, 300.0]},
            index=["s1", "s2", "s3"]
        )
        ribo_df = pd.DataFrame(
            {"gene1": [5.0, 10.0, 15.0], "gene2": [10.0, 20.0, 30.0], "gene3": [10.0, 20.0, 30.0]},
            index=["s1", "s2", "s3"]
        )

        result = ribosome_profiling_integration(rna_df, ribo_df)

        assert len(result) == 3
        assert "gene_id" in result.columns
        assert "rna_level" in result.columns
        assert "ribo_level" in result.columns
        assert "translation_rate" in result.columns
        assert "translationally_regulated" in result.columns

        # gene1: translation rate = 0.5 (normal)
        # gene2: translation rate = 2.0 (high - translationally upregulated)
        # gene3: translation rate = 0.1 (low - translationally downregulated)

    def test_translation_rate_calculation(self) -> None:
        """Test that translation rate is calculated correctly."""
        rna_df = pd.DataFrame({"gene1": [100.0]}, index=["s1"])
        ribo_df = pd.DataFrame({"gene1": [50.0]}, index=["s1"])

        result = ribosome_profiling_integration(rna_df, ribo_df)

        gene1_rate = result[result["gene_id"] == "gene1"]["translation_rate"].values[0]
        assert abs(gene1_rate - 0.5) < 0.01  # Should be 50/100 = 0.5

    def test_zero_rna_handling(self) -> None:
        """Test handling of zero RNA values."""
        rna_df = pd.DataFrame(
            {"gene1": [0.0, 10.0], "gene2": [10.0, 0.0]},
            index=["s1", "s2"]
        )
        ribo_df = pd.DataFrame(
            {"gene1": [5.0, 10.0], "gene2": [5.0, 5.0]},
            index=["s1", "s2"]
        )

        result = ribosome_profiling_integration(rna_df, ribo_df)

        # Should handle gracefully
        assert len(result) == 2

    def test_empty_dataframes(self) -> None:
        """Test with empty DataFrames."""
        result = ribosome_profiling_integration(pd.DataFrame(), pd.DataFrame())
        assert len(result) == 0
        assert "translationally_regulated" in result.columns

    def test_translational_regulation_detection(self) -> None:
        """Test that translationally regulated genes are detected."""
        # Create data with larger set for meaningful z-score calculation
        # Need at least 10+ genes for reliable statistics
        n_normal_genes = 20

        # Normal genes with similar translation rates (ribo/rna ~ 0.1)
        rna_data = {f"gene{i}": [10.0, 20.0, 30.0] for i in range(n_normal_genes)}
        ribo_data = {f"gene{i}": [1.0, 2.0, 3.0] for i in range(n_normal_genes)}

        # Add outlier gene with very high translation rate (ribo/rna ~ 10.0)
        rna_data["outlier_high"] = [10.0, 20.0, 30.0]
        ribo_data["outlier_high"] = [100.0, 200.0, 300.0]

        # Add outlier gene with very low translation rate (ribo/rna ~ 0.01)
        rna_data["outlier_low"] = [100.0, 200.0, 300.0]
        ribo_data["outlier_low"] = [1.0, 2.0, 3.0]

        rna_df = pd.DataFrame(rna_data, index=["s1", "s2", "s3"])
        ribo_df = pd.DataFrame(ribo_data, index=["s1", "s2", "s3"])

        result = ribosome_profiling_integration(rna_df, ribo_df)

        # Verify the outlier has significantly different translation rate
        outlier_high = result[result["gene_id"] == "outlier_high"]
        outlier_low = result[result["gene_id"] == "outlier_low"]
        normal = result[result["gene_id"] == "gene0"]

        # Check translation rates are computed correctly
        assert outlier_high["translation_rate"].values[0] > 5.0  # High rate
        assert outlier_low["translation_rate"].values[0] < 0.05  # Low rate
        assert 0.05 < normal["translation_rate"].values[0] < 0.5  # Normal rate

        # Check that at least one outlier is detected
        regulated_count = result["translationally_regulated"].sum()
        # With z-score > 2, the extreme outliers should be detected
        assert regulated_count >= 1, "Expected at least one translationally regulated gene"


class TestEdgeCases:
    """Edge case tests for RNA analysis functions."""

    def test_nan_handling_in_efficiency(self) -> None:
        """Test NaN handling in translation efficiency calculation."""
        rna_df = pd.DataFrame(
            {"gene1": [10.0, np.nan, 30.0], "gene2": [5.0, 10.0, np.nan]},
            index=["s1", "s2", "s3"]
        )
        protein_df = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0], "gene2": [0.5, np.nan, 1.5]},
            index=["s1", "s2", "s3"]
        )

        result = calculate_translation_efficiency(rna_df, protein_df)
        # Should still produce results for valid pairs
        assert len(result) >= 1

    def test_single_sample(self) -> None:
        """Test with single sample (edge case for correlation)."""
        rna_df = pd.DataFrame({"gene1": [10.0]}, index=["s1"])
        protein_df = pd.DataFrame({"gene1": [1.0]}, index=["s1"])

        # Ratio should work
        result = calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        assert len(result) == 1

        # Correlation requires at least 2 samples
        result = calculate_translation_efficiency(rna_df, protein_df, method="correlation")
        assert len(result) == 0

    def test_many_genes(self) -> None:
        """Test performance with many genes."""
        n_genes = 100
        n_samples = 10

        rna_data = {f"gene{i}": np.random.rand(n_samples) * 100 for i in range(n_genes)}
        protein_data = {f"gene{i}": np.random.rand(n_samples) * 10 for i in range(n_genes)}

        rna_df = pd.DataFrame(rna_data, index=[f"s{i}" for i in range(n_samples)])
        protein_df = pd.DataFrame(protein_data, index=[f"s{i}" for i in range(n_samples)])

        result = calculate_translation_efficiency(rna_df, protein_df)
        assert len(result) == n_genes
