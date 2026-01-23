"""Tests for RNA-Protein integration functions.

This module tests protein integration functionality following NO_MOCKING_POLICY.
All tests use real data structures and calculations.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna import protein_integration


class TestCalculateTranslationEfficiency:
    """Test calculate_translation_efficiency function."""

    def test_calculate_translation_efficiency_empty(self):
        """Test calculate_translation_efficiency with empty DataFrames."""
        rna_df = pd.DataFrame()
        protein_df = pd.DataFrame()

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_calculate_translation_efficiency_no_common_samples(self):
        """Test calculate_translation_efficiency with no common samples."""
        rna_df = pd.DataFrame({"gene1": [1, 2, 3]}, index=["sample1", "sample2", "sample3"])
        protein_df = pd.DataFrame({"gene1": [1, 2, 3]}, index=["sample4", "sample5", "sample6"])

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_calculate_translation_efficiency_ratio_method(self):
        """Test calculate_translation_efficiency with ratio method."""
        rna_df = pd.DataFrame(
            {"gene1": [10, 20, 30], "gene2": [5, 10, 15]},
            index=["sample1", "sample2", "sample3"],
        )
        protein_df = pd.DataFrame(
            {"gene1": [1, 2, 3], "gene2": [0.5, 1, 1.5]},
            index=["sample1", "sample2", "sample3"],
        )

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "gene_id" in result.columns
        assert "efficiency" in result.columns
        assert "method" in result.columns

    def test_calculate_translation_efficiency_correlation_method(self):
        """Test calculate_translation_efficiency with correlation method."""
        rna_df = pd.DataFrame(
            {"gene1": [10, 20, 30, 40, 50]},
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
        )
        protein_df = pd.DataFrame(
            {"gene1": [1, 2, 3, 4, 5]},
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
        )

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df, method="correlation")
        assert isinstance(result, pd.DataFrame)
        if len(result) > 0:
            assert "gene_id" in result.columns
            assert "efficiency" in result.columns


class TestPredictProteinAbundanceFromRna:
    """Test predict_protein_abundance_from_rna function."""

    def test_predict_protein_abundance_empty(self):
        """Test predict_protein_abundance_from_rna with empty DataFrame."""
        rna_df = pd.DataFrame()

        result = protein_integration.predict_protein_abundance_from_rna(rna_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_predict_protein_abundance_simple(self):
        """Test predict_protein_abundance_from_rna with simple prediction."""
        rna_df = pd.DataFrame(
            {"gene1": [10, 20, 30], "gene2": [5, 10, 15]},
            index=["sample1", "sample2", "sample3"],
        )

        result = protein_integration.predict_protein_abundance_from_rna(rna_df, method="linear")
        assert isinstance(result, pd.DataFrame)
        assert result.shape == rna_df.shape
        assert list(result.columns) == list(rna_df.columns)
        assert list(result.index) == list(rna_df.index)

    def test_predict_protein_abundance_with_training(self):
        """Test predict_protein_abundance_from_rna with training data."""
        training_rna = pd.DataFrame(
            {"gene1": [10, 20, 30]},
            index=["train1", "train2", "train3"],
        )
        training_protein = pd.DataFrame(
            {"gene1": [1, 2, 3]},
            index=["train1", "train2", "train3"],
        )
        rna_df = pd.DataFrame(
            {"gene1": [15, 25]},
            index=["test1", "test2"],
        )

        result = protein_integration.predict_protein_abundance_from_rna(
            rna_df,
            training_rna=training_rna,
            training_protein=training_protein,
            method="linear",
        )
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == rna_df.shape[0]


class TestRibosomeProfilingIntegration:
    """Test ribosome_profiling_integration function."""

    def test_ribosome_profiling_integration_empty(self):
        """Test ribosome_profiling_integration with empty DataFrames."""
        rna_df = pd.DataFrame()
        ribo_df = pd.DataFrame()

        result = protein_integration.ribosome_profiling_integration(rna_df, ribo_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_ribosome_profiling_integration_basic(self):
        """Test ribosome_profiling_integration with basic data."""
        rna_df = pd.DataFrame(
            {"gene1": [10, 20, 30], "gene2": [5, 10, 15]},
            index=["sample1", "sample2", "sample3"],
        )
        ribo_df = pd.DataFrame(
            {"gene1": [1, 2, 3], "gene2": [0.5, 1, 1.5]},
            index=["sample1", "sample2", "sample3"],
        )

        result = protein_integration.ribosome_profiling_integration(rna_df, ribo_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "gene_id" in result.columns
        assert "rna_level" in result.columns
        assert "ribo_level" in result.columns
        assert "translation_rate" in result.columns
        assert "translationally_regulated" in result.columns

    def test_ribosome_profiling_integration_no_common(self):
        """Test ribosome_profiling_integration with no common samples."""
        rna_df = pd.DataFrame(
            {"gene1": [10, 20, 30]},
            index=["sample1", "sample2", "sample3"],
        )
        ribo_df = pd.DataFrame(
            {"gene1": [1, 2, 3]},
            index=["sample4", "sample5", "sample6"],
        )

        result = protein_integration.ribosome_profiling_integration(rna_df, ribo_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


class TestProteinIntegrationDocumentation:
    """Test that protein integration functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all protein integration functions have docstrings."""
        functions = [
            protein_integration.calculate_translation_efficiency,
            protein_integration.predict_protein_abundance_from_rna,
            protein_integration.ribosome_profiling_integration,
        ]

        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0


class TestProteinIntegrationEdgeCases:
    """Test edge cases in protein integration functions."""

    def test_calculate_translation_efficiency_handles_zeros(self):
        """Test calculate_translation_efficiency handles zero values."""
        rna_df = pd.DataFrame(
            {"gene1": [0, 10, 20]},
            index=["sample1", "sample2", "sample3"],
        )
        protein_df = pd.DataFrame(
            {"gene1": [0, 1, 2]},
            index=["sample1", "sample2", "sample3"],
        )

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        # Should filter out zeros
        assert isinstance(result, pd.DataFrame)

    def test_calculate_translation_efficiency_handles_nans(self):
        """Test calculate_translation_efficiency handles NaN values."""
        rna_df = pd.DataFrame(
            {"gene1": [np.nan, 10, 20]},
            index=["sample1", "sample2", "sample3"],
        )
        protein_df = pd.DataFrame(
            {"gene1": [1, np.nan, 2]},
            index=["sample1", "sample2", "sample3"],
        )

        result = protein_integration.calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        # Should filter out NaNs
        assert isinstance(result, pd.DataFrame)
