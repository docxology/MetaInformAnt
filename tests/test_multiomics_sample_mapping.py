"""Tests for multiomics sample mapping functionality."""

import pandas as pd
import pytest

from metainformant.multiomics.analysis.integration import integrate_omics_data


class TestSampleMapping:
    """Test sample ID mapping functionality."""

    def test_integrate_omics_with_sample_mapping(self):
        """Test integration with sample ID mapping."""
        # Create test data with different sample IDs
        genomics_data = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0]},
            index=["sample_A", "sample_B", "sample_C"],
        )

        transcriptomics_data = pd.DataFrame(
            {"gene1": [4.0, 5.0, 6.0]},
            index=["sample_1", "sample_2", "sample_3"],
        )

        # Define sample mapping
        sample_mapping = {
            "genomics": {"sample_A": "sample_1", "sample_B": "sample_2", "sample_C": "sample_3"},
        }

        # Integrate with mapping
        integrated = integrate_omics_data(
            {"genomics": genomics_data, "transcriptomics": transcriptomics_data},
            sample_mapping=sample_mapping,
        )

        # Check that samples are properly aligned
        assert len(integrated.samples) == 3
        assert "sample_1" in integrated.samples
        assert "sample_2" in integrated.samples
        assert "sample_3" in integrated.samples

        # Check that data is accessible
        genomics_layer = integrated.get_layer("genomics")
        transcriptomics_layer = integrated.get_layer("transcriptomics")

        assert genomics_layer.shape == (3, 1)
        assert transcriptomics_layer.shape == (3, 1)

    def test_integrate_omics_with_feature_mapping(self):
        """Test integration with feature ID mapping."""
        genomics_data = pd.DataFrame(
            {"GENE1": [1.0, 2.0, 3.0]},
            index=["sample1", "sample2", "sample3"],
        )

        transcriptomics_data = pd.DataFrame(
            {"gene_001": [4.0, 5.0, 6.0]},
            index=["sample1", "sample2", "sample3"],
        )

        # Define feature mapping
        feature_mapping = {
            "genomics": {"GENE1": "gene_001"},
        }

        # Integrate with mapping
        integrated = integrate_omics_data(
            {"genomics": genomics_data, "transcriptomics": transcriptomics_data},
            feature_mapping=feature_mapping,
        )

        # Check that features are properly mapped
        genomics_layer = integrated.get_layer("genomics")
        transcriptomics_layer = integrated.get_layer("transcriptomics")

        assert "gene_001" in genomics_layer.columns
        assert "gene_001" in transcriptomics_layer.columns

    def test_integrate_omics_with_both_mappings(self):
        """Test integration with both sample and feature mapping."""
        genomics_data = pd.DataFrame(
            {"GENE1": [1.0, 2.0, 3.0]},
            index=["sample_A", "sample_B", "sample_C"],
        )

        transcriptomics_data = pd.DataFrame(
            {"gene_001": [4.0, 5.0, 6.0]},
            index=["sample_1", "sample_2", "sample_3"],
        )

        # Define both mappings
        sample_mapping = {
            "genomics": {"sample_A": "sample_1", "sample_B": "sample_2", "sample_C": "sample_3"},
        }

        feature_mapping = {
            "genomics": {"GENE1": "gene_001"},
        }

        # Integrate with both mappings
        integrated = integrate_omics_data(
            {"genomics": genomics_data, "transcriptomics": transcriptomics_data},
            sample_mapping=sample_mapping,
            feature_mapping=feature_mapping,
        )

        # Check that both mappings are applied
        genomics_layer = integrated.get_layer("genomics")
        transcriptomics_layer = integrated.get_layer("transcriptomics")

        assert "gene_001" in genomics_layer.columns
        assert "gene_001" in transcriptomics_layer.columns
        assert len(integrated.samples) == 3

    def test_integrate_omics_with_metadata_mapping(self):
        """Test integration with metadata and sample mapping."""
        genomics_data = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0]},
            index=["sample_A", "sample_B", "sample_C"],
        )

        metadata = pd.DataFrame(
            {"condition": ["control", "treated", "control"]},
            index=["sample_1", "sample_2", "sample_3"],
        )

        # Define sample mapping
        sample_mapping = {
            "genomics": {"sample_A": "sample_1", "sample_B": "sample_2", "sample_C": "sample_3"},
        }

        # Integrate with mapping and metadata
        integrated = integrate_omics_data(
            {"genomics": genomics_data},
            sample_mapping=sample_mapping,
            metadata=metadata,
        )

        # Check that metadata is properly mapped
        assert integrated.metadata is not None
        assert len(integrated.metadata) == 3
        assert "condition" in integrated.metadata.columns

    def test_integrate_omics_mapping_preserves_data(self):
        """Test that sample mapping preserves data integrity."""
        genomics_data = pd.DataFrame(
            {"gene1": [1.0, 2.0, 3.0], "gene2": [4.0, 5.0, 6.0]},
            index=["sample_A", "sample_B", "sample_C"],
        )

        # Define sample mapping
        sample_mapping = {
            "genomics": {"sample_A": "sample_1", "sample_B": "sample_2", "sample_C": "sample_3"},
        }

        # Integrate with mapping
        integrated = integrate_omics_data(
            {"genomics": genomics_data},
            sample_mapping=sample_mapping,
        )

        # Check that data values are preserved
        genomics_layer = integrated.get_layer("genomics")
        assert genomics_layer.loc["sample_1", "gene1"] == 1.0
        assert genomics_layer.loc["sample_2", "gene2"] == 5.0
        assert genomics_layer.loc["sample_3", "gene1"] == 3.0
