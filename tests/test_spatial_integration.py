"""Tests for metainformant.spatial.integration -- scRNA-seq to spatial mapping and imputation.

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.integration.scrna_mapping import (
    ImputationResult,
    MappingResult,
    anchor_based_transfer,
    correlation_mapping,
    impute_spatial_genes,
    map_scrna_to_spatial,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def scrna_data() -> np.ndarray:
    """Synthetic scRNA-seq data: 60 cells x 20 genes, 3 cell types."""
    rng = np.random.RandomState(42)
    # Create expression that differs by type
    data = np.zeros((60, 20), dtype=np.float64)
    # Type A: genes 0-6 high
    data[:20, :7] = rng.exponential(5.0, (20, 7))
    data[:20, 7:] = rng.exponential(0.5, (20, 13))
    # Type B: genes 7-13 high
    data[20:40, :7] = rng.exponential(0.5, (20, 7))
    data[20:40, 7:14] = rng.exponential(5.0, (20, 7))
    data[20:40, 14:] = rng.exponential(0.5, (20, 6))
    # Type C: genes 14-19 high
    data[40:, :14] = rng.exponential(0.5, (20, 14))
    data[40:, 14:] = rng.exponential(5.0, (20, 6))
    return data


@pytest.fixture()
def scrna_labels() -> np.ndarray:
    return np.array(["A"] * 20 + ["B"] * 20 + ["C"] * 20)


@pytest.fixture()
def spatial_data(scrna_data: np.ndarray) -> np.ndarray:
    """Synthetic spatial data: 15 spots x 20 genes resembling scRNA types."""
    rng = np.random.RandomState(99)
    # 5 spots per type, each close to the type mean
    data = np.zeros((15, 20), dtype=np.float64)
    # Type A spots
    data[:5] = scrna_data[:20].mean(axis=0) + rng.standard_normal((5, 20)) * 0.3
    # Type B spots
    data[5:10] = scrna_data[20:40].mean(axis=0) + rng.standard_normal((5, 20)) * 0.3
    # Type C spots
    data[10:] = scrna_data[40:].mean(axis=0) + rng.standard_normal((5, 20)) * 0.3
    return np.abs(data)


@pytest.fixture()
def gene_names() -> list[str]:
    return [f"gene_{i}" for i in range(20)]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestMapScrnaToSpatial:
    def test_correlation_method(
        self,
        scrna_data: np.ndarray,
        spatial_data: np.ndarray,
        scrna_labels: np.ndarray,
    ) -> None:
        result = map_scrna_to_spatial(
            scrna_data,
            spatial_data,
            method="correlation",
            scrna_labels=scrna_labels,
        )
        assert isinstance(result, MappingResult)
        assert len(result.predicted_labels) == 15
        assert result.label_probabilities.shape[0] == 15
        assert len(result.cell_type_names) == 3
        # Prediction scores should be between -1 and 1 (correlations)
        assert all(isinstance(s, (float, np.floating)) for s in result.prediction_scores)

    def test_anchor_method(
        self,
        scrna_data: np.ndarray,
        spatial_data: np.ndarray,
        scrna_labels: np.ndarray,
    ) -> None:
        result = map_scrna_to_spatial(
            scrna_data,
            spatial_data,
            method="anchor",
            scrna_labels=scrna_labels,
        )
        assert isinstance(result, MappingResult)
        assert len(result.predicted_labels) == 15
        assert result.method == "anchor_based"


class TestCorrelationMapping:
    def test_returns_mapping_result(self) -> None:
        rng = np.random.RandomState(42)
        spatial = np.abs(rng.standard_normal((10, 8)))
        profiles = np.abs(rng.standard_normal((3, 8)))
        result = correlation_mapping(spatial, profiles, cell_type_names=["A", "B", "C"])
        assert isinstance(result, MappingResult)
        assert len(result.predicted_labels) == 10
        assert result.label_probabilities.shape == (10, 3)
        # Probabilities should sum to 1 per spot
        np.testing.assert_allclose(result.label_probabilities.sum(axis=1), 1.0, atol=1e-6)


class TestAnchorBasedTransfer:
    def test_basic_operation(
        self,
        scrna_data: np.ndarray,
        spatial_data: np.ndarray,
        scrna_labels: np.ndarray,
    ) -> None:
        result = anchor_based_transfer(
            spatial_data,
            scrna_data,
            anchors=None,
            scrna_labels=scrna_labels,
            n_pcs=10,
            k_anchor=3,
            seed=42,
        )
        assert isinstance(result, MappingResult)
        assert len(result.predicted_labels) == 15
        assert result.label_probabilities.shape == (15, 3)


class TestImputeSpatialGenes:
    def test_returns_imputation_result(
        self,
        scrna_data: np.ndarray,
        spatial_data: np.ndarray,
        gene_names: list[str],
    ) -> None:
        # Spatial has only first 15 genes, impute genes 15-19
        spatial_subset = spatial_data[:, :15]
        spatial_gene_names = gene_names[:15]
        scrna_gene_names = gene_names

        result = impute_spatial_genes(
            spatial_subset,
            scrna_data,
            genes=gene_names[15:],  # impute last 5 genes
            spatial_genes=spatial_gene_names,
            scrna_genes=scrna_gene_names,
            n_neighbors=5,
            n_pcs=10,
            seed=42,
        )
        assert isinstance(result, ImputationResult)
        assert result.imputed_expression.shape == (15, 5)
        assert len(result.gene_names) == 5
        assert len(result.confidence) == 5
        assert all(0.0 <= c <= 1.0 for c in result.confidence)
        assert result.method == "knn_weighted"
