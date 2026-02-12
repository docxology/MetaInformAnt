"""Tests for metagenomics comparative submodule.

Tests differential abundance, CLR transform, indicator species, biomarkers.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

import math

import pytest

from metainformant.metagenomics.comparative.differential_abundance import (
    biomarker_discovery,
    clr_transform,
    differential_abundance,
    effect_size_analysis,
    indicator_species,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def count_matrix_two_groups() -> tuple[list[list[int]], list[int], list[str]]:
    """Count matrix with clear differential abundance between two groups.

    Group 0 (samples 0-3): taxon_A dominant.
    Group 1 (samples 4-7): taxon_B dominant.
    """
    counts = [
        # Group 0: high taxon_A, low taxon_B
        [100, 5, 20, 10],
        [90, 8, 18, 12],
        [110, 3, 22, 9],
        [95, 6, 19, 11],
        # Group 1: low taxon_A, high taxon_B
        [5, 100, 15, 12],
        [8, 95, 17, 10],
        [3, 110, 14, 13],
        [6, 105, 16, 11],
    ]
    groups = [0, 0, 0, 0, 1, 1, 1, 1]
    taxa_names = ["taxon_A", "taxon_B", "taxon_C", "taxon_D"]
    return counts, groups, taxa_names


@pytest.fixture
def simple_count_matrix() -> tuple[list[list[int]], list[int], list[str]]:
    """Minimal count matrix for basic function checks."""
    counts = [
        [10, 20, 30],
        [15, 25, 35],
        [5, 50, 10],
        [8, 45, 12],
    ]
    groups = [0, 0, 1, 1]
    taxa_names = ["sp_X", "sp_Y", "sp_Z"]
    return counts, groups, taxa_names


# ---------------------------------------------------------------------------
# Tests: differential_abundance
# ---------------------------------------------------------------------------


class TestDifferentialAbundance:
    """Tests for differential abundance analysis."""

    def test_aldex2_like_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="aldex2_like")
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "log2fc" in r
            assert "p_value" in r
            assert "adjusted_p" in r
            assert "effect_size" in r

    def test_ancom_like_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="ancom_like")
        assert len(results) == len(taxa)

    def test_simple_deseq_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="simple_deseq")
        assert len(results) == len(taxa)

    def test_sorted_by_adjusted_p(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa)
        adjusted_ps = [r["adjusted_p"] for r in results]
        assert adjusted_ps == sorted(adjusted_ps)

    def test_invalid_method_raises(
        self,
        simple_count_matrix: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = simple_count_matrix
        with pytest.raises(ValueError, match="Invalid method"):
            differential_abundance(counts, groups, taxa, method="bogus")

    def test_mismatched_groups_raises(self) -> None:
        with pytest.raises(ValueError, match="groups length"):
            differential_abundance([[1, 2]], [0, 0], ["a", "b"])

    def test_three_groups_raises(self) -> None:
        with pytest.raises(ValueError, match="exactly 2 unique"):
            differential_abundance(
                [[1, 2], [3, 4], [5, 6]],
                [0, 1, 2],
                ["a", "b"],
            )


# ---------------------------------------------------------------------------
# Tests: clr_transform
# ---------------------------------------------------------------------------


class TestClrTransform:
    """Tests for centered log-ratio transformation."""

    def test_rows_are_centered(self) -> None:
        counts = [[10, 20, 30], [5, 50, 10]]
        result = clr_transform(counts)
        assert len(result) == 2
        for row in result:
            assert len(row) == 3
            # CLR rows should sum to approximately 0
            assert abs(sum(row)) < 1e-6

    def test_pseudocount_avoids_log_zero(self) -> None:
        counts = [[0, 0, 100]]
        result = clr_transform(counts, pseudocount=0.5)
        assert len(result) == 1
        # No NaN or Inf values
        for v in result[0]:
            assert math.isfinite(v)

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            clr_transform([])


# ---------------------------------------------------------------------------
# Tests: indicator_species
# ---------------------------------------------------------------------------


class TestIndicatorSpecies:
    """Tests for indicator species analysis."""

    def test_returns_species_group_associations(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = indicator_species(counts, groups, taxa, n_permutations=99, seed=42)
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "indicator_value" in r
            assert "p_value" in r
            assert "associated_group" in r
            assert 0.0 <= r["indicator_value"] <= 1.0
            assert 0.0 <= r["p_value"] <= 1.0

    def test_sorted_by_indicator_value(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = indicator_species(counts, groups, taxa, n_permutations=49, seed=42)
        indvals = [r["indicator_value"] for r in results]
        assert indvals == sorted(indvals, reverse=True)


# ---------------------------------------------------------------------------
# Tests: effect_size_analysis
# ---------------------------------------------------------------------------


class TestEffectSizeAnalysis:
    """Tests for LEfSe-style effect size ranking."""

    def test_returns_lda_scores(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = effect_size_analysis(counts, groups, taxa)
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "lda_score" in r
            assert "p_value" in r
            assert "direction" in r

    def test_sorted_by_abs_lda(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = effect_size_analysis(counts, groups, taxa)
        lda_scores = [abs(r["lda_score"]) for r in results]
        assert lda_scores == sorted(lda_scores, reverse=True)


# ---------------------------------------------------------------------------
# Tests: biomarker_discovery
# ---------------------------------------------------------------------------


class TestBiomarkerDiscovery:
    """Tests for ML-based biomarker discovery."""

    def test_returns_expected_keys(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        result = biomarker_discovery(counts, groups, taxa, seed=42)
        assert isinstance(result, dict)
        assert "selected_taxa" in result
        assert "importances" in result
        assert "cv_accuracy" in result
        assert "model_summary" in result
        assert isinstance(result["selected_taxa"], list)
        assert isinstance(result["importances"], dict)
        assert all(t in taxa for t in result["importances"])

    def test_importances_are_nonneg(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        result = biomarker_discovery(counts, groups, taxa, seed=42)
        for imp in result["importances"].values():
            assert imp >= 0.0

    def test_invalid_method_raises(
        self,
        simple_count_matrix: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = simple_count_matrix
        with pytest.raises(ValueError, match="Invalid method"):
            biomarker_discovery(counts, groups, taxa, method="xgboost")
