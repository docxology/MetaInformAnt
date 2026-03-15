"""Tests for ecology functional diversity module.

Covers:
    - functional_richness (FRic) in 1D, 2D, and 3D+
    - functional_evenness (FEve) via MST
    - functional_divergence (FDiv) centroid-based
    - functional_dispersion (FDis) weighted centroid distance
    - raos_quadratic_entropy (Rao's Q)
    - community_weighted_mean (CWM)
    - functional_redundancy (Simpson - Rao's Q)
    - trait_distance_matrix (Euclidean, Gower)
    - functional_beta_diversity (turnover and nestedness decomposition)
    - functional_diversity_suite (convenience wrapper)

Uses real implementations only (NO mocking per project policy).
"""

from __future__ import annotations

import math
from typing import List

import pytest

from metainformant.ecology.analysis.functional import (
    community_weighted_mean,
    functional_beta_diversity,
    functional_dispersion,
    functional_divergence,
    functional_diversity_suite,
    functional_evenness,
    functional_redundancy,
    functional_richness,
    raos_quadratic_entropy,
    trait_distance_matrix,
)

# ---------------------------------------------------------------------------
# Shared test data fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def trait_1d() -> List[List[float]]:
    """Five species, one trait each."""
    return [[0], [2], [5], [8], [10]]


@pytest.fixture
def trait_2d() -> List[List[float]]:
    """Five species, two traits each."""
    return [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [0.5, 0.5],
    ]


@pytest.fixture
def trait_3d() -> List[List[float]]:
    """Four species, three traits each (for hyper-volume FRic)."""
    return [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 4.0],
    ]


@pytest.fixture
def equal_abundances_5() -> List[float]:
    """Perfectly even abundance vector (5 species)."""
    return [10, 10, 10, 10, 10]


@pytest.fixture
def unequal_abundances_5() -> List[float]:
    """Uneven abundance vector (5 species)."""
    return [50, 20, 15, 10, 5]


@pytest.fixture
def equal_abundances_4() -> List[float]:
    """Perfectly even abundance vector (4 species)."""
    return [10, 10, 10, 10]


# ============================================================================
# Trait distance matrix tests
# ============================================================================


class TestTraitDistanceMatrix:
    """Test trait distance matrix computation."""

    def test_euclidean_known_value(self) -> None:
        """Euclidean distance between (0,0) and (3,4) = 5."""
        dm = trait_distance_matrix([[0, 0], [3, 4]], method="euclidean")
        assert abs(dm[0][1] - 5.0) < 1e-10
        assert abs(dm[1][0] - 5.0) < 1e-10
        assert dm[0][0] == 0.0

    def test_euclidean_symmetric(self, trait_2d: List[List[float]]) -> None:
        """Euclidean distance matrix must be symmetric with zero diagonal."""
        dm = trait_distance_matrix(trait_2d, method="euclidean")
        n = len(trait_2d)
        for i in range(n):
            assert dm[i][i] == 0.0
            for j in range(n):
                assert abs(dm[i][j] - dm[j][i]) < 1e-12

    def test_gower_range_normalized(self) -> None:
        """Gower distance normalises by range, so identical values give 0."""
        dm = trait_distance_matrix(
            [[0, 0], [10, 10], [5, 5]],
            method="gower",
        )
        # Same proportional position -> Gower(0,2) should be 0.5
        assert abs(dm[0][2] - 0.5) < 1e-10

    def test_gower_identical_species(self) -> None:
        """Gower distance between identical trait vectors should be 0."""
        dm = trait_distance_matrix([[1, 2, 3], [1, 2, 3]], method="gower")
        assert abs(dm[0][1]) < 1e-10

    def test_unsupported_method_raises(self) -> None:
        """Unknown method raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported distance method"):
            trait_distance_matrix([[1, 2], [3, 4]], method="hamming")

    def test_empty_raises(self) -> None:
        """Empty trait matrix raises ValueError."""
        with pytest.raises((ValueError, Exception)):
            trait_distance_matrix([], method="euclidean")


# ============================================================================
# Functional Richness (FRic) tests
# ============================================================================


class TestFunctionalRichness:
    """Test Functional Richness (FRic)."""

    def test_fric_1d_is_range(self, trait_1d: List[List[float]]) -> None:
        """In 1D, FRic equals the range of trait values."""
        fric = functional_richness(trait_1d)
        # Traits: 0, 2, 5, 8, 10 -> range = 10
        assert abs(fric - 10.0) < 1e-10

    def test_fric_2d_convex_hull(self) -> None:
        """In 2D, FRic equals the convex hull area."""
        # Right triangle with legs of length 1 -> area = 0.5
        traits = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        fric = functional_richness(traits)
        assert abs(fric - 0.5) < 1e-10

    def test_fric_2d_unit_square(self, trait_2d: List[List[float]]) -> None:
        """Unit square corners + midpoint -> convex hull area = 1.0."""
        fric = functional_richness(trait_2d)
        assert abs(fric - 1.0) < 1e-10

    def test_fric_3d_bounding_box(self, trait_3d: List[List[float]]) -> None:
        """In 3D, FRic approximates the bounding-box volume = 2 * 3 * 4 = 24."""
        fric = functional_richness(trait_3d)
        assert abs(fric - 24.0) < 1e-10

    def test_fric_filters_zero_abundance(self) -> None:
        """Species with zero abundance are excluded."""
        traits = [[0], [5], [10]]
        fric = functional_richness(traits, abundances=[1, 0, 1])
        # Only species 0 and 2 remain -> range = 10
        assert abs(fric - 10.0) < 1e-10

    def test_fric_single_species_is_zero(self) -> None:
        """A single species has zero FRic."""
        assert functional_richness([[1.0, 2.0]]) == 0.0

    def test_fric_empty_raises(self) -> None:
        """Empty trait matrix raises ValueError."""
        with pytest.raises((ValueError, Exception)):
            functional_richness([])

    def test_fric_abundance_mismatch_raises(self) -> None:
        """Mismatched abundances length raises ValueError."""
        with pytest.raises(ValueError, match="abundances length"):
            functional_richness([[1], [2], [3]], abundances=[1, 2])


# ============================================================================
# Functional Evenness (FEve) tests
# ============================================================================


class TestFunctionalEvenness:
    """Test Functional Evenness (FEve)."""

    def test_feve_perfectly_even(self) -> None:
        """Equal abundances with uniform spacing yield FEve = 1.0."""
        traits = [[0], [5], [10]]
        feve = functional_evenness(traits, [10, 10, 10])
        assert abs(feve - 1.0) < 0.05  # Slight tolerance for MST-based calculation

    def test_feve_range_0_to_1(
        self,
        trait_2d: List[List[float]],
        unequal_abundances_5: List[float],
    ) -> None:
        """FEve is always in [0, 1]."""
        feve = functional_evenness(trait_2d, unequal_abundances_5)
        assert 0.0 <= feve <= 1.0

    def test_feve_highly_uneven(self) -> None:
        """Highly uneven abundances should yield lower FEve."""
        traits = [[0], [5], [10]]
        feve_even = functional_evenness(traits, [10, 10, 10])
        feve_uneven = functional_evenness(traits, [100, 1, 1])
        assert feve_even >= feve_uneven - 0.01

    def test_feve_single_species_zero(self) -> None:
        """Fewer than 2 present species gives FEve = 0."""
        feve = functional_evenness([[1], [2], [3]], [0, 0, 5])
        assert feve == 0.0

    def test_feve_mismatch_raises(self) -> None:
        """Mismatched abundances raises ValueError."""
        with pytest.raises(ValueError, match="abundances length"):
            functional_evenness([[1], [2]], [1])


# ============================================================================
# Functional Divergence (FDiv) tests
# ============================================================================


class TestFunctionalDivergence:
    """Test Functional Divergence (FDiv)."""

    def test_fdiv_range_0_to_1(
        self,
        trait_2d: List[List[float]],
        equal_abundances_5: List[float],
    ) -> None:
        """FDiv is in [0, 1]."""
        fdiv = functional_divergence(trait_2d, equal_abundances_5)
        assert 0.0 <= fdiv <= 1.0

    def test_fdiv_extreme_abundance_at_edges(self) -> None:
        """When most abundance is at extreme traits, FDiv should be high."""
        traits = [[0], [5], [10]]
        fdiv_extreme = functional_divergence(traits, [50, 1, 50])
        fdiv_center = functional_divergence(traits, [1, 50, 1])
        assert fdiv_extreme > fdiv_center

    def test_fdiv_single_species_zero(self) -> None:
        """Fewer than 2 present species gives FDiv = 0."""
        fdiv = functional_divergence([[1, 2], [3, 4]], [0, 5])
        assert fdiv == 0.0

    def test_fdiv_mismatch_raises(self) -> None:
        """Mismatched abundances raises ValueError."""
        with pytest.raises(ValueError, match="abundances length"):
            functional_divergence([[1], [2], [3]], [1, 2])


# ============================================================================
# Functional Dispersion (FDis) tests
# ============================================================================


class TestFunctionalDispersion:
    """Test Functional Dispersion (FDis)."""

    def test_fdis_known_value(self) -> None:
        """Two equal-abundance species at 0 and 10: centroid=5, FDis=5."""
        fdis = functional_dispersion([[0], [10]], [1, 1])
        assert abs(fdis - 5.0) < 1e-10

    def test_fdis_zero_for_identical_traits(self) -> None:
        """Identical trait values give FDis = 0."""
        fdis = functional_dispersion([[5], [5], [5]], [1, 1, 1])
        assert abs(fdis) < 1e-10

    def test_fdis_nonnegative(
        self,
        trait_2d: List[List[float]],
        unequal_abundances_5: List[float],
    ) -> None:
        """FDis is always non-negative."""
        fdis = functional_dispersion(trait_2d, unequal_abundances_5)
        assert fdis >= -1e-10

    def test_fdis_single_species_zero(self) -> None:
        """Fewer than 2 present species gives FDis = 0."""
        fdis = functional_dispersion([[1], [2]], [0, 5])
        assert fdis == 0.0


# ============================================================================
# Rao's Quadratic Entropy tests
# ============================================================================


class TestRaosQuadraticEntropy:
    """Test Rao's Quadratic Entropy (Q)."""

    def test_raos_q_known_value(self) -> None:
        """Two equal species at trait distance 10: Q = 10 * 0.5 * 0.5 * 2 = 5."""
        q = raos_quadratic_entropy([[0], [10]], [1, 1])
        assert abs(q - 5.0) < 1e-10

    def test_raos_q_zero_for_one_species(self) -> None:
        """Single-species community: Q = 0."""
        q = raos_quadratic_entropy([[5]], [10])
        assert q == 0.0

    def test_raos_q_nonnegative(
        self,
        trait_2d: List[List[float]],
        equal_abundances_5: List[float],
    ) -> None:
        """Rao's Q is always non-negative."""
        q = raos_quadratic_entropy(trait_2d, equal_abundances_5)
        assert q >= -1e-10

    def test_raos_q_increases_with_trait_divergence(self) -> None:
        """More divergent traits yield higher Rao's Q."""
        q_close = raos_quadratic_entropy([[0], [1]], [1, 1])
        q_far = raos_quadratic_entropy([[0], [100]], [1, 1])
        assert q_far > q_close


# ============================================================================
# Community Weighted Mean (CWM) tests
# ============================================================================


class TestCommunityWeightedMean:
    """Test Community Weighted Mean trait values."""

    def test_cwm_equal_abundances(self) -> None:
        """Equal abundances yield the arithmetic mean of traits."""
        cwm = community_weighted_mean([[2, 4], [8, 6]], [1, 1])
        assert abs(cwm[0] - 5.0) < 1e-10
        assert abs(cwm[1] - 5.0) < 1e-10

    def test_cwm_weighted_by_abundance(self) -> None:
        """CWM weighted toward the dominant species."""
        cwm = community_weighted_mean([[0], [10]], [9, 1])
        # CWM = 0*0.9 + 10*0.1 = 1.0
        assert abs(cwm[0] - 1.0) < 1e-10

    def test_cwm_single_trait(self) -> None:
        """CWM with one trait returns a single-element list."""
        cwm = community_weighted_mean([[5], [15]], [1, 1])
        assert len(cwm) == 1
        assert abs(cwm[0] - 10.0) < 1e-10

    def test_cwm_multiple_traits(self) -> None:
        """CWM returns one value per trait."""
        cwm = community_weighted_mean(
            [[1, 2, 3], [4, 5, 6]],
            [1, 1],
        )
        assert len(cwm) == 3
        assert abs(cwm[0] - 2.5) < 1e-10
        assert abs(cwm[1] - 3.5) < 1e-10
        assert abs(cwm[2] - 4.5) < 1e-10

    def test_cwm_zero_abundance_excluded(self) -> None:
        """Species with zero abundance are excluded from CWM."""
        cwm = community_weighted_mean([[0], [100], [50]], [0, 1, 1])
        # Only species 1 (100) and 2 (50) -> mean = 75
        assert abs(cwm[0] - 75.0) < 1e-10


# ============================================================================
# Functional Redundancy tests
# ============================================================================


class TestFunctionalRedundancy:
    """Test Functional Redundancy (Simpson - Rao's Q)."""

    def test_redundancy_identical_traits(self) -> None:
        """Identical traits yield Simpson minus 0 = Simpson."""
        fr = functional_redundancy([[0], [0]], [5, 5])
        # Simpson = 1 - (0.5^2 + 0.5^2) = 0.5
        # Rao's Q = 0 (distance = 0)
        assert abs(fr - 0.5) < 1e-10

    def test_redundancy_diverse_traits(self) -> None:
        """With diverse traits, redundancy should be less than Simpson."""
        fr = functional_redundancy([[0], [10]], [5, 5])
        # Simpson = 0.5, Rao's Q = 5.0 -> FR = 0.5 - 5.0 = -4.5
        # (negative when distances are un-normalised)
        assert fr < 0.5

    def test_redundancy_single_species_zero(self) -> None:
        """Fewer than 2 present species gives redundancy = 0."""
        fr = functional_redundancy([[1], [2]], [0, 5])
        assert fr == 0.0


# ============================================================================
# Functional Beta Diversity tests
# ============================================================================


class TestFunctionalBetaDiversity:
    """Test Functional Beta Diversity decomposition."""

    def test_beta_div_output_keys(self) -> None:
        """Output has total, turnover, and nestedness."""
        result = functional_beta_diversity(
            [[0], [5]],
            [1, 1],
            [[10], [15]],
            [1, 1],
        )
        for key in ("total", "turnover", "nestedness"):
            assert key in result

    def test_beta_div_total_nonnegative(self) -> None:
        """Total functional beta diversity is non-negative."""
        result = functional_beta_diversity(
            [[0], [5]],
            [1, 1],
            [[10], [15]],
            [1, 1],
        )
        assert result["total"] >= -1e-10

    def test_beta_div_identical_communities_zero(self) -> None:
        """Identical communities have zero beta diversity."""
        result = functional_beta_diversity(
            [[0], [5], [10]],
            [1, 1, 1],
            [[0], [5], [10]],
            [1, 1, 1],
        )
        # With identical communities, total beta should be 0 or near-0
        assert result["total"] < 0.1

    def test_beta_div_turnover_plus_nestedness_leq_total(self) -> None:
        """Turnover + nestedness should equal total."""
        result = functional_beta_diversity(
            [[0], [5]],
            [1, 1],
            [[10], [15]],
            [1, 1],
        )
        assert abs(result["turnover"] + result["nestedness"] - result["total"]) < 1e-10


# ============================================================================
# Functional Diversity Suite (wrapper) tests
# ============================================================================


class TestFunctionalDiversitySuite:
    """Test the convenience wrapper that computes all metrics."""

    def test_suite_returns_all_keys(self) -> None:
        """Suite returns fric, feve, fdiv, fdis, raos_q, cwm, redundancy."""
        traits = [[0], [5], [10]]
        ab = [10, 10, 10]
        suite = functional_diversity_suite(traits, ab)
        expected_keys = {"fric", "feve", "fdiv", "fdis", "raos_q", "cwm", "redundancy"}
        assert set(suite.keys()) == expected_keys

    def test_suite_values_consistent_with_individual(self) -> None:
        """Suite values match individual function calls."""
        traits = [[0, 0], [1, 0], [0, 1], [1, 1]]
        ab = [10, 20, 30, 40]
        suite = functional_diversity_suite(traits, ab)

        fric_direct = functional_richness(traits, ab)
        feve_direct = functional_evenness(traits, ab)
        fdiv_direct = functional_divergence(traits, ab)
        fdis_direct = functional_dispersion(traits, ab)
        rq_direct = raos_quadratic_entropy(traits, ab)
        cwm_direct = community_weighted_mean(traits, ab)
        fr_direct = functional_redundancy(traits, ab)

        assert abs(suite["fric"] - fric_direct) < 1e-10
        assert abs(suite["feve"] - feve_direct) < 1e-10
        assert abs(suite["fdiv"] - fdiv_direct) < 1e-10
        assert abs(suite["fdis"] - fdis_direct) < 1e-10
        assert abs(suite["raos_q"] - rq_direct) < 1e-10
        assert suite["cwm"] == cwm_direct
        assert abs(suite["redundancy"] - fr_direct) < 1e-10

    def test_suite_with_many_species(self) -> None:
        """Suite handles larger communities without error."""
        n = 20
        # Use non-collinear 2D traits: vary both dimensions independently
        traits = [[i * 0.5, (i % 5) * 1.2 + (i // 5) * 0.7] for i in range(n)]
        ab = [max(1, (i + 1) * 3) for i in range(n)]
        suite = functional_diversity_suite(traits, ab)
        assert suite["fric"] > 0.0
        assert len(suite["cwm"]) == 2
