"""Tests for metagenomics diversity submodule.

Tests alpha/beta diversity, rarefaction, ordination.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

import math

import pytest

from metainformant.metagenomics.diversity.metrics import (
    alpha_diversity,
    beta_diversity,
    ordination,
    permanova,
    rarefaction_curve,
    rarefy,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def equal_community() -> list[int]:
    """Perfectly even community: 5 species, 20 individuals each."""
    return [20, 20, 20, 20, 20]


@pytest.fixture
def single_species() -> list[int]:
    """Community dominated by one species."""
    return [100, 0, 0, 0, 0]


@pytest.fixture
def diverse_community() -> list[int]:
    """Realistically diverse community."""
    return [50, 30, 10, 5, 3, 1, 1]


@pytest.fixture
def two_sample_matrix() -> list[list[float]]:
    """Two samples for beta diversity."""
    return [
        [10.0, 20.0, 30.0, 0.0],
        [10.0, 20.0, 30.0, 0.0],
    ]


@pytest.fixture
def three_sample_matrix() -> list[list[float]]:
    """Three different samples for beta diversity."""
    return [
        [50.0, 30.0, 10.0, 5.0, 5.0],
        [5.0, 5.0, 10.0, 30.0, 50.0],
        [20.0, 20.0, 20.0, 20.0, 20.0],
    ]


# ---------------------------------------------------------------------------
# Tests: alpha_diversity -- Shannon
# ---------------------------------------------------------------------------


class TestAlphaDiversityShannon:
    """Tests for Shannon diversity."""

    def test_equal_community_equals_ln_s(self, equal_community: list[int]) -> None:
        """For a perfectly even community, Shannon H = ln(S)."""
        result = alpha_diversity(equal_community, metric="shannon")
        expected = math.log(5)
        assert result["value"] == pytest.approx(expected, abs=1e-6)
        assert result["metric"] == "shannon"
        assert result["n_species"] == 5

    def test_single_species_is_zero(self, single_species: list[int]) -> None:
        """Single-species community has Shannon H = 0."""
        result = alpha_diversity(single_species, metric="shannon")
        assert result["value"] == pytest.approx(0.0, abs=1e-6)

    def test_positive_for_diverse(self, diverse_community: list[int]) -> None:
        result = alpha_diversity(diverse_community, metric="shannon")
        assert result["value"] > 0.0

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            alpha_diversity([], metric="shannon")


# ---------------------------------------------------------------------------
# Tests: alpha_diversity -- Simpson
# ---------------------------------------------------------------------------


class TestAlphaDiversitySimpson:
    """Tests for Simpson diversity (dominance form D = sum(p^2))."""

    def test_single_species_dominance_one(self, single_species: list[int]) -> None:
        """Single species: Simpson D = 1.0 (dominance = all probability in one species)."""
        result = alpha_diversity(single_species, metric="simpson")
        # Simpson index returns p^2 dominance. For single species p=1, D=1.
        assert result["value"] == pytest.approx(1.0, abs=1e-6)

    def test_equal_community_dominance(self, equal_community: list[int]) -> None:
        """Equal community: Simpson D = 1/S for equal proportions."""
        result = alpha_diversity(equal_community, metric="simpson")
        assert result["value"] == pytest.approx(1.0 / 5.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Tests: alpha_diversity -- Observed species
# ---------------------------------------------------------------------------


class TestAlphaDiversityObserved:
    """Tests for observed species count."""

    def test_count_nonzero(self, equal_community: list[int]) -> None:
        result = alpha_diversity(equal_community, metric="observed")
        assert result["value"] == 5.0

    def test_single_species_count(self, single_species: list[int]) -> None:
        result = alpha_diversity(single_species, metric="observed")
        assert result["value"] == 1.0


# ---------------------------------------------------------------------------
# Tests: alpha_diversity -- other metrics
# ---------------------------------------------------------------------------


class TestAlphaDiversityOther:
    """Tests for chao1, ace, invsimpson, fisher_alpha, pielou_evenness."""

    def test_chao1(self, diverse_community: list[int]) -> None:
        result = alpha_diversity(diverse_community, metric="chao1")
        # Chao1 >= observed
        assert result["value"] >= result["n_species"]

    def test_ace(self, diverse_community: list[int]) -> None:
        result = alpha_diversity(diverse_community, metric="ace")
        assert result["value"] > 0.0

    def test_invsimpson(self, equal_community: list[int]) -> None:
        result = alpha_diversity(equal_community, metric="invsimpson")
        # For equal proportions, 1/D = S
        assert result["value"] == pytest.approx(5.0, abs=1e-4)

    def test_fisher_alpha(self, diverse_community: list[int]) -> None:
        result = alpha_diversity(diverse_community, metric="fisher_alpha")
        assert result["value"] > 0.0

    def test_pielou_evenness(self, equal_community: list[int]) -> None:
        result = alpha_diversity(equal_community, metric="pielou_evenness")
        # Perfectly even = 1.0
        assert result["value"] == pytest.approx(1.0, abs=1e-6)

    def test_invalid_metric_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid metric"):
            alpha_diversity([10, 20], metric="bogus")


# ---------------------------------------------------------------------------
# Tests: beta_diversity
# ---------------------------------------------------------------------------


class TestBetaDiversity:
    """Tests for beta diversity distance matrices."""

    def test_bray_curtis_identical_is_zero(self, two_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(two_sample_matrix, metric="bray_curtis")
        assert result["distance_matrix"][0][1] == pytest.approx(0.0, abs=1e-6)
        assert result["metric"] == "bray_curtis"
        assert result["n_samples"] == 2

    def test_bray_curtis_range(self, three_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(three_sample_matrix, metric="bray_curtis")
        dm = result["distance_matrix"]
        for i in range(3):
            for j in range(3):
                assert 0.0 <= dm[i][j] <= 1.0

    def test_jaccard_identical_is_zero(self, two_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(two_sample_matrix, metric="jaccard")
        assert result["distance_matrix"][0][1] == pytest.approx(0.0, abs=1e-6)

    def test_jaccard_range(self, three_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(three_sample_matrix, metric="jaccard")
        dm = result["distance_matrix"]
        for i in range(3):
            for j in range(3):
                assert 0.0 <= dm[i][j] <= 1.0

    def test_aitchison(self, three_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(three_sample_matrix, metric="aitchison")
        dm = result["distance_matrix"]
        # Diagonal should be zero
        for i in range(3):
            assert dm[i][i] == pytest.approx(0.0, abs=1e-6)

    def test_symmetric(self, three_sample_matrix: list[list[float]]) -> None:
        result = beta_diversity(three_sample_matrix, metric="bray_curtis")
        dm = result["distance_matrix"]
        for i in range(3):
            for j in range(3):
                assert dm[i][j] == pytest.approx(dm[j][i], abs=1e-10)

    def test_too_few_samples_raises(self) -> None:
        with pytest.raises(ValueError, match="At least 2"):
            beta_diversity([[1.0, 2.0]], metric="bray_curtis")

    def test_invalid_metric_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid metric"):
            beta_diversity([[1.0], [2.0]], metric="fake")


# ---------------------------------------------------------------------------
# Tests: rarefaction_curve
# ---------------------------------------------------------------------------


class TestRarefactionCurve:
    """Tests for rarefaction curves."""

    def test_increasing_richness(self, diverse_community: list[int]) -> None:
        result = rarefaction_curve(diverse_community, n_iterations=10, seed=42)
        assert "depths" in result
        assert "mean_species" in result
        assert len(result["depths"]) == len(result["mean_species"])
        # Mean species should be generally non-decreasing
        species = result["mean_species"]
        assert species[-1] >= species[0]

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            rarefaction_curve([])

    def test_negative_values_raises(self) -> None:
        with pytest.raises(ValueError, match="negative"):
            rarefaction_curve([10, -5, 3])


# ---------------------------------------------------------------------------
# Tests: rarefy
# ---------------------------------------------------------------------------


class TestRarefy:
    """Tests for rarefaction subsampling."""

    def test_sum_equals_depth(self, diverse_community: list[int]) -> None:
        total = sum(diverse_community)
        depth = total // 2
        result = rarefy(diverse_community, depth, seed=42)
        assert sum(result) == depth
        assert len(result) == len(diverse_community)

    def test_full_depth_preserves_total(self, equal_community: list[int]) -> None:
        total = sum(equal_community)
        result = rarefy(equal_community, total, seed=42)
        assert sum(result) == total

    def test_zero_depth(self, equal_community: list[int]) -> None:
        result = rarefy(equal_community, 0)
        assert all(v == 0 for v in result)

    def test_exceeds_total_raises(self, equal_community: list[int]) -> None:
        total = sum(equal_community)
        with pytest.raises(ValueError, match="exceeds total"):
            rarefy(equal_community, total + 1)


# ---------------------------------------------------------------------------
# Tests: permanova
# ---------------------------------------------------------------------------


class TestPermanova:
    """Tests for PERMANOVA."""

    def test_basic_permanova(self) -> None:
        # Two groups with distinct distances
        dm = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.85],
            [0.8, 0.7, 0.0, 0.1],
            [0.9, 0.85, 0.1, 0.0],
        ]
        groups = ["A", "A", "B", "B"]
        result = permanova(dm, groups, n_permutations=99, seed=42)
        assert "pseudo_f" in result
        assert "p_value" in result
        assert "r_squared" in result
        assert 0.0 <= result["p_value"] <= 1.0
        assert result["r_squared"] >= 0.0

    def test_one_group_raises(self) -> None:
        dm = [[0.0, 0.5], [0.5, 0.0]]
        with pytest.raises(ValueError, match="At least 2 groups"):
            permanova(dm, ["A", "A"])


# ---------------------------------------------------------------------------
# Tests: ordination
# ---------------------------------------------------------------------------


class TestOrdination:
    """Tests for PCoA and NMDS ordination."""

    def test_pcoa_returns_coordinates(self) -> None:
        dm = [
            [0.0, 0.5, 0.8],
            [0.5, 0.0, 0.6],
            [0.8, 0.6, 0.0],
        ]
        result = ordination(dm, method="pcoa", n_components=2)
        assert "coordinates" in result
        assert len(result["coordinates"]) == 3
        assert all(len(row) == 2 for row in result["coordinates"])
        assert "eigenvalues" in result

    def test_nmds_returns_coordinates(self) -> None:
        dm = [
            [0.0, 0.5, 0.8],
            [0.5, 0.0, 0.6],
            [0.8, 0.6, 0.0],
        ]
        result = ordination(dm, method="nmds", n_components=2)
        assert "coordinates" in result
        assert len(result["coordinates"]) == 3

    def test_invalid_method_raises(self) -> None:
        dm = [[0.0, 0.5], [0.5, 0.0]]
        with pytest.raises(ValueError, match="Invalid method"):
            ordination(dm, method="fake")
