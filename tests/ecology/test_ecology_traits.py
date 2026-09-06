"""Tests for ecology traits subpackage (numpy-based trait metrics).

Covers community_weighted_mean, functional_richness, rao_quadratic_entropy,
and trait_diversity, including zero-sum abundance guards (regression:
all-zero abundances previously produced NaN instead of zeros).

Uses real implementations only (real-implementation policy).
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.ecology.traits.functional import (
    TraitData,
    community_weighted_mean,
    functional_richness,
    rao_quadratic_entropy,
    trait_diversity,
)


class TestCommunityWeightedMean:
    def test_cwm_weighted_by_abundance(self) -> None:
        abundances = np.array([3.0, 1.0])
        trait_values = np.array([[10.0], [0.0]])

        cwm = community_weighted_mean(abundances, trait_values)

        assert cwm["Trait_0"] == pytest.approx(7.5)  # (3*10 + 1*0) / 4

    def test_cwm_zero_sum_abundances_returns_zeros(self) -> None:
        """All-zero abundances yield zero CWM instead of NaN (regression)."""
        abundances = np.zeros(3)
        trait_values = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

        cwm = community_weighted_mean(abundances, trait_values)

        assert all(value == 0.0 for value in cwm.values())

    def test_cwm_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="match"):
            community_weighted_mean(np.array([1.0, 2.0]), np.array([[1.0]]))


class TestRaoQuadraticEntropy:
    def test_rao_two_species(self) -> None:
        # Full double-counted loop: Q = sum_ij d_ij p_i p_j = 2 * (10 * 0.25) = 5.0
        q = rao_quadratic_entropy(np.array([1.0, 1.0]), np.array([[0.0], [10.0]]))

        assert q == pytest.approx(5.0)

    def test_rao_zero_sum_abundances_is_zero(self) -> None:
        """All-zero abundances yield Q = 0 instead of NaN (regression)."""
        q = rao_quadratic_entropy(np.zeros(2), np.array([[0.0], [10.0]]))

        assert q == 0.0


class TestFunctionalRichness:
    def test_fric_1d_is_range(self) -> None:
        fric = functional_richness(np.array([[1.0], [5.0], [3.0]]))

        assert fric == pytest.approx(4.0)

    def test_fric_degenerate_config_is_zero(self) -> None:
        """Co-linear 2D points have zero hull area, not an error."""
        fric = functional_richness(np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]))

        assert fric == pytest.approx(0.0)

    def test_fric_2d_hull_area(self) -> None:
        # Unit square corners: hull area = 1.0
        points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        fric = functional_richness(points)

        assert fric == pytest.approx(1.0)


class TestTraitDiversity:
    def test_trait_diversity_metrics_populated(self) -> None:
        abundances = np.array([4.0, 3.0, 2.0, 1.0])
        trait_data = TraitData(
            species=["s1", "s2", "s3", "s4"],
            trait_names=["body", "length"],
            values=np.array(
                [
                    [10.0, 1.0],
                    [8.0, 2.0],
                    [6.0, 3.0],
                    [4.0, 4.0],
                ]
            ),
        )

        result = trait_diversity(abundances, trait_data)

        assert result.fric >= 0.0
        assert 0.0 <= result.feve <= 1.0
        assert result.fdiv >= 0.0
        assert result.rao_q >= 0.0
        assert result.cwm["body"] == pytest.approx(np.average([10, 8, 6, 4], weights=abundances))
