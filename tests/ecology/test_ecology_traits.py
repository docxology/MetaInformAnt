"""Regression tests for the consolidated ecology functional diversity API.

The numpy-based duplicate ``metainformant.ecology.traits.functional`` was
consolidated into ``metainformant.ecology.analysis.functional`` (the
Villeger et al. 2008 MST-based FEve/FDiv implementation). These tests keep
the zero-sum abundance regression guards and pin CWM, Rao's Q and FRic on
the consolidated API. Real implementations only.
"""

from __future__ import annotations

import pytest

from metainformant.ecology.analysis.functional import (
    community_weighted_mean,
    functional_diversity_suite,
    functional_richness,
    raos_quadratic_entropy,
)


class TestCommunityWeightedMean:
    def test_cwm_weighted_by_abundance(self) -> None:
        cwm = community_weighted_mean([[10.0], [0.0]], [3.0, 1.0])

        assert cwm == pytest.approx([7.5])  # (3*10 + 1*0) / 4

    def test_cwm_zero_sum_abundances_returns_zeros(self) -> None:
        """All-zero abundances yield zero CWM instead of NaN (regression)."""
        cwm = community_weighted_mean(
            [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], [0.0, 0.0, 0.0]
        )

        assert cwm == [0.0, 0.0]

    def test_cwm_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="abundances length"):
            community_weighted_mean([[1.0], [2.0]], [1.0])


class TestRaoQuadraticEntropy:
    def test_rao_two_species(self) -> None:
        # Double-counted pairs: Q = 2 * d * p1 * p2 = 2 * 10 * 0.25 = 5.0
        q = raos_quadratic_entropy([[0.0], [10.0]], [1.0, 1.0])

        assert q == pytest.approx(5.0)

    def test_rao_zero_sum_abundances_is_zero(self) -> None:
        """All-zero abundances yield Q = 0 instead of NaN (regression)."""
        q = raos_quadratic_entropy([[0.0], [10.0]], [0.0, 0.0])

        assert q == 0.0


class TestFunctionalRichness:
    def test_fric_1d_is_range(self) -> None:
        fric = functional_richness([[1.0], [5.0], [3.0]])

        assert fric == pytest.approx(4.0)

    def test_fric_degenerate_config_is_zero(self) -> None:
        """Co-linear 2D points have zero hull area, not an error."""
        fric = functional_richness([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])

        assert fric == pytest.approx(0.0)

    def test_fric_2d_hull_area(self) -> None:
        # Unit square corners: hull area = 1.0
        fric = functional_richness([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])

        assert fric == pytest.approx(1.0)


class TestFunctionalDiversitySuite:
    def test_suite_metrics_populated(self) -> None:
        """functional_diversity_suite replaces the removed traits.trait_diversity."""
        suite = functional_diversity_suite(
            [[10.0, 1.0], [8.0, 2.0], [6.0, 3.0], [4.0, 4.0]],
            [4.0, 3.0, 2.0, 1.0],
        )

        assert suite["cwm"] == pytest.approx([8.0, 2.0])
        assert suite["fric"] >= 0.0
        assert 0.0 <= suite["feve"] <= 1.0
        assert 0.0 <= suite["fdiv"] <= 1.0
        assert suite["raos_q"] >= 0.0
