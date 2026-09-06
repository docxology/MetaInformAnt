"""Tests for math selection module.

All tests follow real-implementation policy and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.math.population_genetics.selection import (
    kin_selection_response,
    multilevel_selection_decomposition,
    mutation_selection_balance_dominant,
    mutation_selection_balance_recessive,
    mutation_update,
    relative_fitness,
    selection_differential,
    selection_gradient,
    selection_intensity,
    selection_update,
)


class TestKinSelectionResponse:
    """Test Hamilton's rule calculation."""

    def test_kin_selection_favored(self):
        """Test when kin selection favors trait (r*b > c)."""
        # High relatedness (0.5), high benefit (0.4), low cost (0.1)
        result = kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.1)
        assert result == 0.1  # 0.5 * 0.4 - 0.1 = 0.2 - 0.1 = 0.1
        assert result > 0  # Selection favors

    def test_kin_selection_not_favored(self):
        """Test when kin selection does not favor trait (r*b < c)."""
        # Low relatedness (0.25), moderate benefit (0.3), high cost (0.2)
        result = kin_selection_response(relatedness=0.25, benefit=0.3, cost=0.2)
        assert result == -0.125  # 0.25 * 0.3 - 0.2 = 0.075 - 0.2 = -0.125
        assert result < 0  # Selection against

    def test_kin_selection_neutral(self):
        """Test when kin selection is neutral (r*b = c)."""
        # Relatedness 0.5, benefit 0.4, cost 0.2
        result = kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.2)
        assert result == 0.0  # 0.5 * 0.4 - 0.2 = 0.2 - 0.2 = 0.0

    def test_kin_selection_zero_relatedness(self):
        """Test with zero relatedness."""
        result = kin_selection_response(relatedness=0.0, benefit=1.0, cost=0.1)
        assert result == -0.1  # Always negative (cost only)

    def test_kin_selection_zero_cost(self):
        """Test with zero cost."""
        result = kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.0)
        assert result == 0.2  # 0.5 * 0.4 = 0.2


class TestMultilevelSelectionDecomposition:
    """Test multilevel selection decomposition."""

    def test_multilevel_selection_basic(self):
        """Test basic multilevel selection decomposition."""
        group_means = [1.0, 2.0, 3.0]
        individual_deviations = [0.1, -0.1, 0.2, -0.2, 0.0, 0.1]
        group_selection = 0.5
        individual_selection = 0.3

        between, within, total = multilevel_selection_decomposition(
            group_means=group_means,
            individual_deviations=individual_deviations,
            selection_strength_group=group_selection,
            selection_strength_individual=individual_selection,
        )

        assert isinstance(between, float)
        assert isinstance(within, float)
        assert isinstance(total, float)
        # Total should be sum of between and within components
        assert abs(total - (between + within)) < 1e-10

    def test_multilevel_selection_single_group(self):
        """Test with single group."""
        group_means = [1.0]
        individual_deviations = [0.1, -0.1, 0.0]
        group_selection = 0.5
        individual_selection = 0.3

        between, within, total = multilevel_selection_decomposition(
            group_means=group_means,
            individual_deviations=individual_deviations,
            selection_strength_group=group_selection,
            selection_strength_individual=individual_selection,
        )

        # With single group, between-group component should be zero
        assert between == 0.0
        assert within >= 0
        assert total == within

    def test_multilevel_selection_equal_groups(self):
        """Test with equal group means."""
        group_means = [2.0, 2.0, 2.0]
        individual_deviations = [0.1, -0.1, 0.2, -0.2, 0.0, 0.1]
        group_selection = 0.5
        individual_selection = 0.3

        between, within, total = multilevel_selection_decomposition(
            group_means=group_means,
            individual_deviations=individual_deviations,
            selection_strength_group=group_selection,
            selection_strength_individual=individual_selection,
        )

        # With equal group means, between-group component should be zero
        assert between == 0.0
        assert total == within

    def test_multilevel_selection_zero_individual_selection(self):
        """Test with zero individual-level selection."""
        group_means = [1.0, 2.0, 3.0]
        individual_deviations = [0.1, -0.1, 0.2, -0.2, 0.0, 0.1]
        group_selection = 0.5
        individual_selection = 0.0

        between, within, total = multilevel_selection_decomposition(
            group_means=group_means,
            individual_deviations=individual_deviations,
            selection_strength_group=group_selection,
            selection_strength_individual=individual_selection,
        )

        # With zero individual selection, within-group component should be zero
        assert within == 0.0
        assert total == between


class TestMutationUpdate:
    def test_known_update(self):
        # delta_p = -u*p + v*(1-p) = -0.1*0.9 + 0 = -0.09
        assert mutation_update(0.9, 0.1, 0.0) == pytest.approx(0.81)

    def test_clamps_to_unit_interval(self):
        assert mutation_update(1.0, 0.5, 0.0) == 0.5
        assert mutation_update(0.0, 0.0, 0.5) == 0.5

    def test_validation(self):
        with pytest.raises(ValueError, match="between 0 and 1"):
            mutation_update(1.5, 0.01, 0.01)
        with pytest.raises(ValueError, match="required"):
            mutation_update(0.5, mutation_rate_forward=0.01)


class TestSelectionUpdate:
    def test_additive_selection_known_value(self):
        # delta_p = s*p*q/2 = 0.5*0.2*0.8/2 = 0.04
        assert selection_update(0.2, selection_coefficient=0.5, dominance_coefficient=0.5) == pytest.approx(0.24)

    def test_fitness_mode_known_value(self):
        # wbar = 0.04*1.2 + 0.32*1.1 + 0.64*1.0 = 1.04
        p_new = selection_update(0.2, fitness_AA=1.2, fitness_Aa=1.1, fitness_aa=1.0)
        assert p_new == pytest.approx((0.04 * 1.2 + 0.5 * 0.32 * 1.1) / 1.04)

    def test_zero_mean_fitness_returns_original_frequency(self):
        assert selection_update(0.3, fitness_AA=0.0, fitness_Aa=0.0, fitness_aa=0.0) == 0.3

    def test_validation(self):
        with pytest.raises(ValueError, match="between 0 and 1"):
            selection_update(-0.1, 0.5)


class TestQuantitativeSelectionMetrics:
    def test_selection_differential_known_value(self):
        # Cov(trait, fitness)/mean(fitness) = (20/3)/2 = 10/3
        assert selection_differential([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]) == pytest.approx(10.0 / 3.0)

    def test_selection_gradient_known_value(self):
        # beta = Cov(trait, relative fitness)/Var(trait) = (10/3)/(200/3) = 0.05
        assert selection_gradient([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]) == pytest.approx(0.05)

    def test_selection_intensity_known_value(self):
        # i = S / sigma(trait) = (10/3)/sqrt(200/3)
        expected = (10.0 / 3.0) / (200.0 / 3.0) ** 0.5
        assert selection_intensity([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]) == pytest.approx(expected)

    def test_zero_trait_variance_is_handled(self):
        assert selection_gradient([1.0, 2.0, 3.0], [5.0, 5.0, 5.0]) == 0.0
        assert selection_intensity([1.0, 2.0, 3.0], [5.0, 5.0, 5.0]) == 0.0

    def test_length_mismatch_and_empty_raise(self):
        with pytest.raises(ValueError, match="same length"):
            selection_differential([1.0, 2.0], [1.0])
        with pytest.raises(ValueError, match="empty"):
            selection_gradient([], [])

    def test_relative_fitness_zero_mean(self):
        assert relative_fitness([0.0, 0.0]) == [0.0, 0.0]
        assert relative_fitness([]) == []


class TestMutationSelectionBalance:
    def test_recessive_known_value(self):
        assert mutation_selection_balance_recessive(1e-5, 0.1) == pytest.approx(0.01)

    def test_dominant_known_value(self):
        assert mutation_selection_balance_dominant(1e-5, 0.1) == pytest.approx(1e-4)

    def test_validation(self):
        with pytest.raises(ValueError, match="positive"):
            mutation_selection_balance_recessive(1e-5, 0.0)
        with pytest.raises(ValueError, match="negative"):
            mutation_selection_balance_recessive(-1e-5, 0.1)
        with pytest.raises(ValueError, match="positive"):
            mutation_selection_balance_dominant(1e-5, 0.0)
        with pytest.raises(ValueError, match="negative"):
            mutation_selection_balance_dominant(-1e-5, 0.1)
