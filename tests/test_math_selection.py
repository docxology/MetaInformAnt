"""Tests for math selection module.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.math.selection import (
    kin_selection_response,
    multilevel_selection_decomposition,
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
            selection_strength_individual=individual_selection
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
            selection_strength_individual=individual_selection
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
            selection_strength_individual=individual_selection
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
            selection_strength_individual=individual_selection
        )
        
        # With zero individual selection, within-group component should be zero
        assert within == 0.0
        assert total == between

