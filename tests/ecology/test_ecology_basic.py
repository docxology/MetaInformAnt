"""Basic tests for ecology module functionality."""

import pytest

from metainformant.ecology.analysis import community


class TestCommunityEcology:
    """Test community ecology analysis functions."""

    def test_shannon_diversity_empty(self):
        """Test Shannon diversity with empty abundances."""
        result = community.shannon_diversity([])
        assert result == 0.0

    def test_shannon_diversity_single_species(self):
        """Test Shannon diversity with single species."""
        result = community.shannon_diversity([1])
        assert result == 0.0

    def test_shannon_diversity_equal_abundances(self):
        """Test Shannon diversity with equal abundances."""
        import math

        result = community.shannon_diversity([1, 1, 1, 1])
        expected = -4 * (0.25 * math.log(0.25))  # -sum(p * ln(p)) = -4 * (0.25 * ln(0.25))
        assert abs(result - expected) < 1e-10

    def test_shannon_diversity_unequal_abundances(self):
        """Test Shannon diversity with unequal abundances."""
        result = community.shannon_diversity([1, 2, 3, 4])
        # Should be between 0 and max possible
        assert 0 <= result <= 2.0

    def test_shannon_diversity_zero_abundances(self):
        """Test Shannon diversity with zero total abundance."""
        result = community.shannon_diversity([0, 0, 0])
        assert result == 0.0

    def test_simpson_diversity_basic(self):
        """Test Simpson diversity index calculation."""
        result = community.simpson_diversity([1, 1, 1, 1])
        expected = 1 - sum((0.25) ** 2 for _ in range(4))  # 1 - 4*(0.25)^2
        assert abs(result - expected) < 1e-10

    def test_simpson_diversity_empty(self):
        """Test Simpson diversity with empty abundances."""
        result = community.simpson_diversity([])
        assert result == 0.0

    def test_species_richness_basic(self):
        """Test species richness calculation."""
        abundances = [1, 2, 0, 3, 0, 1]
        result = community.species_richness(abundances)
        assert result == 4  # Only 4 non-zero abundances

    def test_species_richness_all_zeros(self):
        """Test species richness with all zero abundances."""
        result = community.species_richness([0, 0, 0])
        assert result == 0

    def test_species_richness_empty(self):
        """Test species richness with empty list."""
        result = community.species_richness([])
        assert result == 0

    def test_pielou_evenness_perfect_evenness(self):
        """Test Pielou evenness with perfect evenness."""
        abundances = [1, 1, 1, 1]  # Perfect evenness
        result = community.pielou_evenness(abundances)

        # Should be 1.0 for perfect evenness
        assert abs(result - 1.0) < 1e-10

    def test_pielou_evenness_no_evenness(self):
        """Test Pielou evenness with no evenness."""
        abundances = [1, 0, 0, 0]  # One species dominates
        result = community.pielou_evenness(abundances)

        # When S=1, pielou_evenness returns 1.0 (edge case: single species)
        assert abs(result - 1.0) < 1e-10

    def test_pielou_evenness_zero_richness(self):
        """Test Pielou evenness with zero richness."""
        result = community.pielou_evenness([])
        assert result == 0.0

    def test_chao1_estimator_basic(self):
        """Test Chao1 species richness estimator."""
        abundances = [10, 8, 6, 4, 2, 1, 1, 1]
        result = community.chao1_estimator(abundances)

        # Chao1 should be >= observed richness
        observed_richness = community.species_richness(abundances)
        assert result >= observed_richness

    def test_chao1_estimator_all_singletons(self):
        """Test Chao1 with all singletons."""
        abundances = [1, 1, 1, 1, 1]
        result = community.chao1_estimator(abundances)
        # Chao1 formula: S_obs + (f1^2)/(2*f2) where f1=5 singletons, f2=0 doubletons
        # When f2=0: S_obs + (f1*(f1-1))/2 = 5 + (5*4)/2 = 5 + 10 = 15.0
        assert result == 15.0

    def test_chao1_estimator_empty(self):
        """Test Chao1 with empty abundances."""
        result = community.chao1_estimator([])
        assert result == 0.0

    def test_community_metrics_basic(self):
        """Test comprehensive community metrics."""
        abundances = [10, 8, 6, 4, 2, 1, 1, 1]

        metrics = community.community_metrics(abundances)

        # Check all metrics are present
        expected_keys = ["shannon", "simpson", "richness", "pielou", "chao1"]
        for key in expected_keys:
            assert key in metrics
            assert isinstance(metrics[key], (int, float))

        # Check values are reasonable
        assert metrics["richness"] == 8
        assert 0 <= metrics["shannon"] <= 3.0
        assert 0 <= metrics["simpson"] <= 1.0
        assert 0 <= metrics["pielou"] <= 1.0
        assert metrics["chao1"] >= metrics["richness"]

    def test_community_metrics_empty(self):
        """Test community metrics with empty data."""
        metrics = community.community_metrics([])
        assert metrics["richness"] == 0
        assert metrics["shannon"] == 0.0
        assert metrics["simpson"] == 0.0
        assert metrics["pielou"] == 0.0
        assert metrics["chao1"] == 0.0

    def test_community_metrics_all_zeros(self):
        """Test community metrics with all zero abundances."""
        metrics = community.community_metrics([0, 0, 0])
        assert metrics["richness"] == 0
        assert metrics["shannon"] == 0.0
        assert metrics["simpson"] == 0.0
        assert metrics["pielou"] == 0.0
        assert metrics["chao1"] == 0.0
