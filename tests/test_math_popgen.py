"""Tests for mathematical population genetics functionality."""

from __future__ import annotations

import pytest

from metainformant.math import (
    fixation_probability,
    hardy_weinberg_genotype_freqs,
    mutation_update,
    popgen,
    selection_update,
)


class TestMathPopgen:
    """Test enhanced mathematical population genetics functionality."""

    def test_effective_population_size_estimation(self):
        """Test effective population size estimation from heterozygosity."""
        # Test with reasonable heterozygosity values and mutation rate
        mu = 1e-8  # Typical mutation rate

        ne_50 = popgen.effective_population_size_from_heterozygosity(0.5, mutation_rate=mu)
        assert ne_50 > 0

        ne_90 = popgen.effective_population_size_from_heterozygosity(0.9, mutation_rate=mu)
        assert ne_90 > ne_50  # Higher heterozygosity should give higher Ne estimate

        ne_10 = popgen.effective_population_size_from_heterozygosity(0.1, mutation_rate=mu)
        assert ne_10 < ne_50  # Lower heterozygosity should give lower Ne estimate

        # Test edge cases - should raise ValueError
        with pytest.raises(ValueError):
            popgen.effective_population_size_from_heterozygosity(0.0, mutation_rate=mu)

        with pytest.raises(ValueError):
            popgen.effective_population_size_from_heterozygosity(1.0, mutation_rate=mu)

        # Test with invalid mutation rate
        with pytest.raises(ValueError):
            popgen.effective_population_size_from_heterozygosity(0.5, mutation_rate=0.0)

    def test_inbreeding_from_fst(self):
        """Test inbreeding coefficient estimation from F_ST."""
        # Test with reasonable F_ST values
        f_01 = popgen.inbreeding_coefficient_from_fst(0.1)
        f_05 = popgen.inbreeding_coefficient_from_fst(0.5)
        assert f_05 > f_01  # Higher F_ST should give higher inbreeding estimate

        # Test edge cases
        f_0 = popgen.inbreeding_coefficient_from_fst(0.0)
        assert f_0 == 0.0

        f_1 = popgen.inbreeding_coefficient_from_fst(1.0)
        assert f_1 == float("inf")

    def test_linkage_disequilibrium_decay(self):
        """Test linkage disequilibrium decay distance calculation."""
        # Test with reasonable values
        r2_50 = 0.5
        recomb_rate = 1e-8  # 1 cM/Mb
        distance_50 = popgen.linkage_disequilibrium_decay_distance(r2_50, recomb_rate)

        r2_10 = 0.1
        distance_10 = popgen.linkage_disequilibrium_decay_distance(r2_10, recomb_rate)

        # Higher r² means stronger LD (loci closer together), so smaller distance
        # Formula: distance = -ln(r²) / (2 * recombination_rate)
        # Since -ln is decreasing, higher r² gives smaller distance
        assert distance_50 < distance_10

        # Test edge cases
        distance_0 = popgen.linkage_disequilibrium_decay_distance(0.0, recomb_rate)
        assert distance_0 == 0.0

        distance_100 = popgen.linkage_disequilibrium_decay_distance(1.0, recomb_rate)
        assert distance_100 == float("inf")

        # Test with zero recombination rate
        distance_inf = popgen.linkage_disequilibrium_decay_distance(0.5, 0.0)
        assert distance_inf == float("inf")

    def test_coalescent_time_calculation(self):
        """Test coalescent time to most recent common ancestor."""
        # Test with reasonable values
        tmrca_10 = popgen.coalescent_time_to_mrca(10, 1000)
        tmrca_100 = popgen.coalescent_time_to_mrca(100, 1000)

        # Larger sample size should take longer to coalesce
        assert tmrca_100 > tmrca_10

        # Test edge cases
        tmrca_1 = popgen.coalescent_time_to_mrca(1, 1000)
        assert tmrca_1 == 0.0  # Single sample has no coalescence time

        tmrca_2 = popgen.coalescent_time_to_mrca(2, 1000)
        assert tmrca_2 > 0  # Two samples should have coalescence time

    def test_input_validation(self):
        """Test input validation for all new functions."""
        # Test invalid heterozygosity values
        with pytest.raises(ValueError):
            popgen.effective_population_size_from_heterozygosity(-0.1, mutation_rate=1e-8)

        with pytest.raises(ValueError):
            popgen.effective_population_size_from_heterozygosity(1.1, mutation_rate=1e-8)

        # Test invalid F_ST values
        with pytest.raises(ValueError):
            popgen.inbreeding_coefficient_from_fst(-0.1)

        with pytest.raises(ValueError):
            popgen.inbreeding_coefficient_from_fst(1.1)

        # Test invalid r² values
        assert popgen.linkage_disequilibrium_decay_distance(-0.1, 1e-8) == 0.0  # r² <= 0 returns 0.0

        assert popgen.linkage_disequilibrium_decay_distance(1.1, 1e-8) == float("inf")  # r² >= 1 returns inf

    def test_mathematical_consistency(self):
        """Test mathematical consistency of calculations."""
        # Test that different approaches give consistent results where applicable
        # For example, inbreeding coefficient should be consistent across methods

        fst = 0.1
        inbreeding1 = popgen.inbreeding_coefficient_from_fst(fst)

        # Test with different subpopulations (shouldn't affect result for 2-pop case)
        inbreeding2 = popgen.inbreeding_coefficient_from_fst(fst, subpopulations=5)
        assert inbreeding1 == inbreeding2  # Same result for different subpopulation counts in this method

    def test_realistic_biological_values(self):
        """Test with realistic biological parameter values."""
        # Test with typical population genetics values
        # Effective population size for humans: ~10,000
        # Human mutation rate ~1e-8 per base pair per generation
        ne_human = popgen.effective_population_size_from_heterozygosity(0.001, mutation_rate=1e-8)
        assert ne_human > 0  # Should be positive
        # For very low H, Ne estimate should be reasonable
        assert ne_human < 1e6  # Shouldn't be unreasonably large

        # F_ST values typically 0.01-0.1 for human populations
        fst_typical = 0.05
        inbreeding_typical = popgen.inbreeding_coefficient_from_fst(fst_typical)
        assert 0 < inbreeding_typical < 1  # Should be reasonable

        # LD decay distances for typical recombination rates
        # Human recombination rate ~1e-8 per base pair
        distance_typical = popgen.linkage_disequilibrium_decay_distance(0.5, 1e-8)
        assert distance_typical > 1000  # Should be thousands of base pairs

    def test_numerical_stability(self):
        """Test numerical stability of calculations."""
        # Test with very small values
        tiny_fst = 1e-10
        inbreeding_tiny = popgen.inbreeding_coefficient_from_fst(tiny_fst)
        assert inbreeding_tiny > 0

        # Test with very large values
        large_fst = 0.999999
        inbreeding_large = popgen.inbreeding_coefficient_from_fst(large_fst)
        assert inbreeding_large > 1000  # Should be large but finite

        # Test coalescent time with large population
        tmrca_large = popgen.coalescent_time_to_mrca(1000, 1e6)
        assert tmrca_large > 0

    def test_hardy_weinberg_genotype_freqs_basic(self):
        """Test Hardy-Weinberg genotype frequency calculations."""
        p2, two_pq, q2 = hardy_weinberg_genotype_freqs(0.25)
        assert abs(p2 - 0.0625) < 1e-9
        assert abs(two_pq - 0.375) < 1e-9
        assert abs(q2 - 0.5625) < 1e-9

    def test_selection_update_balancing_and_directional(self):
        """Test selection update calculations."""
        # Directional selection favoring A (w_AA > w_aa)
        p_next = selection_update(0.2, fitness_AA=1.2, fitness_Aa=1.1, fitness_aa=1.0)
        assert p_next > 0.2

        # No selection: mean fitness cancels out
        p_same = selection_update(0.2, 1.0, 1.0, 1.0)
        assert abs(p_same - 0.2) < 1e-12

    def test_mutation_update_forward_and_back(self):
        """Test mutation update calculations."""
        p_next = mutation_update(0.5, mu=0.01, nu=0.02)
        # Expected: 0.5*(0.99) + 0.5*0.02 = 0.495 + 0.01 = 0.505
        assert abs(p_next - 0.505) < 1e-12

    def test_fixation_probability_limits(self):
        """Test fixation probability calculations."""
        # Neutral case equals initial frequency
        assert abs(fixation_probability(0.1, 1000, 0.0) - 0.1) < 1e-12
        # Certain boundaries
        assert fixation_probability(0.0, 1000, 0.1) == 0.0
        assert fixation_probability(1.0, 1000, -0.1) == 1.0
