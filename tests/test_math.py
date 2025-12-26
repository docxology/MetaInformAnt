"""Tests for mathematical utilities."""

import pytest
from metainformant.math import popgen


class TestMath:
    """Test mathematical functionality."""

    def test_correlation_coefficient(self):
        """Test Pearson correlation coefficient calculation."""
        from metainformant.math import correlation_coefficient
        
        # Perfect positive correlation
        x = [1, 2, 3, 4, 5]
        y = [2, 4, 6, 8, 10]
        corr = correlation_coefficient(x, y)
        assert abs(corr - 1.0) < 1e-10
        
        # Perfect negative correlation
        y_neg = [10, 8, 6, 4, 2]
        corr_neg = correlation_coefficient(x, y_neg)
        assert abs(corr_neg + 1.0) < 1e-10
        
        # No correlation
        y_random = [1, 5, 2, 8, 3]
        corr_random = correlation_coefficient(x, y_random)
        assert abs(corr_random) < 0.5  # Should be close to zero

    def test_linear_regression(self):
        """Test linear regression calculation."""
        from metainformant.math import linear_regression
        
        # Simple linear relationship: y = 2x + 1
        x = [1, 2, 3, 4, 5]
        y = [3, 5, 7, 9, 11]
        
        slope, intercept, r_squared = linear_regression(x, y)
        
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept - 1.0) < 1e-10
        assert abs(r_squared - 1.0) < 1e-10

    def test_fisher_exact_test(self):
        """Test Fisher's exact test calculation."""
        from metainformant.math import fisher_exact_test
        
        # 2x2 contingency table
        a, b, c, d = 10, 5, 3, 12
        odds_ratio, p_value = fisher_exact_test(a, b, c, d)
        
        assert odds_ratio > 0
        assert 0 <= p_value <= 1

    def test_shannon_entropy(self):
        """Test Shannon entropy calculation."""
        from metainformant.math import shannon_entropy
        
        # Uniform distribution
        uniform = [0.25, 0.25, 0.25, 0.25]
        entropy_uniform = shannon_entropy(uniform)
        assert abs(entropy_uniform - 2.0) < 1e-10  # log2(4) = 2
        
        # Non-uniform distribution
        non_uniform = [0.5, 0.3, 0.2]
        entropy_non = shannon_entropy(non_uniform)
        assert entropy_non < entropy_uniform

    def test_jensen_shannon_divergence(self):
        """Test Jensen-Shannon divergence calculation."""
        from metainformant.math import jensen_shannon_divergence
        
        # Identical distributions
        p = [0.3, 0.7]
        q = [0.3, 0.7]
        jsd = jensen_shannon_divergence(p, q)
        assert jsd == 0.0
        
        # Different distributions
        p_diff = [0.5, 0.5]
        q_diff = [0.3, 0.7]
        jsd_diff = jensen_shannon_divergence(p_diff, q_diff)
        assert jsd_diff > 0

    def test_effective_population_size_estimation(self):
        """Test effective population size estimation."""
        # Test with reasonable heterozygosity values and mutation rate
        mu = 1e-8  # Typical mutation rate
        ne_50 = popgen.effective_population_size_from_heterozygosity(0.5, mutation_rate=mu)
        assert ne_50 > 0

        ne_90 = popgen.effective_population_size_from_heterozygosity(0.9, mutation_rate=mu)
        assert ne_90 > ne_50  # Higher heterozygosity should give higher Ne estimate

        ne_10 = popgen.effective_population_size_from_heterozygosity(0.1, mutation_rate=mu)
        assert ne_10 < ne_50  # Lower heterozygosity should give lower Ne estimate

    def test_inbreeding_from_fst(self):
        """Test inbreeding coefficient estimation from F_ST."""
        f_01 = popgen.inbreeding_coefficient_from_fst(0.1)
        f_05 = popgen.inbreeding_coefficient_from_fst(0.5)
        assert f_05 > f_01  # Higher F_ST should give higher inbreeding estimate

    def test_linkage_disequilibrium_decay(self):
        """Test linkage disequilibrium decay distance calculation."""
        r2_50 = 0.5
        recomb_rate = 1e-8
        distance_50 = popgen.linkage_disequilibrium_decay_distance(r2_50, recomb_rate)
    
        r2_10 = 0.1
        distance_10 = popgen.linkage_disequilibrium_decay_distance(r2_10, recomb_rate)
    
        # Higher r² means stronger LD, which typically means shorter distance
        # (loci are closer together, less recombination has occurred)
        # So distance_50 < distance_10
        assert distance_50 < distance_10

    def test_coalescent_time_calculation(self):
        """Test coalescent time to most recent common ancestor."""
        tmrca_10 = popgen.coalescent_time_to_mrca(10, 1000)
        tmrca_100 = popgen.coalescent_time_to_mrca(100, 1000)

        # Larger sample size should take longer to coalesce
        assert tmrca_100 > tmrca_10

    def test_input_validation(self):
        """Test input validation for all new functions."""
        from metainformant.math import correlation_coefficient, linear_regression, jensen_shannon_divergence
        
        # Test invalid correlation inputs
        with pytest.raises(ValueError):
            correlation_coefficient([1, 2], [1, 2, 3])  # Different lengths
            
        # Too short lists return 0.0, not ValueError
        assert correlation_coefficient([1], [2]) == 0.0

        # Test invalid regression inputs
        with pytest.raises(ValueError):
            linear_regression([1, 2], [1, 2, 3])  # Different lengths

        # Test invalid JSD inputs
        with pytest.raises(ValueError):
            jensen_shannon_divergence([0.5, 0.5], [0.3, 0.3, 0.4])  # Different lengths

    def test_mathematical_consistency(self):
        """Test mathematical consistency of calculations."""
        from metainformant.math import correlation_coefficient, linear_regression
        
        # Test that correlation and regression are consistent
        x = [1, 2, 3, 4, 5]
        y = [2, 4, 6, 8, 10]  # Perfect correlation
        
        corr = correlation_coefficient(x, y)
        slope, intercept, r_squared = linear_regression(x, y)
        
        # Should have perfect correlation
        assert abs(corr - 1.0) < 1e-10
        assert abs(r_squared - 1.0) < 1e-10
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept - 0.0) < 1e-10

    def test_realistic_biological_values(self):
        """Test with realistic biological parameter values."""
        # Test with typical population genetics values
        # Effective population size for humans: ~10,000
        # Human mutation rate ~1e-8 per base pair per generation
        ne_human = popgen.effective_population_size_from_heterozygosity(0.001, mutation_rate=1e-8)
        assert ne_human > 0  # Should be positive

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
        from metainformant.math import correlation_coefficient, linear_regression
        
        # Test with very small values
        x_small = [1e-10, 2e-10, 3e-10]
        y_small = [1e-10, 2e-10, 3e-10]
        corr_small = correlation_coefficient(x_small, y_small)
        assert abs(corr_small - 1.0) < 1e-6  # Should still detect perfect correlation
        
        # Test with very large values
        x_large = [1e10, 2e10, 3e10]
        y_large = [2e10, 4e10, 6e10]
        slope_large, intercept_large, r_squared_large = linear_regression(x_large, y_large)
        assert abs(slope_large - 2.0) < 1e-6
        assert abs(r_squared_large - 1.0) < 1e-6
