"""Tests for population genetics statistical methods."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.math.population_genetics.statistics import (
    bottleneck_effective_size,
    bootstrap_confidence_interval,
    calculate_confidence_intervals,
    compare_population_statistic,
    compare_statistics,
    deterministic_replicate_seeds,
    detect_outliers,
    equilibrium_heterozygosity_infinite_alleles,
    expected_pairwise_diversity,
    expected_segregating_sites,
    fixation_probability,
    kurtosis,
    normal_approximate_power,
    permutation_test,
    sample_size_for_power,
    skewness,
    standard_deviation,
    tajimas_d_outliers,
    variance,
)


def test_bootstrap_confidence_interval():
    """Test bootstrap confidence interval calculation."""
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    result = bootstrap_confidence_interval(data, np.mean, n_bootstrap=100, random_state=42)

    assert "statistic" in result
    assert "ci_lower" in result
    assert "ci_upper" in result
    assert "confidence_level" in result
    assert result["confidence_level"] == 0.95
    assert result["ci_lower"] <= result["statistic"] <= result["ci_upper"]


def test_permutation_test():
    """Test permutation test."""
    group1 = [1.0, 2.0, 3.0]
    group2 = [4.0, 5.0, 6.0]

    result = permutation_test(group1, group2, n_permutations=1000, random_state=42)

    assert "statistic" in result
    assert "p_value" in result
    assert 0.0 <= result["p_value"] <= 1.0
    assert result["p_value"] <= 0.1  # Should be significant (allowing boundary)


def test_detect_outliers():
    """Test outlier detection."""
    values = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0]

    result = detect_outliers(values, method="zscore", threshold=2.0)

    assert "outlier_indices" in result
    assert "outlier_values" in result
    assert len(result["outlier_indices"]) > 0
    assert 100.0 in result["outlier_values"]


def test_tajimas_d_outliers():
    """Test Tajima's D outlier detection."""
    d_values = [-0.5, -0.3, 0.1, 5.0, 0.2]

    result = tajimas_d_outliers(d_values, threshold=1.0)  # Lower threshold to detect 5.0

    assert "outlier_indices" in result
    assert len(result["outlier_indices"]) > 0


def test_compare_statistics():
    """Test statistical comparison."""
    stat1 = [1.0, 2.0, 3.0]
    stat2 = [4.0, 5.0, 6.0]

    result = compare_statistics(stat1, stat2, test_type="mannwhitney")

    assert "test_statistic" in result
    assert "p_value" in result
    assert "test_type" in result
    assert result["test_type"] == "mannwhitney"


def test_test_population_difference():
    """Test population difference testing."""
    pop1_stats = {"pi": 0.01, "theta": 0.01}
    pop2_stats = {"pi": 0.02, "theta": 0.02}

    result = compare_population_statistic(pop1_stats, pop2_stats, "pi")

    assert "test_statistic" in result
    assert "p_value" in result
    assert "statistic_name" in result
    assert result["statistic_name"] == "pi"


def test_calculate_confidence_intervals():
    """Test confidence interval calculation."""
    statistics = {"pi": 0.01, "theta": 0.01}

    result = calculate_confidence_intervals(statistics, method="normal")

    assert "pi" in result
    assert "theta" in result
    assert "ci_lower" in result["pi"]
    assert "ci_upper" in result["pi"]


def test_bootstrap_with_empty_data():
    """Test bootstrap with empty data."""
    data = []
    result = bootstrap_confidence_interval(data, np.mean, n_bootstrap=100)

    assert "statistic" in result
    assert result["n_bootstrap"] == 0
    # Empty data should return NaN for statistic
    assert np.isnan(result["statistic"])


def test_permutation_test_with_empty_groups():
    """Test permutation test with empty groups."""
    group1 = []
    group2 = [1.0, 2.0]

    result = permutation_test(group1, group2, n_permutations=100)

    assert result["p_value"] == 1.0
    assert result["n_permutations"] == 0


def test_expected_pairwise_diversity_and_segregating_sites_degenerate_n():
    assert expected_pairwise_diversity(1, 5.0) == 0.0
    assert expected_segregating_sites(0, 5.0) == 0.0
    # n=2: E[pi] = theta * 2/1 * 1 = 2*theta; E[S] = theta
    assert expected_pairwise_diversity(2, 0.5) == pytest.approx(1.0)
    assert expected_segregating_sites(2, 0.5) == pytest.approx(0.5)


def test_equilibrium_heterozygosity_infinite_alleles():
    assert equilibrium_heterozygosity_infinite_alleles(1.0) == pytest.approx(0.5)
    assert equilibrium_heterozygosity_infinite_alleles(0.0) == pytest.approx(0.0)
    # mutation_rate form: theta = 4*Ne*mu
    assert equilibrium_heterozygosity_infinite_alleles(1000.0, mutation_rate=1e-6) == pytest.approx(
        (4 * 1000.0 * 1e-6) / (1 + 4 * 1000.0 * 1e-6)
    )


def test_bottleneck_effective_size_known_value():
    # Ne = 1 / (1/Ni + d/Nb + 1/Nf)
    expected = 1.0 / (1.0 / 10000.0 + 10.0 / 100.0 + 1.0 / 10000.0)
    assert bottleneck_effective_size(10000, 100, 10, 10000) == pytest.approx(expected)


def test_variance_standard_deviation_skewness_kurtosis_edges():
    assert variance([]) == 0.0
    assert variance([3.0]) == 0.0
    assert standard_deviation([]) == 0.0
    assert standard_deviation([3.0]) == 0.0
    assert variance([1.0, 2.0, 3.0]) == pytest.approx(2.0 / 3.0)
    assert standard_deviation([1.0, 2.0, 3.0]) == pytest.approx(np.sqrt(2.0 / 3.0))
    # Too few points for higher moments -> 0.0
    assert skewness([1.0, 2.0]) == 0.0
    assert kurtosis([1.0, 2.0, 3.0]) == 0.0
    assert skewness([1.0, 2.0, 3.0, 4.0]) == pytest.approx(0.0, abs=1e-12)
    assert kurtosis([1.0, 2.0, 3.0, 4.0]) == pytest.approx(-1.36)


def test_detect_outliers_zscore_method():
    result = detect_outliers([1.0, 1.1, 0.9, 1.05, 10.0], method="zscore", threshold=1.9)
    assert result["outlier_indices"] == [4]
    assert result["outlier_values"] == [10.0]


def test_detect_outliers_unknown_method_raises():
    with pytest.raises(ValueError, match="Unknown outlier detection method"):
        detect_outliers([1.0, 2.0], method="mad")


def test_compare_population_statistic_missing_statistic_raises():
    with pytest.raises(ValueError, match="not found"):
        compare_population_statistic({"pi": 0.1}, {"theta": 0.2}, "pi")


def test_compare_population_statistic_significance_fields():
    result = compare_population_statistic({"pi": 0.5, "variance": 0.01}, {"pi": 0.1, "variance": 0.01}, "pi")
    assert result["difference"] == pytest.approx(0.4)
    assert result["significant"] is True
    assert 0.0 <= result["p_value"] <= 1.0


def test_compare_statistics_vector_input():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    b = [2.0, 3.0, 4.0, 5.0, 6.0]
    result = compare_statistics(a, b, test_type="ttest")
    assert result["test_type"] == "ttest"
    assert 0.0 <= result["p_value"] <= 1.0


def test_calculate_confidence_intervals_normal_and_unknown_method():
    lower, upper = calculate_confidence_intervals([1.0, 2.0, 3.0, 4.0], method="normal")
    assert lower < 2.5 < upper
    with pytest.raises(ValueError, match="Unknown method"):
        calculate_confidence_intervals([1.0, 2.0], method="jackknife")


def test_calculate_confidence_intervals_empty_raises():
    with pytest.raises(ValueError, match="empty"):
        calculate_confidence_intervals([])


def test_normal_approximate_power_validation_and_value():
    # Huge effect size with n=100 per group -> power ~ 1
    assert normal_approximate_power(1.0, 100) > 0.99
    # Zero effect -> power equals half alpha under the one-tail formula
    assert normal_approximate_power(0.0, 100) == pytest.approx(0.025, abs=1e-4)
    with pytest.raises(ValueError, match="non-negative"):
        normal_approximate_power(-0.5, 50)
    with pytest.raises(ValueError, match=">= 2"):
        normal_approximate_power(0.5, 1)
    with pytest.raises(ValueError, match="alpha"):
        normal_approximate_power(0.5, 50, alpha=1.5)


def test_sample_size_for_power_inverse_of_power():
    n = sample_size_for_power(0.5, target_power=0.8)
    assert n >= 2
    assert normal_approximate_power(0.5, n) >= 0.8
    assert normal_approximate_power(0.5, n - 1) < 0.8
    with pytest.raises(ValueError, match="> 0"):
        sample_size_for_power(0.0)
    with pytest.raises(ValueError, match="target_power"):
        sample_size_for_power(0.5, target_power=0.01)


def test_deterministic_replicate_seeds_order_independent():
    seeds = deterministic_replicate_seeds(42, 5)
    assert len(seeds) == 5
    assert len(set(seeds)) == 5
    assert deterministic_replicate_seeds(42, 5) == seeds


def test_fixation_probability_aliases_and_limits():
    # Keyword alias conventions
    assert fixation_probability(p0=0.1, Ne=1000, selection_coefficient=0.0) == pytest.approx(0.1)
    assert fixation_probability(selection_coefficient=0.0, population_size=100) == pytest.approx(0.005)
    # Missing population size raises
    with pytest.raises(ValueError, match="population_size"):
        fixation_probability(initial_frequency=0.5)


def test_bootstrap_confidence_interval_default_statistic():
    result = bootstrap_confidence_interval([1.0, 2.0, 3.0], n_bootstraps=50, random_state=0)
    assert result["statistic"] == pytest.approx(2.0)
    assert result["n_bootstrap"] == 50
