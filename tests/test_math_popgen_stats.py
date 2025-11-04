"""Tests for population genetics statistical methods."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.math.popgen_stats import (
    bootstrap_confidence_interval,
    calculate_confidence_intervals,
    compare_statistics,
    detect_outliers,
    permutation_test,
    tajimas_d_outliers,
    test_population_difference,
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
    assert result["p_value"] < 0.1  # Should be significant


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
    
    result = tajimas_d_outliers(d_values, threshold=2.0)
    
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
    
    result = test_population_difference(pop1_stats, pop2_stats, "pi")
    
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


def test_permutation_test_with_empty_groups():
    """Test permutation test with empty groups."""
    group1 = []
    group2 = [1.0, 2.0]
    
    result = permutation_test(group1, group2, n_permutations=100)
    
    assert result["p_value"] == 1.0
    assert result["n_permutations"] == 0

