"""Tests for advanced population genetics visualization functions."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.dna.population_viz import (
    plot_allele_frequency_spectrum,
    plot_bootstrap_distribution,
    plot_fst_matrix,
    plot_hardy_weinberg_test,
    plot_heterozygosity_distribution,
    plot_neutrality_test_suite,
    plot_outlier_detection,
    plot_pi_vs_theta,
    plot_permutation_test,
    plot_statistic_correlation_matrix,
    plot_statistic_distribution,
)


def test_plot_allele_frequency_spectrum():
    """Test allele frequency spectrum plot."""
    sfs = [10, 5, 3, 2, 1]
    fig = plot_allele_frequency_spectrum(sfs)
    
    assert fig is not None


def test_plot_heterozygosity_distribution():
    """Test heterozygosity distribution plot."""
    het_values = [0.1, 0.2, 0.3, 0.4, 0.5]
    fig = plot_heterozygosity_distribution(het_values)
    
    assert fig is not None


def test_plot_statistic_distribution():
    """Test statistic distribution plot."""
    stats = {
        "scenario1": [1.0, 2.0, 3.0],
        "scenario2": [4.0, 5.0, 6.0],
    }
    fig = plot_statistic_distribution(stats, plot_type="histogram")
    
    assert fig is not None


def test_plot_pi_vs_theta():
    """Test π vs θ plot."""
    pi_values = [0.01, 0.02, 0.03]
    theta_values = [0.01, 0.02, 0.03]
    fig = plot_pi_vs_theta(pi_values, theta_values)
    
    assert fig is not None


def test_plot_statistic_correlation_matrix():
    """Test statistic correlation matrix plot."""
    stats = {
        "pi": [0.01, 0.02, 0.03],
        "theta": [0.01, 0.02, 0.03],
        "d": [-0.5, 0.0, 0.5],
    }
    fig = plot_statistic_correlation_matrix(stats)
    
    assert fig is not None


def test_plot_fst_matrix():
    """Test Fst matrix plot."""
    fst_matrix = [[0.0, 0.1, 0.2], [0.1, 0.0, 0.15], [0.2, 0.15, 0.0]]
    fig = plot_fst_matrix(fst_matrix)
    
    assert fig is not None


def test_plot_neutrality_test_suite():
    """Test neutrality test suite plot."""
    test_results = {
        "tajimas_d": {"statistic": -0.5},
        "fu_and_li_d": {"statistic": -0.3},
        "fay_wu_h": {"statistic": -0.2},
    }
    fig = plot_neutrality_test_suite(test_results)
    
    assert fig is not None


def test_plot_hardy_weinberg_test():
    """Test Hardy-Weinberg test plot."""
    hwe_results = {
        "chi_square": 2.5,
        "p_value": 0.1,
        "degrees_of_freedom": 1,
        "hwe_deviated": False,
    }
    fig = plot_hardy_weinberg_test(hwe_results)
    
    assert fig is not None


def test_plot_bootstrap_distribution():
    """Test bootstrap distribution plot."""
    bootstrap_values = [1.0, 1.5, 2.0, 2.5, 3.0]
    fig = plot_bootstrap_distribution(bootstrap_values, observed_value=2.0)
    
    assert fig is not None


def test_plot_permutation_test():
    """Test permutation test plot."""
    permuted_values = [0.1, 0.2, 0.3, 0.4, 0.5]
    observed_value = 0.8
    fig = plot_permutation_test(permuted_values, observed_value, p_value=0.01)
    
    assert fig is not None


def test_plot_outlier_detection():
    """Test outlier detection plot."""
    statistic_values = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0]
    outlier_indices = [3]
    fig = plot_outlier_detection(statistic_values, outlier_indices)
    
    assert fig is not None


