"""Population genetics statistical functions.

This module provides mathematical functions for population genetics statistics,
coalescent theory, and evolutionary computations.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def expected_pairwise_diversity(sample_size: int, theta: float) -> float:
    """Calculate expected pairwise nucleotide diversity (π).

    Args:
        sample_size: Number of sequences (n)
        theta: Population mutation parameter (4Neμ)

    Returns:
        Expected pairwise diversity
    """
    if sample_size < 2:
        return 0.0

    # E[π] = θ * (n/(n-1)) * sum(1/i for i in 1 to n-1)
    harmonic_sum = sum(1.0 / i for i in range(1, sample_size))
    return theta * (sample_size / (sample_size - 1)) * harmonic_sum


def expected_segregating_sites(sample_size: int, theta: float) -> float:
    """Calculate expected number of segregating sites.

    Args:
        sample_size: Number of sequences (n)
        theta: Population mutation parameter (4Neμ)

    Returns:
        Expected number of segregating sites
    """
    if sample_size < 2:
        return 0.0

    # E[S] = θ * sum(1/i for i in 1 to n-1)
    harmonic_sum = sum(1.0 / i for i in range(1, sample_size))
    return theta * harmonic_sum


def expected_coalescent_waiting_times(sample_size: int, Ne: float) -> List[float]:
    """Calculate expected waiting times for coalescent events.

    Args:
        sample_size: Starting number of lineages
        Ne: Effective population size

    Returns:
        List of expected waiting times for each coalescent interval
    """
    waiting_times = []
    current_lineages = sample_size

    while current_lineages > 1:
        # Expected waiting time for next coalescence
        # Time is in units of 4Ne generations
        rate = current_lineages * (current_lineages - 1) / 2
        expected_time = 1.0 / rate
        waiting_times.append(expected_time)
        current_lineages -= 1

    return waiting_times


def expected_r2_from_Ne_c(recombination_rate: float, Ne: float, distance_bp: float) -> float:
    """Calculate expected LD (r²) from effective population size and recombination.

    Args:
        recombination_rate: Recombination rate per bp per generation
        Ne: Effective population size
        distance_bp: Distance between sites in base pairs

    Returns:
        Expected r² value
    """
    # r² decays as 1/(1 + 4Ne*c*d) where c is recombination rate, d is distance
    c = recombination_rate
    d = distance_bp

    denominator = 1 + 4 * Ne * c * d
    return 1.0 / denominator


def equilibrium_heterozygosity_infinite_alleles(theta: float) -> float:
    """Calculate equilibrium heterozygosity under infinite alleles model.

    Args:
        theta: Population mutation parameter (4Neμ)

    Returns:
        Equilibrium heterozygosity (H)
    """
    # H = θ / (1 + θ) for infinite alleles model
    return theta / (1 + theta)


def fixation_probability(
    initial_frequency: float | None = None,
    population_size: int | None = None,
    *,
    selection_coefficient: float = 0.0,
    Ne: int | None = None,
    p0: float | None = None,
) -> float:
    """Calculate fixation probability under selection.

    Supports multiple calling conventions:
    - fixation_probability(p0, Ne, selection_coefficient=s)
    - fixation_probability(selection_coefficient=s, population_size=N)

    Args:
        initial_frequency: Initial allele frequency (p0)
        population_size: Population size (N)
        selection_coefficient: Selection coefficient (s)
        Ne: Alias for population_size
        p0: Alias for initial_frequency

    Returns:
        Probability of fixation
    """
    # Handle parameter aliases
    if Ne is not None and population_size is None:
        population_size = Ne
    if p0 is not None and initial_frequency is None:
        initial_frequency = p0

    # Default values
    if initial_frequency is None:
        initial_frequency = 1.0 / (2 * population_size) if population_size else 0.5
    if population_size is None:
        raise ValueError("population_size (or Ne) must be provided")

    s = selection_coefficient
    N = population_size
    p = initial_frequency

    # For neutral allele, fixation probability equals initial frequency
    if s == 0:
        return p

    # Kimura's diffusion approximation
    # P(fix) = (1 - exp(-4*N*s*p)) / (1 - exp(-4*N*s))
    try:
        exp_term = 4 * N * s
        numerator = 1 - math.exp(-exp_term * p)
        denominator = 1 - math.exp(-exp_term)

        if abs(denominator) < 1e-10:
            return p  # Very weak selection, nearly neutral

        return numerator / denominator
    except OverflowError:
        # For very strong selection
        if s > 0:
            return 1.0 if p > 0 else 0.0
        else:
            return 0.0


def bottleneck_effective_size(
    initial_size: int, bottleneck_size: int, bottleneck_duration: int, final_size: int
) -> float:
    """Calculate effective population size after bottleneck.

    Args:
        initial_size: Initial population size
        bottleneck_size: Population size during bottleneck
        bottleneck_duration: Duration of bottleneck in generations
        final_size: Final population size

    Returns:
        Effective population size
    """
    # Simplified approximation
    # Ne = 1 / (1/N_initial + 1/N_bottleneck * duration + 1/N_final)
    Ne = 1.0 / (1.0 / initial_size + bottleneck_duration * 1.0 / bottleneck_size + 1.0 / final_size)
    return Ne


def effective_size_from_family_size_variance(family_sizes: List[int]) -> float:
    """Estimate effective population size from variance in family sizes.

    Args:
        family_sizes: List of offspring counts per parent

    Returns:
        Estimated effective population size
    """
    if not family_sizes:
        return 0.0

    n = len(family_sizes)
    mean_family_size = sum(family_sizes) / n
    variance = sum((x - mean_family_size) ** 2 for x in family_sizes) / n

    # Ne = k / (V_k - 1) where k is mean family size, V_k is variance
    if variance > 1:
        return mean_family_size / (variance - 1)
    else:
        return float("inf")  # No variance means infinite Ne


def bootstrap_confidence_interval(
    values: List[float], confidence_level: float = 0.95, n_bootstraps: int = 1000
) -> Tuple[float, float]:
    """Calculate bootstrap confidence interval.

    Args:
        values: Original sample values
        confidence_level: Confidence level (0.95 for 95% CI)
        n_bootstraps: Number of bootstrap resamples

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if not values:
        return (0.0, 0.0)

    np.random.seed(42)  # For reproducibility
    bootstrap_means = []

    for _ in range(n_bootstraps):
        # Resample with replacement
        sample = np.random.choice(values, size=len(values), replace=True)
        bootstrap_means.append(np.mean(sample))

    # Calculate confidence interval
    lower_percentile = (1 - confidence_level) / 2 * 100
    upper_percentile = (1 + confidence_level) / 2 * 100

    lower_bound = np.percentile(bootstrap_means, lower_percentile)
    upper_bound = np.percentile(bootstrap_means, upper_percentile)

    return (lower_bound, upper_bound)


def calculate_confidence_intervals(
    values: List[float], confidence_level: float = 0.95, method: str = "bootstrap"
) -> tuple[float, float]:
    """Calculate confidence intervals for population genetic statistics.

    Args:
        values: List of observed values
        confidence_level: Confidence level (default 0.95)
        method: Method for CI calculation ("bootstrap", "normal", "percentile")

    Returns:
        Tuple of (lower_bound, upper_bound)

    Examples:
        >>> values = [0.1, 0.15, 0.12, 0.18, 0.09]
        >>> lower, upper = calculate_confidence_intervals(values)
        >>> print(f"95% CI: [{lower:.3f}, {upper:.3f}]")
    """
    if not values:
        raise ValueError("Values list cannot be empty")

    if not (0 < confidence_level < 1):
        raise ValueError("Confidence level must be between 0 and 1")

    values_arr = np.array(values)
    n = len(values_arr)

    if method == "bootstrap":
        # Bootstrap confidence interval
        n_bootstrap = min(1000, n * 10)  # Reasonable bootstrap sample size
        bootstrap_means = []

        np.random.seed(42)  # For reproducibility
        for _ in range(n_bootstrap):
            sample = np.random.choice(values_arr, size=n, replace=True)
            bootstrap_means.append(np.mean(sample))

        lower_percentile = (1 - confidence_level) / 2 * 100
        upper_percentile = (1 + confidence_level) / 2 * 100

        lower_bound = np.percentile(bootstrap_means, lower_percentile)
        upper_bound = np.percentile(bootstrap_means, upper_percentile)

    elif method == "normal":
        # Normal approximation confidence interval
        mean_val = np.mean(values_arr)
        se = np.std(values_arr, ddof=1) / np.sqrt(n)

        # z-score for confidence level
        if confidence_level == 0.95:
            z = 1.96
        elif confidence_level == 0.99:
            z = 2.576
        else:
            # General case using normal distribution
            from scipy import stats

            z = stats.norm.ppf((1 + confidence_level) / 2)

        lower_bound = mean_val - z * se
        upper_bound = mean_val + z * se

    elif method == "percentile":
        # Percentile method
        lower_percentile = (1 - confidence_level) / 2 * 100
        upper_percentile = (1 + confidence_level) / 2 * 100

        lower_bound = np.percentile(values_arr, lower_percentile)
        upper_bound = np.percentile(values_arr, upper_percentile)

    else:
        raise ValueError(f"Unknown method: {method}. Use 'bootstrap', 'normal', or 'percentile'")

    return float(lower_bound), float(upper_bound)


def compare_population_statistic(
    pop1_stats: Dict[str, Any], pop2_stats: Dict[str, Any], statistic_name: str
) -> Dict[str, Any]:
    """Compare a population genetic statistic between two populations.

    Args:
        pop1_stats: Statistics for population 1
        pop2_stats: Statistics for population 2
        statistic_name: Name of the statistic to compare

    Returns:
        Dictionary with comparison results
    """
    if statistic_name not in pop1_stats or statistic_name not in pop2_stats:
        raise ValueError(f"Statistic '{statistic_name}' not found in both populations")

    stat1 = pop1_stats[statistic_name]
    stat2 = pop2_stats[statistic_name]

    # Simple t-test approximation (for demonstration)
    # In practice, you'd want proper statistical tests
    diff = stat1 - stat2
    pooled_var = (pop1_stats.get("variance", 1.0) + pop2_stats.get("variance", 1.0)) / 2
    se = math.sqrt(pooled_var * 2 / 10)  # Approximate SE

    if se > 0:
        t_stat = abs(diff) / se
        # Approximate p-value using normal distribution
        p_value = 2 * (1 - 0.5 * (1 + math.erf(t_stat / math.sqrt(2))))
    else:
        t_stat = 0.0
        p_value = 1.0

    return {
        "statistic_name": statistic_name,
        "population1_value": stat1,
        "population2_value": stat2,
        "difference": diff,
        "test_statistic": t_stat,
        "p_value": p_value,
        "significant": p_value < 0.05,
    }


def compare_statistics(stats1: Dict[str, float], stats2: Dict[str, float]) -> Dict[str, Any]:
    """Compare multiple population statistics between two populations.

    Args:
        stats1: Statistics dictionary for population 1
        stats2: Statistics dictionary for population 2

    Returns:
        Dictionary with comparison results for each statistic
    """
    results = {}
    common_stats = set(stats1.keys()) & set(stats2.keys())

    for stat_name in common_stats:
        results[stat_name] = compare_population_statistic(
            {"variance": 1.0, **stats1}, {"variance": 1.0, **stats2}, stat_name  # Add default variance
        )

    return results


def detect_outliers(values: List[float], method: str = "iqr", threshold: float = 1.5) -> Dict[str, Any]:
    """Detect outliers in population genetic data.

    Args:
        values: List of values to check for outliers
        method: Outlier detection method ("iqr", "zscore")
        threshold: Threshold for outlier detection

    Returns:
        Dictionary with outlier detection results
    """
    if not values:
        return {"outliers": [], "outlier_indices": []}

    values_array = np.array(values)

    if method == "iqr":
        q1 = np.percentile(values_array, 25)
        q3 = np.percentile(values_array, 75)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr

        outlier_mask = (values_array < lower_bound) | (values_array > upper_bound)

    elif method == "zscore":
        z_scores = np.abs((values_array - np.mean(values_array)) / np.std(values_array))
        outlier_mask = z_scores > threshold

    else:
        raise ValueError(f"Unknown outlier detection method: {method}")

    outlier_indices = np.where(outlier_mask)[0].tolist()
    outliers = values_array[outlier_mask].tolist()

    return {
        "outliers": outliers,
        "outlier_indices": outlier_indices,
        "method": method,
        "threshold": threshold,
        "total_values": len(values),
        "outlier_percentage": (len(outliers) / len(values)) * 100,
    }


def permutation_test(values1: List[float], values2: List[float], n_permutations: int = 1000) -> Dict[str, Any]:
    """Perform permutation test to compare two distributions.

    Args:
        values1: Values from first distribution
        values2: Values from second distribution
        n_permutations: Number of permutations to perform

    Returns:
        Dictionary with permutation test results
    """
    if not values1 or not values2:
        return {"p_value": 1.0, "test_statistic": 0.0}

    # Observed difference in means
    observed_diff = abs(np.mean(values1) - np.mean(values2))

    # Combine data
    combined = np.array(values1 + values2)
    n1 = len(values1)

    # Count permutations with difference >= observed
    count = 0
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_diff = abs(np.mean(combined[:n1]) - np.mean(combined[n1:]))
        if perm_diff >= observed_diff:
            count += 1

    p_value = count / n_permutations

    return {
        "test_statistic": observed_diff,
        "p_value": p_value,
        "n_permutations": n_permutations,
        "significant": p_value < 0.05,
    }


def tajimas_d_outliers(d_values: List[float], window_size: int = 10) -> Dict[str, Any]:
    """Detect outlier Tajima's D values using sliding window analysis.

    Args:
        d_values: List of Tajima's D values along a chromosome/sequence
        window_size: Size of sliding window for outlier detection

    Returns:
        Dictionary with outlier detection results
    """
    if len(d_values) < window_size:
        return {"outliers": [], "outlier_indices": []}

    outliers = []
    outlier_indices = []

    # Calculate local statistics using sliding window
    for i in range(len(d_values) - window_size + 1):
        window = d_values[i : i + window_size]
        window_mean = np.mean(window)
        window_std = np.std(window)

        if window_std > 0:
            # Check each value in window for being outlier
            for j, value in enumerate(window):
                z_score = abs(value - window_mean) / window_std
                if z_score > 3.0:  # 3 standard deviations
                    global_idx = i + j
                    if global_idx not in outlier_indices:
                        outliers.append(value)
                        outlier_indices.append(global_idx)

    return {
        "outliers": outliers,
        "outlier_indices": outlier_indices,
        "window_size": window_size,
        "total_windows": len(d_values) - window_size + 1,
        "outlier_percentage": (len(outliers) / len(d_values)) * 100 if d_values else 0,
    }


def variance(values: List[float]) -> float:
    """Calculate variance of a list of values."""
    if not values:
        return 0.0
    return float(np.var(values, ddof=1)) if len(values) > 1 else 0.0


def standard_deviation(values: List[float]) -> float:
    """Calculate standard deviation of a list of values."""
    if not values:
        return 0.0
    return float(np.std(values, ddof=1)) if len(values) > 1 else 0.0


def skewness(values: List[float]) -> float:
    """Calculate skewness of a list of values."""
    if not values or len(values) < 3:
        return 0.0
    from scipy.stats import skew

    return float(skew(values))


def kurtosis(values: List[float]) -> float:
    """Calculate kurtosis of a list of values."""
    if not values or len(values) < 4:
        return 0.0
    from scipy.stats import kurtosis as scipy_kurtosis

    return float(scipy_kurtosis(values))
