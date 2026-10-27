"""Population genetics statistical functions.

This module provides mathematical functions for population genetics statistics,
coalescent theory, and evolutionary computations.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Sequence

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


def expected_r2_from_Ne_c(recombination_rate: float, Ne: float | None = None, distance_bp: float = 1.0) -> float:
    """Calculate expected LD (r²) from effective population size and recombination.

    Args:
        recombination_rate: Recombination rate per bp per generation
        Ne: Effective population size
        distance_bp: Distance between sites in base pairs

    Returns:
        Expected r² value
    """
    # r² decays as 1/(1 + 4Ne*c*d) where c is recombination rate, d is distance
    if Ne is None:
        raise ValueError("Ne must be provided")
    if distance_bp == 1.0 and recombination_rate > 1 and Ne <= 0.5:
        Ne, recombination_rate = recombination_rate, Ne
    c = recombination_rate
    d = distance_bp

    denominator = 1 + 4 * Ne * c * d
    return 1.0 / denominator


def equilibrium_heterozygosity_infinite_alleles(theta: float, mutation_rate: float | None = None) -> float:
    """Calculate equilibrium heterozygosity under infinite alleles model.

    Args:
        theta: Population mutation parameter (4Neμ)

    Returns:
        Equilibrium heterozygosity (H)
    """
    if mutation_rate is not None:
        theta = 4.0 * theta * mutation_rate
    # H = θ / (1 + θ) for infinite alleles model
    return theta / (1 + theta)


def fixation_probability(
    initial_frequency: float | None = None,
    population_size: int | None = None,
    selection_coefficient: float = 0.0,
    *,
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
    """Calculate effective population size across a bottleneck.

    Ne is the harmonic mean of the census sizes (Wright 1938) over one generation
    at the initial size, ``bottleneck_duration`` generations at the bottleneck
    size and one generation at the final size:
    Ne = (d + 2) / (1/N_initial + d/N_bottleneck + 1/N_final).

    Args:
        initial_size: Initial population size
        bottleneck_size: Population size during bottleneck
        bottleneck_duration: Duration of bottleneck in generations (d)
        final_size: Final population size

    Returns:
        Effective population size

    Raises:
        ValueError: If a size is not positive or the duration is negative
    """
    if min(initial_size, bottleneck_size, final_size) <= 0:
        raise ValueError("Population sizes must be positive")
    if bottleneck_duration < 0:
        raise ValueError("Bottleneck duration cannot be negative")
    generations = bottleneck_duration + 2
    return generations / (1.0 / initial_size + bottleneck_duration / bottleneck_size + 1.0 / final_size)


def effective_size_from_family_size_variance(family_sizes: List[int] | float, Vk: float | None = None) -> float:
    """Estimate effective population size from variance in family sizes.

    Two forms:

    - With ``Vk``, ``family_sizes`` is the census size N (an int or a float) and the
      Crow-Denniston approximation Ne = (4N - 2) / (Vk + 2) is returned.
    - Without ``Vk``, ``family_sizes`` is a list (or tuple or array) of offspring counts, one
      per parent.

    Args:
        family_sizes: Census population size N when ``Vk`` is given, otherwise the offspring
            counts per parent
        Vk: Variance in family size, given together with the census size N

    Returns:
        Estimated effective population size

    Raises:
        TypeError: If ``Vk`` is given but ``family_sizes`` is not a number, or ``Vk`` is omitted
            but ``family_sizes`` is not a list, tuple or array
        ValueError: If N is not positive and finite, or Vk is negative or not finite
    """
    if Vk is not None:
        if isinstance(family_sizes, bool) or not isinstance(family_sizes, (int, float, np.integer, np.floating)):
            raise TypeError(
                f"family_sizes must be the census size N (an int or a float) when Vk is given, "
                f"got {type(family_sizes).__name__}"
            )
        N = float(family_sizes)
        if not math.isfinite(N) or N <= 0:
            raise ValueError(f"The census size N must be positive and finite, got {family_sizes}")
        if not math.isfinite(Vk) or Vk < 0:
            raise ValueError(f"Vk must be non-negative and finite, got {Vk}")
        return (4.0 * N - 2.0) / (Vk + 2.0)

    if not isinstance(family_sizes, (list, tuple, np.ndarray)):
        raise TypeError(
            f"family_sizes must be a list of offspring counts when Vk is not given, got {type(family_sizes).__name__}"
        )
    sizes = list(family_sizes)

    if not sizes:
        return 0.0

    n = len(sizes)
    mean_family_size = sum(sizes) / n
    variance = sum((x - mean_family_size) ** 2 for x in sizes) / n

    # Ne = k / (V_k - 1) where k is mean family size, V_k is variance
    if variance > 1:
        return mean_family_size / (variance - 1)
    else:
        return float("inf")  # No variance means infinite Ne


def bootstrap_confidence_interval(
    values: List[float],
    statistic_func: Any | None = None,
    confidence_level: float = 0.95,
    n_bootstrap: int | None = None,
    n_bootstraps: int = 1000,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Calculate bootstrap confidence interval.

    Args:
        values: Original sample values
        confidence_level: Confidence level (0.95 for 95% CI)
        n_bootstraps: Number of bootstrap resamples

    Returns:
        Dictionary with keys ``statistic``, ``ci_lower``, ``ci_upper``,
        ``confidence_level`` and ``n_bootstrap`.
    """
    if n_bootstrap is not None:
        n_bootstraps = n_bootstrap
    if statistic_func is None:
        statistic_func = np.mean
    if not values:
        return {
            "statistic": float("nan"),
            "ci_lower": float("nan"),
            "ci_upper": float("nan"),
            "confidence_level": confidence_level,
            "n_bootstrap": 0,
        }

    rng = np.random.default_rng(random_state)
    bootstrap_means = []

    for _ in range(n_bootstraps):
        sample = rng.choice(values, size=len(values), replace=True)
        bootstrap_means.append(float(statistic_func(sample)))

    # Calculate confidence interval
    lower_percentile = (1 - confidence_level) / 2 * 100
    upper_percentile = (1 + confidence_level) / 2 * 100

    lower_bound = float(np.percentile(bootstrap_means, lower_percentile))
    upper_bound = float(np.percentile(bootstrap_means, upper_percentile))

    return {
        "statistic": float(statistic_func(np.asarray(values))),
        "ci_lower": lower_bound,
        "ci_upper": upper_bound,
        "confidence_level": confidence_level,
        "n_bootstrap": n_bootstraps,
    }


def calculate_confidence_intervals(
    values: Sequence[float] | np.ndarray | Dict[str, Any], confidence_level: float = 0.95, method: str = "bootstrap"
) -> tuple[float, float] | Dict[str, Dict[str, float]]:
    """Calculate confidence intervals for population genetic statistics.

    Args:
        values: Observed values (a list or array), giving one interval; or a dictionary that maps
            each statistic name to its own observations (a list or array), giving one interval
            per name. An interval is estimated from the observations: a dictionary of single point
            values has nothing to estimate it from and raises ``ValueError``
        confidence_level: Confidence level (default 0.95)
        method: Method for CI calculation ("bootstrap", "normal", "percentile")

    Returns:
        Tuple of (lower_bound, upper_bound) for a list or array, or
        ``{name: {"ci_lower": ..., "ci_upper": ...}}`` for a dictionary of observations

    Raises:
        ValueError: If there are no values, ``confidence_level`` is outside (0, 1), ``method`` is
            unknown, or a dictionary entry is a single number rather than observations

    Examples:
        >>> values = [0.1, 0.15, 0.12, 0.18, 0.09]
        >>> lower, upper = calculate_confidence_intervals(values)
        >>> print(f"95% CI: [{lower:.3f}, {upper:.3f}]")
    """
    if isinstance(values, dict):
        if not values:
            raise ValueError("Values dictionary cannot be empty")
        intervals: Dict[str, Dict[str, float]] = {}
        for name, observations in values.items():
            if np.ndim(observations) == 0:
                raise ValueError(
                    f"Cannot build a confidence interval for '{name}' from the single value {observations!r}: "
                    "a point estimate carries no sampling information. Pass the observations behind it "
                    "(a list or array per statistic) instead."
                )
            lower, upper = _confidence_interval_of_sample(observations, confidence_level, method)
            intervals[name] = {"ci_lower": lower, "ci_upper": upper}
        return intervals

    return _confidence_interval_of_sample(values, confidence_level, method)


def _confidence_interval_of_sample(values: Any, confidence_level: float, method: str) -> tuple[float, float]:
    """Interval for one sample of observations; see :func:`calculate_confidence_intervals`."""
    values_arr = np.asarray(values, dtype=float)
    if values_arr.ndim != 1:
        raise ValueError("Values must be a one-dimensional list or array")
    if values_arr.size == 0:
        raise ValueError("Values list cannot be empty")

    if not (0 < confidence_level < 1):
        raise ValueError("Confidence level must be between 0 and 1")

    n = len(values_arr)

    if method == "bootstrap":
        # Bootstrap confidence interval
        n_bootstrap = min(1000, n * 10)  # Reasonable bootstrap sample size
        bootstrap_means = []

        # Fixed seed for reproducibility. A private RandomState draws the same stream as seeding
        # the global generator with 42 did, without resetting the caller's global numpy state.
        bootstrap_rng = np.random.RandomState(42)
        for _ in range(n_bootstrap):
            sample = bootstrap_rng.choice(values_arr, size=n, replace=True)
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

    # INFERENTIAL (post-evidence-freeze use): approximate t-test with normal-approx
    # p-value. Not for cross-species release outputs while the evidence manifest is
    # unfrozen (descriptive-only boundary, SHARED_CONTEXT Round 3+).
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


COMPARISON_TESTS = ("ttest", "mannwhitney", "kruskal")


def compare_statistics(stats1: Any, stats2: Any, test_type: str = "ttest") -> Dict[str, Any]:
    """Compare two groups of observations, or two dictionaries of named statistics.

    With two sequences of observations the requested test is run on them:

    - ``"ttest"`` (default): Welch's two-sample t-test for unequal variances
      (``scipy.stats.ttest_ind(equal_var=False)``)
    - ``"mannwhitney"``: two-sided Mann-Whitney U test (``scipy.stats.mannwhitneyu``)
    - ``"kruskal"``: Kruskal-Wallis H test of the two groups (``scipy.stats.kruskal``)

    Args:
        stats1: Observations for population 1, or a dictionary of statistics for it
        stats2: Observations for population 2, or a dictionary of statistics for it
        test_type: ``"ttest"``, ``"mannwhitney"`` or ``"kruskal"``

    Returns:
        For sequences, a dictionary with ``test_statistic``, ``p_value`` and the ``test_type``
        that was run. For two dictionaries, a dictionary with the comparison results for each
        statistic they share (see :func:`compare_population_statistic`); that comparison is an
        approximate normal-theory one for point values, so only ``test_type="ttest"`` applies
        to it.

    Raises:
        ValueError: If ``test_type`` is not one of the three above; if dictionaries are given with
            a rank test, which needs the observations themselves; or if the observations cannot
            be tested: not one-dimensional, not finite, too few (two per group for the t-test,
            one for the rank tests), or all identical so that the test is undefined

    INFERENTIAL (post-evidence-freeze use): t-test/Mann-Whitney/Kruskal-Wallis p-values are not
    for cross-species release outputs while the evidence manifest is unfrozen
    (descriptive-only boundary, SHARED_CONTEXT Round 3+).
    """
    if test_type not in COMPARISON_TESTS:
        raise ValueError(f"Unknown test_type {test_type!r}. Use one of: {', '.join(repr(t) for t in COMPARISON_TESTS)}")

    if not isinstance(stats1, dict) or not isinstance(stats2, dict):
        return _compare_observations(stats1, stats2, test_type)

    if test_type != "ttest":
        raise ValueError(
            f"test_type {test_type!r} ranks observations, but dictionaries hold only summary statistics; "
            "pass the observations as two sequences, or use test_type='ttest' for the approximate comparison"
        )

    results = {}
    common_stats = set(stats1.keys()) & set(stats2.keys())

    for stat_name in common_stats:
        results[stat_name] = compare_population_statistic(
            {"variance": 1.0, **stats1}, {"variance": 1.0, **stats2}, stat_name  # Add default variance
        )

    return results


def _compare_observations(stats1: Any, stats2: Any, test_type: str) -> Dict[str, Any]:
    """Run the ``test_type`` test on two samples; see :func:`compare_statistics`.

    A degenerate sample raises instead of yielding a made-up p-value.
    """
    from scipy import stats as scipy_stats

    a = np.asarray(stats1, dtype=float)
    b = np.asarray(stats2, dtype=float)
    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("The two groups must be one-dimensional sequences of observations")
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        raise ValueError("The two groups must hold only finite observations")
    minimum = 2 if test_type == "ttest" else 1
    if a.size < minimum or b.size < minimum:
        raise ValueError(
            f"The {test_type!r} comparison needs at least {minimum} observation(s) in each group, "
            f"got {a.size} and {b.size}"
        )

    result: Any  # a scipy result object (ttest_ind, mannwhitneyu or kruskal) with .statistic and .pvalue
    if test_type == "ttest":
        if np.ptp(a) == 0 and np.ptp(b) == 0:
            raise ValueError("Welch's t-test is undefined when every observation in both groups is constant")
        result = scipy_stats.ttest_ind(a, b, equal_var=False)
    else:
        if np.ptp(np.concatenate([a, b])) == 0:
            raise ValueError(f"The {test_type!r} comparison is undefined when all observations are identical")
        if test_type == "mannwhitney":
            result = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
        else:
            result = scipy_stats.kruskal(a, b)

    statistic, p_value = float(result.statistic), float(result.pvalue)
    if math.isnan(statistic) or math.isnan(p_value):
        raise ValueError(f"The {test_type!r} comparison is undefined for these observations")
    return {"test_statistic": statistic, "p_value": p_value, "test_type": test_type}


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
        "outlier_values": outliers,
        "method": method,
        "threshold": threshold,
        "total_values": len(values),
        "outlier_percentage": (len(outliers) / len(values)) * 100,
    }


def permutation_test(
    values1: List[float], values2: List[float], n_permutations: int = 1000, random_state: int | None = None
) -> Dict[str, Any]:
    """Perform permutation test to compare two distributions.

    Args:
        values1: Values from first distribution
        values2: Values from second distribution
        n_permutations: Number of permutations to perform

    Returns:
        Dictionary with permutation test results

    INFERENTIAL (post-evidence-freeze use): permutation p-values are not for
    cross-species release outputs while the evidence manifest is unfrozen
    (descriptive-only boundary, SHARED_CONTEXT Round 3+).
    """
    if not values1 or not values2:
        return {"p_value": 1.0, "test_statistic": 0.0, "statistic": 0.0, "n_permutations": 0}

    # Observed difference in means
    observed_diff = abs(np.mean(values1) - np.mean(values2))

    # Combine data
    combined = np.array(values1 + values2)
    n1 = len(values1)

    rng = np.random.default_rng(random_state)
    count = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(combined)
        perm_diff = abs(np.mean(permuted[:n1]) - np.mean(permuted[n1:]))
        if perm_diff >= observed_diff:
            count += 1

    p_value = count / n_permutations

    return {
        "test_statistic": observed_diff,
        "statistic": observed_diff,
        "p_value": p_value,
        "n_permutations": n_permutations,
        "significant": p_value < 0.05,
    }


def tajimas_d_outliers(d_values: List[float], window_size: int = 10, threshold: float = 3.0) -> Dict[str, Any]:
    """Detect outlier Tajima's D values using sliding window analysis.

    Args:
        d_values: List of Tajima's D values along a chromosome/sequence
        window_size: Size of sliding window for outlier detection

    Returns:
        Dictionary with outlier detection results
    """
    if len(d_values) < window_size:
        return detect_outliers(d_values, method="zscore", threshold=threshold)

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
                if z_score > threshold:
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
    return float(np.var(values, ddof=0)) if len(values) > 1 else 0.0


def standard_deviation(values: List[float]) -> float:
    """Calculate standard deviation of a list of values."""
    if not values:
        return 0.0
    return float(np.std(values, ddof=0)) if len(values) > 1 else 0.0


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


def normal_approximate_power(
    effect_size: float, sample_size_per_group: int, alpha: float = 0.05, two_sided: bool = True
) -> float:
    """Approximate two-sample power for detecting ``effect_size`` (in SD units).

    Uses the standard normal-approximation power formula for a two-sample
    z-test with equal group sizes. Intended as a descriptive planning aid for
    campaign design documents (e.g. RNA-seq group-size justification), not as
    a substitute for exact t-based power when df matter.

    Args:
        effect_size: Cohen's-d-style standardized effect (mean difference / SD).
        sample_size_per_group: Number of samples in each of the two groups.
        alpha: Significance level.
        two_sided: Two-sided alternative when True (default).

    Returns:
        Approximate statistical power in [0, 1].
    """
    from scipy.stats import norm

    if effect_size < 0:
        raise ValueError("effect_size must be non-negative; use |effect| with two_sided")
    if sample_size_per_group < 2:
        raise ValueError("sample_size_per_group must be >= 2")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0, 1)")

    if two_sided:
        z_alpha = float(norm.ppf(1.0 - alpha / 2.0))
    else:
        z_alpha = float(norm.ppf(1.0 - alpha))
    ncp = abs(effect_size) * math.sqrt(sample_size_per_group / 2.0)
    return float(norm.cdf(ncp - z_alpha))


def sample_size_for_power(
    effect_size: float, target_power: float = 0.8, alpha: float = 0.05, two_sided: bool = True
) -> int:
    """Smallest equal-per-group sample size reaching ``target_power`` (normal approx).

    Inverse of :func:`normal_approximate_power`: returns the minimal integer
    per-group sample size n such that the two-sample normal-approximate power
    is at least ``target_power``. Descriptive planning aid only.

    Args:
        effect_size: Cohen's-d-style standardized effect.
        target_power: Desired power in (alpha, 1).
        alpha: Significance level.
        two_sided: Two-sided alternative when True (default).

    Returns:
        Minimal per-group sample size (>= 2) achieving the target power.

    Raises:
        ValueError: If ``effect_size`` is 0 (power is capped at alpha for any
            finite n under the two-sided null) or arguments are out of range.
    """
    from scipy.stats import norm

    if effect_size <= 0:
        raise ValueError("effect_size must be > 0 to solve for sample size")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0, 1)")
    if not (alpha < target_power < 1.0):
        raise ValueError("target_power must be in (alpha, 1)")

    if two_sided:
        z_alpha = float(norm.ppf(1.0 - alpha / 2.0))
    else:
        z_alpha = float(norm.ppf(1.0 - alpha))
    z_power = float(norm.ppf(target_power))
    n_continuous = 2.0 * ((z_alpha + z_power) / effect_size) ** 2
    return max(2, int(math.ceil(n_continuous)))


def deterministic_replicate_seeds(base_seed: int, n_replicates: int) -> List[int]:
    """Stable, order-independent replicate seeds (canonical home: core.utils.seeds).

    Re-exported here for backwards compatibility; the implementation lives in
    ``metainformant.core.utils.seeds`` so other domains can import it without
    crossing the math domain boundary.
    """
    from metainformant.core.utils.seeds import deterministic_replicate_seeds as _impl

    return _impl(base_seed, n_replicates)
