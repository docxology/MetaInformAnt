"""Mathematical and theoretical biology utilities.

This subpackage provides quantitative biology primitives used across domains:
- Price equation and selection decomposition
- Kin and multilevel selection
- Drift–Diffusion models (DDM)
- Population genetics models
- Coalescent theory
- Epidemiological models
- Linkage disequilibrium
- Quantitative genetics

These are deliberately lightweight, dependency-minimal, and intended for
composition with `core` utilities and `dna`/`rna`/`protein` domain modules.

Submodules:
- coalescent: Coalescent theory and neutrality tests (Tajima's D, Fu & Li's tests)
- price: Price equation and selection metrics (covariance, selection gradients)
- popgen: Population genetics models (Hardy-Weinberg, selection, mutation)
- ld: Linkage disequilibrium (LD coefficients, r², mapping functions)
- epidemiology: Disease dynamics models (SIR, SEIR, SIS, R₀ calculations)
- dynamics: Population dynamics (logistic map, Lotka-Volterra)
- ddm: Decision-making models (drift-diffusion model)
- selection: Kin and multilevel selection (Hamilton's rule)
- quantgen: Quantitative genetics (heritability, breeder's equation)
- popgen_stats: Statistical testing (bootstrap, permutation tests, outliers)
- effective_size: Effective population size calculations
- demography: Demographic models (bottlenecks, growth, two-epoch)
- fst: Population differentiation (Fst calculations)
- egt: Evolutionary game theory (replicator dynamics)
- selection_experiments: Natural selection simulations

Utility Functions (in this module):
- correlation_coefficient: Pearson correlation coefficient
- linear_regression: Least-squares linear regression
- fisher_exact_test: Fisher's exact test for 2x2 contingency tables
- shannon_entropy: Information entropy
- jensen_shannon_divergence: Distribution divergence measure

Dependencies:
- numpy: Required for popgen_stats module
- scipy: Optional, used for advanced statistics (hardy_weinberg_test,
         fisher_exact_test, etc.). Available via optional dependency groups:
         pip install metainformant[scientific] or pip install scipy

Examples:
    >>> from metainformant.math import price_equation, tajimas_D
    >>> # Price equation analysis
    >>> fitness = [1.0, 1.2, 0.9]
    >>> traits = [0.2, 0.4, 0.1]
    >>> cov, trans, total = price_equation(fitness, traits)
    >>> # Correlation analysis
    >>> from metainformant.math import correlation_coefficient
    >>> r = correlation_coefficient([1, 2, 3], [2, 4, 6])
    >>> r
    1.0
"""

from .coalescent import (
    CoalescentSummary,
    expected_coalescent_waiting_times,
    expected_pairwise_diversity,
    expected_pairwise_diversity_from_theta,
    expected_segregating_sites,
    expected_sfs_counts,
    expected_time_to_mrca,
    expected_total_branch_length,
    ewens_watterson_test,
    fay_wu_h,
    fu_and_li_d_star,
    fu_and_li_f_star,
    hardy_weinberg_test,
    site_frequency_spectrum_counts,
    tajima_constants,
    tajimas_D,
    watterson_theta,
    wattersons_theta,
)
from .ddm import ddm_analytic_accuracy, ddm_mean_decision_time
from .dynamics import logistic_map, lotka_volterra_step
from .effective_size import (
    effective_size_from_family_size_variance,
    effective_size_sex_ratio,
    harmonic_mean_effective_size,
)
from .demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
    two_epoch_effective_size,
)
from .egt import replicator_derivative, replicator_step
from .epidemiology import (
    basic_reproduction_number,
    effective_reproduction_number,
    herd_immunity_threshold,
    seir_step,
    sir_step,
    sis_step,
)
from .fst import fst_from_allele_freqs, fst_from_heterozygosity
from .ld import (
    expected_r2_from_Ne_c,
    haldane_c_to_d,
    haldane_d_to_c,
    kosambi_c_to_d,
    kosambi_d_to_c,
    ld_coefficients,
    ld_decay_r2,
    r_squared,
)
from .popgen import (
    coalescent_time_to_mrca,
    effective_population_size_from_heterozygosity,
    equilibrium_heterozygosity_infinite_alleles,
    fixation_probability,
    hardy_weinberg_genotype_freqs,
    heterozygosity_decay,
    inbreeding_coefficient,
    inbreeding_coefficient_from_fst,
    island_model_update,
    linkage_disequilibrium_decay_distance,
    mutation_selection_balance_dominant,
    mutation_selection_balance_recessive,
    mutation_update,
    selection_update,
)
from .price import (
    correlation,
    covariance,
    delta_mean_trait,
    expectation,
    price_equation,
    relative_fitness,
    selection_differential,
    selection_gradient,
    selection_intensity,
    standard_deviation,
    variance,
    weighted_correlation,
    weighted_covariance,
    weighted_variance,
)
from .quantgen import (
    breeders_equation_response,
    lande_equation_response,
    narrow_sense_heritability,
    realized_heritability,
)
from .popgen_stats import (
    bootstrap_confidence_interval,
    calculate_confidence_intervals,
    compare_statistics,
    detect_outliers,
    permutation_test,
    tajimas_d_outliers,
    compare_population_statistic,
)
from .selection import kin_selection_response, multilevel_selection_decomposition

__all__ = [
    "price_equation",
    "covariance",
    "expectation",
    "variance",
    "correlation",
    "standard_deviation",
    "relative_fitness",
    "selection_differential",
    "selection_gradient",
    "selection_intensity",
    "delta_mean_trait",
    "kin_selection_response",
    "multilevel_selection_decomposition",
    "ddm_analytic_accuracy",
    "ddm_mean_decision_time",
    # population genetics
    "hardy_weinberg_genotype_freqs",
    "selection_update",
    "mutation_update",
    "fixation_probability",
    "heterozygosity_decay",
    "inbreeding_coefficient",
    "inbreeding_coefficient_from_fst",
    "equilibrium_heterozygosity_infinite_alleles",
    "island_model_update",
    "mutation_selection_balance_recessive",
    "mutation_selection_balance_dominant",
    "effective_population_size_from_heterozygosity",
    "linkage_disequilibrium_decay_distance",
    "coalescent_time_to_mrca",
    # linkage disequilibrium
    "ld_coefficients",
    "r_squared",
    "ld_decay_r2",
    "haldane_d_to_c",
    "haldane_c_to_d",
    "kosambi_d_to_c",
    "kosambi_c_to_d",
    "expected_r2_from_Ne_c",
    # coalescent
    "expected_time_to_mrca",
    "expected_total_branch_length",
    "expected_pairwise_diversity",
    "expected_pairwise_diversity_from_theta",
    "tajima_constants",
    "tajimas_D",
    "expected_segregating_sites",
    "expected_sfs_counts",
    "expected_coalescent_waiting_times",
    "site_frequency_spectrum_counts",
    "watterson_theta",
    "wattersons_theta",
    "CoalescentSummary",
    # quantitative genetics
    "narrow_sense_heritability",
    "breeders_equation_response",
    "lande_equation_response",
    "realized_heritability",
    # simple dynamics
    "logistic_map",
    "lotka_volterra_step",
    # evolutionary game theory
    "replicator_step",
    "replicator_derivative",
    # epidemiology
    "sir_step",
    "basic_reproduction_number",
    "seir_step",
    "herd_immunity_threshold",
    "sis_step",
    "effective_reproduction_number",
    # population structure
    "fst_from_heterozygosity",
    "fst_from_allele_freqs",
    # effective size
    "harmonic_mean_effective_size",
    "effective_size_sex_ratio",
    "effective_size_from_family_size_variance",
    # demography
    "exponential_growth_effective_size",
    "bottleneck_effective_size",
    "two_epoch_effective_size",
    # additional neutrality tests
    "fu_and_li_d_star",
    "fu_and_li_f_star",
    "fay_wu_h",
    "hardy_weinberg_test",
    "ewens_watterson_test",
    # statistical testing
    "bootstrap_confidence_interval",
    "permutation_test",
    "detect_outliers",
    "tajimas_d_outliers",
    "compare_statistics",
    "compare_population_statistic",
    "calculate_confidence_intervals",
    # utility functions
    "correlation_coefficient",
    "linear_regression",
    "fisher_exact_test",
    "shannon_entropy",
    "jensen_shannon_divergence",
]


# Statistical utilities
def correlation_coefficient(x: list[float], y: list[float]) -> float:
    """Calculate Pearson correlation coefficient between two lists.
    
    Computes the linear correlation between two variables, measuring the
    strength and direction of their linear relationship.
    
    Args:
        x: First list of numeric values
        y: Second list of numeric values (must match length of x)
        
    Returns:
        Pearson correlation coefficient in [-1, 1]. Returns 0.0 if:
        - Lists have different lengths
        - Lists have fewer than 2 elements
        - Either list has zero variance
        
    Raises:
        ValueError: If lists have different lengths (when len(x) != len(y))
        
    Examples:
        >>> correlation_coefficient([1, 2, 3], [2, 4, 6])
        1.0
        >>> correlation_coefficient([1, 2, 3], [3, 2, 1])
        -1.0
        
    Note:
        This is a convenience function. For more advanced correlation
        analysis, use the `correlation` function from `price.py`.
    """
    if len(x) != len(y):
        raise ValueError(f"Lists must have same length: {len(x)} != {len(y)}")
    if len(x) < 2:
        return 0.0
        
    import math
    
    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi * xi for xi in x)
    sum_y2 = sum(yi * yi for yi in y)
    
    numerator = n * sum_xy - sum_x * sum_y
    denominator = math.sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y))
    
    if denominator == 0:
        return 0.0
        
    return numerator / denominator


def linear_regression(x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Calculate linear regression: y = mx + b.
    
    Performs least-squares linear regression to fit a line through data points
    and calculates the goodness of fit.
    
    Args:
        x: Independent variable values
        y: Dependent variable values (must match length of x)
        
    Returns:
        Tuple of (slope, intercept, r_squared):
        - slope: Regression coefficient (m)
        - intercept: Y-axis intercept (b)
        - r_squared: Coefficient of determination (0 to 1)
        
        Returns (0.0, 0.0, 0.0) if:
        - Lists have fewer than 2 elements
        - X has zero variance
        
    Raises:
        ValueError: If lists have different lengths
        
    Examples:
        >>> slope, intercept, r2 = linear_regression([1, 2, 3], [2, 4, 6])
        >>> abs(slope - 2.0) < 0.01
        True
        >>> abs(r2 - 1.0) < 0.01  # Perfect fit
        True
    """
    if len(x) != len(y):
        raise ValueError(f"Lists must have same length: {len(x)} != {len(y)}")
    if len(x) < 2:
        return 0.0, 0.0, 0.0
        
    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi * xi for xi in x)
    sum_y2 = sum(yi * yi for yi in y)
    
    # Calculate slope and intercept
    denominator = n * sum_x2 - sum_x * sum_x
    if denominator == 0:
        return 0.0, 0.0, 0.0
        
    slope = (n * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y - slope * sum_x) / n
    
    # Calculate R-squared
    y_mean = sum_y / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    return slope, intercept, r_squared


try:
    from scipy.stats import fisher_exact as _scipy_fisher_exact
    
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    _scipy_fisher_exact = None


def fisher_exact_test(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Calculate Fisher's exact test for 2x2 contingency table.
    
    Computes odds ratio and p-value for a 2x2 contingency table using
    the hypergeometric distribution. Useful for testing independence
    between two categorical variables.
    
    Args:
        a: Top-left cell count
        b: Top-right cell count
        c: Bottom-left cell count
        d: Bottom-right cell count
        
    Returns:
        Tuple of (odds_ratio, p_value):
        - odds_ratio: (a × d) / (b × c). Returns inf if b or c is zero.
        - p_value: Two-tailed p-value from Fisher's exact test.
          Returns 1.0 if scipy is unavailable (fallback).
          
    Raises:
        ImportError: If scipy is not available. Install with:
            pip install metainformant[scientific]
            or
            pip install scipy
            
    Examples:
        >>> odds, p = fisher_exact_test(10, 2, 3, 15)
        >>> odds
        25.0
        >>> 0.0 <= p <= 1.0
        True
        
    Note:
        This function uses scipy.stats.fisher_exact when available for
        accurate p-value calculation. If scipy is not installed, the
        function will raise ImportError. For manual calculation without
        scipy, consider using scipy.stats.fisher_exact directly.
        
    References:
        Fisher, R. A. (1922). On the interpretation of χ² from contingency
        tables, and the calculation of P. Journal of the Royal Statistical
        Society, 85(1), 87-94.
    """
    import math
    
    # Calculate odds ratio
    odds_ratio = (a * d) / (b * c) if b > 0 and c > 0 else float('inf')
    
    # Calculate p-value using scipy if available
    if SCIPY_AVAILABLE and _scipy_fisher_exact is not None:
        # Create 2x2 contingency table
        table = [[a, b], [c, d]]
        odds_ratio_scipy, p_value = _scipy_fisher_exact(table, alternative='two-sided')
        # Use scipy's odds ratio for consistency
        return float(odds_ratio_scipy), float(p_value)
    else:
        # Raise ImportError if scipy not available
        raise ImportError(
            "scipy is required for fisher_exact_test(). "
            "Install with: pip install metainformant[scientific] or pip install scipy"
        )


def shannon_entropy(values: list[float]) -> float:
    """Calculate Shannon entropy of a probability distribution.
    
    Measures the information content or uncertainty in a probability distribution.
    Higher entropy indicates greater uncertainty/randomness.
    
    Args:
        values: List of probability values (should sum to 1, but function
            handles non-normalized inputs by only considering positive values)
        
    Returns:
        Shannon entropy in bits (using log base 2). Formula:
        H = -Σ p_i × log2(p_i) for p_i > 0
        
    Examples:
        >>> shannon_entropy([0.5, 0.5])  # Maximum entropy for 2 outcomes
        1.0
        >>> shannon_entropy([1.0, 0.0])  # Certainty (zero entropy)
        0.0
        >>> shannon_entropy([0.25, 0.25, 0.25, 0.25])  # Maximum for 4 outcomes
        2.0
        
    References:
        Shannon, C. E. (1948). A mathematical theory of communication.
        Bell System Technical Journal, 27(3), 379-423.
    """
    import math
    
    entropy = 0.0
    for p in values:
        if p > 0:
            entropy -= p * math.log2(p)
    
    return entropy


def jensen_shannon_divergence(p: list[float], q: list[float]) -> float:
    """Calculate Jensen-Shannon divergence between two probability distributions.
    
    JS divergence is a symmetrized version of Kullback-Leibler divergence,
    measuring the distance between two probability distributions. Always
    bounded in [0, 1] when using log base 2.
    
    Args:
        p: First probability distribution (list of probabilities)
        q: Second probability distribution (must have same length as p)
        
    Returns:
        Jensen-Shannon divergence. Returns 0.0 if distributions are identical
        or if either distribution sums to zero.
        Formula: JS(P||Q) = H(M) - [H(P) + H(Q)]/2 where M = (P+Q)/2
        
    Examples:
        >>> p = [0.5, 0.5]
        >>> q = [0.5, 0.5]
        >>> jensen_shannon_divergence(p, q)
        0.0  # Identical distributions
        
        >>> p = [1.0, 0.0]
        >>> q = [0.0, 1.0]
        >>> js = jensen_shannon_divergence(p, q)
        >>> js > 0.0  # Divergent distributions
        True
        
    Raises:
        ValueError: If distributions have different lengths
        
    References:
        Lin, J. (1991). Divergence measures based on the Shannon entropy.
        IEEE Transactions on Information Theory, 37(1), 145-151.
    """
    import math
    
    if len(p) != len(q):
        raise ValueError("Distributions must have same length")
    
    # Normalize distributions
    p_sum = sum(p)
    q_sum = sum(q)
    
    if p_sum == 0 or q_sum == 0:
        return 0.0
    
    p_norm = [pi / p_sum for pi in p]
    q_norm = [qi / q_sum for qi in q]
    
    # Calculate midpoint distribution
    m = [(pi + qi) / 2 for pi, qi in zip(p_norm, q_norm)]
    
    # Calculate entropies
    h_p = shannon_entropy(p_norm)
    h_q = shannon_entropy(q_norm)
    h_m = shannon_entropy(m)
    
    return h_m - (h_p + h_q) / 2
