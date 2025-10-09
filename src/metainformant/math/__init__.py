"""Mathematical and theoretical biology utilities.

This subpackage provides quantitative biology primitives used across domains:
- Price equation and selection decomposition
- Kin and multilevel selection toy calculators
- Drift–Diffusion models (DDM)

These are deliberately lightweight, dependency-minimal, and intended for
composition with `core` utilities and `dna`/`rna`/`protein` domain modules.
"""

from .coalescent import (
    expected_pairwise_diversity,
    expected_segregating_sites,
    expected_time_to_mrca,
    expected_total_branch_length,
    tajima_constants,
    tajimas_D,
)
from .ddm import ddm_analytic_accuracy, ddm_mean_decision_time
from .dynamics import logistic_map, lotka_volterra_step
from .effective_size import (
    effective_size_from_family_size_variance,
    effective_size_sex_ratio,
    harmonic_mean_effective_size,
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
    watterson_theta,
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
    "watterson_theta",
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
    "tajima_constants",
    "tajimas_D",
    "expected_segregating_sites",
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
]


# Statistical utilities
def correlation_coefficient(x: list[float], y: list[float]) -> float:
    """Calculate Pearson correlation coefficient between two lists.
    
    Args:
        x: First list of values
        y: Second list of values
        
    Returns:
        Correlation coefficient (-1 to 1)
    """
    if len(x) != len(y) or len(x) < 2:
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
    
    Args:
        x: Independent variable values
        y: Dependent variable values
        
    Returns:
        Tuple of (slope, intercept, r_squared)
    """
    if len(x) != len(y) or len(x) < 2:
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


def fisher_exact_test(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Calculate Fisher's exact test for 2x2 contingency table.
    
    Args:
        a: Top-left cell
        b: Top-right cell  
        c: Bottom-left cell
        d: Bottom-right cell
        
    Returns:
        Tuple of (odds_ratio, p_value)
    """
    import math
    
    # Calculate odds ratio
    odds_ratio = (a * d) / (b * c) if b > 0 and c > 0 else float('inf')
    
    # Calculate p-value using hypergeometric distribution
    # This is a simplified calculation
    total = a + b + c + d
    p_value = 1.0  # Placeholder - full implementation would use scipy.stats.fisher_exact
    
    return odds_ratio, p_value


def shannon_entropy(values: list[float]) -> float:
    """Calculate Shannon entropy of a probability distribution.
    
    Args:
        values: List of probability values (should sum to 1)
        
    Returns:
        Shannon entropy value
    """
    import math
    
    entropy = 0.0
    for p in values:
        if p > 0:
            entropy -= p * math.log2(p)
    
    return entropy


def jensen_shannon_divergence(p: list[float], q: list[float]) -> float:
    """Calculate Jensen-Shannon divergence between two probability distributions.
    
    Args:
        p: First probability distribution
        q: Second probability distribution
        
    Returns:
        Jensen-Shannon divergence
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
