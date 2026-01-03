"""Mathematical biology and theoretical modeling module for METAINFORMANT.

This module provides mathematical and computational tools for biological modeling,
including population genetics theory, evolutionary dynamics, epidemiology,
coalescent theory, and quantitative genetics.
"""

from __future__ import annotations

# Import all mathematical biology submodules
from . import (
    coalescent,
    ddm,
    demography,
    dynamics,
    effective_size,
    egt,
    epidemiology,
    fst,
    ld,
    popgen,
    popgen_stats,
    price,
    quantgen,
    selection,
    utilities,
    visualization,
)

# Import selection_experiments submodule
from . import selection_experiments

# Direct imports of commonly used functions
from .coalescent import expected_time_to_mrca, watterson_theta, expected_total_branch_length, tajima_constants, tajimas_D
from .ddm import ddm_analytic_accuracy, ddm_mean_decision_time, island_model_update
from .demography import exponential_growth_model, logistic_growth_model, age_structure_model
from .dynamics import logistic_map, lotka_volterra_step
from .epidemiology import basic_reproduction_number, effective_reproduction_number, herd_immunity_threshold, sir_step, seir_step, sis_step
from .egt import replicator_derivative, replicator_step
from .ld import ld_coefficients, ld_decay_r2, haldane_c_to_d, haldane_d_to_c, kosambi_c_to_d, kosambi_d_to_c
from .popgen import hardy_weinberg_genotype_freqs, heterozygosity_decay, inbreeding_coefficient, mutation_selection_balance_dominant, mutation_selection_balance_recessive
from .popgen_stats import expected_pairwise_diversity, expected_segregating_sites, expected_coalescent_waiting_times, expected_r2_from_Ne_c, equilibrium_heterozygosity_infinite_alleles, fixation_probability, bottleneck_effective_size, effective_size_from_family_size_variance, bootstrap_confidence_interval, calculate_confidence_intervals
from .selection import breeders_equation_response, kin_selection_response, mutation_update, selection_update, selection_differential, selection_gradient
from .price import (price_equation, delta_mean_trait, expectation, relative_fitness,
                    selection_intensity, standard_deviation, variance, weighted_variance,
                    weighted_covariance, weighted_correlation)
from .effective_size import effective_size_sex_ratio, harmonic_mean_effective_size
from .quantgen import narrow_sense_heritability, realized_heritability, lande_equation_response
from .selection import multilevel_selection_decomposition
from .utilities import correlation, correlation_coefficient, linear_regression, r_squared, fisher_exact_test, covariance, shannon_entropy, jensen_shannon_divergence
from .fst import fst_from_allele_freqs, fst_from_heterozygosity

# Optional imports with graceful fallbacks
try:
    from . import selection_experiments
except ImportError:
    selection_experiments = None

try:
    from . import ddm
except ImportError:
    ddm = None

try:
    from . import egt
except ImportError:
    egt = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Population genetics
    "popgen",
    "popgen_stats",
    "fst",
    "ld",
    "effective_size",

    # Evolutionary theory
    "coalescent",
    "selection",
    "price",

    # Population dynamics and ecology
    "demography",
    "dynamics",
    "epidemiology",

    # Decision theory and behavior
    "ddm",
    "island_model_update",
    "egt",

    # Quantitative genetics
    "quantgen",

    # Specialized experiments
    "selection_experiments",

    # Visualization
    "visualization",

    # Direct function exports
    "expected_time_to_mrca",
    "watterson_theta",
    "expected_total_branch_length",
    "tajima_constants",
    "tajimas_D",
    "exponential_growth_model",
    "logistic_growth_model",
    "age_structure_model",
    "logistic_map",
    "lotka_volterra_step",
    "basic_reproduction_number",
    "effective_reproduction_number",
    "herd_immunity_threshold",
    "sir_step",
    "seir_step",
    "sis_step",
    "ld_coefficients",
    "ld_decay_r2",
    "haldane_c_to_d",
    "haldane_d_to_c",
    "hardy_weinberg_genotype_freqs",
    "mutation_selection_balance_dominant",
    "mutation_selection_balance_recessive",
    "heterozygosity_decay",
    "inbreeding_coefficient",
    "expected_pairwise_diversity",
    "expected_segregating_sites",
    "expected_coalescent_waiting_times",
    "expected_r2_from_Ne_c",
    "equilibrium_heterozygosity_infinite_alleles",
    "fixation_probability",
    "bottleneck_effective_size",
    "fst_from_heterozygosity",
    "harmonic_mean_effective_size",
    "replicator_derivative",
    "replicator_step",
    "effective_size_from_family_size_variance",
    "bootstrap_confidence_interval",
    "breeders_equation_response",
    "kin_selection_response",
    "mutation_update",
    "selection_update",
    "selection_differential",
    "selection_gradient",
    "price_equation",
    "delta_mean_trait",
    "expectation",
    "relative_fitness",
    "selection_intensity",
    "standard_deviation",
    "variance",
    "weighted_variance",
    "weighted_covariance",
    "weighted_correlation",
    "effective_size_sex_ratio",
    "narrow_sense_heritability",
    "realized_heritability",
    "multilevel_selection_decomposition",
    "lande_equation_response",
    "correlation",
    "correlation_coefficient",
    "linear_regression",
    "r_squared",
    "shannon_entropy",
    "jensen_shannon_divergence",
]



