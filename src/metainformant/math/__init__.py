"""Mathematical biology and theoretical modeling module for METAINFORMANT.

This module provides mathematical and computational tools for biological modeling.
"""

from __future__ import annotations

# Subpackages
from . import (
    bayesian,
    decision_theory,
    epidemiology,
    evolutionary_dynamics,
    population_genetics,
    quantitative_genetics,
)

# Create alias for backward compatibility
popgen = population_genetics
# Bayesian inference
from .bayesian.inference import (
    abc_rejection,
    compute_bayes_factor,
    compute_dic,
    compute_waic,
    conjugate_beta_binomial,
    conjugate_normal,
    metropolis_hastings,
)

# Re-export core utilities for backward compatibility
from .core.utilities import (
    correlation,
    correlation_coefficient,
    covariance,
    fisher_exact_test,
    jensen_shannon_divergence,
    linear_regression,
    r_squared,
    shannon_entropy,
)

# Re-export Decision Theory
from .decision_theory.ddm import (
    ddm_analytic_accuracy,
    ddm_mean_decision_time,
)
from .epidemiology import models as epimodels

# Re-export Epidemiology
from .epidemiology.models import (
    basic_reproduction_number,
    effective_reproduction_number,
    herd_immunity_threshold,
    seir_step,
    sir_step,
    sis_step,
)

# Re-export Evolutionary Dynamics
from .evolutionary_dynamics.core import (
    logistic_map,
    lotka_volterra_step,
    replicator_derivative,
    replicator_step,
)

# Re-export coalescent from population_genetics for backward compatibility
from .population_genetics import coalescent, demography
from .population_genetics import statistics as popgen_stats
from .population_genetics.coalescent import (
    expected_pairwise_diversity,
    expected_segregating_sites,
    expected_time_to_mrca,
    expected_total_branch_length,
    tajima_constants,
    tajimas_D,
    watterson_theta,
)

# Re-export core functions
from .population_genetics.core import (
    hardy_weinberg_genotype_freqs,
    heterozygosity_decay,
    inbreeding_coefficient,
)
from .population_genetics.demography import island_model_update
from .population_genetics.effective_size import (
    effective_size_sex_ratio,
    harmonic_mean_effective_size,
)
from .population_genetics.fst import (
    fst_from_allele_freqs,
    fst_from_heterozygosity,
)

# Re-export LD functions
from .population_genetics.ld import (
    haldane_c_to_d,
    haldane_d_to_c,
    kosambi_c_to_d,
    kosambi_d_to_c,
    ld_coefficients,
    ld_decay_r2,
)
from .population_genetics.selection import (
    breeders_equation_response,
    kin_selection_response,
    mutation_selection_balance_dominant,
    mutation_selection_balance_recessive,
    mutation_update,
    relative_fitness,
    selection_differential,
    selection_gradient,
    selection_intensity,
    selection_update,
)
from .population_genetics.statistics import (
    effective_size_from_family_size_variance,
    equilibrium_heterozygosity_infinite_alleles,
    expected_r2_from_Ne_c,
    fixation_probability,
    kurtosis,
    skewness,
    standard_deviation,
    variance,
)
from .quantitative_genetics.core import (
    lande_equation_response,
    narrow_sense_heritability,
    realized_heritability,
)
from .quantitative_genetics.price import (
    delta_mean_trait,
    expectation,
    price_equation,
    weighted_correlation,
    weighted_covariance,
    weighted_variance,
)

__all__ = [
    "coalescent",
    "decision_theory",
    "epidemiology",
    "evolutionary_dynamics",
    "population_genetics",
    "quantitative_genetics",
    "demography",
    "popgen",
    "popgen_stats",
    "epimodels",
    # Coalescent
    "expected_time_to_mrca",
    "watterson_theta",
    "expected_total_branch_length",
    "expected_pairwise_diversity",
    "tajimas_D",
    "tajima_constants",
    "expected_segregating_sites",
    # Core
    "hardy_weinberg_genotype_freqs",
    "fisher_exact_test",
    "correlation",
    # PopGen
    "fixation_probability",
    "expected_r2_from_Ne_c",
    "mutation_update",
    "selection_update",
    "kin_selection_response",
    "harmonic_mean_effective_size",
    "island_model_update",
    "heterozygosity_decay",
    "inbreeding_coefficient",
    "fst_from_allele_freqs",
    "fst_from_heterozygosity",
    # QuantGen
    "breeders_equation_response",
    "realized_heritability",
    "lande_equation_response",
    "narrow_sense_heritability",
    "price_equation",
    "delta_mean_trait",
    "expectation",
    "weighted_correlation",
    "weighted_covariance",
    "weighted_variance",
    # Utilities
    "covariance",
    "shannon_entropy",
    "jensen_shannon_divergence",
    "r_squared",
    "correlation_coefficient",
    "linear_regression",
    "relative_fitness",
    "selection_intensity",
    "selection_differential",
    "selection_gradient",
    "mutation_selection_balance_dominant",
    "mutation_selection_balance_recessive",
    # Epidemiology
    "basic_reproduction_number",
    "effective_reproduction_number",
    "herd_immunity_threshold",
    "kurtosis",
    "skewness",
    "standard_deviation",
    "variance",
    # PopGen Extensions
    "effective_size_from_family_size_variance",
    "equilibrium_heterozygosity_infinite_alleles",
    "effective_size_sex_ratio",
    # Evolutionary Dynamics
    "logistic_map",
    "lotka_volterra_step",
    "replicator_derivative",
    "replicator_step",
    # Decision Theory
    "ddm_analytic_accuracy",
    "ddm_mean_decision_time",
    # LD
    "ld_coefficients",
    "ld_decay_r2",
    "haldane_c_to_d",
    "haldane_d_to_c",
    "kosambi_c_to_d",
    "kosambi_d_to_c",
    # Bayesian
    "bayesian",
    "metropolis_hastings",
    "abc_rejection",
    "compute_bayes_factor",
    "conjugate_beta_binomial",
    "conjugate_normal",
    "compute_dic",
    "compute_waic",
]
