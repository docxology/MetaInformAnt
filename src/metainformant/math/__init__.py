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
    equilibrium_heterozygosity_infinite_alleles,
    fixation_probability,
    hardy_weinberg_genotype_freqs,
    heterozygosity_decay,
    inbreeding_coefficient,
    island_model_update,
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
    "equilibrium_heterozygosity_infinite_alleles",
    "island_model_update",
    "mutation_selection_balance_recessive",
    "mutation_selection_balance_dominant",
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
