"""Mathematical and theoretical biology utilities.

This subpackage provides quantitative biology primitives used across domains:
- Price equation and selection decomposition
- Kin and multilevel selection toy calculators
- Drift–Diffusion models (DDM)

These are deliberately lightweight, dependency-minimal, and intended for
composition with `core` utilities and `dna`/`rna`/`protein` domain modules.
"""

from .price import (
    price_equation,
    covariance,
    expectation,
    variance,
    correlation,
    standard_deviation,
    weighted_variance,
    weighted_covariance,
    weighted_correlation,
    relative_fitness,
    selection_differential,
    selection_gradient,
    selection_intensity,
    delta_mean_trait,
)
from .selection import kin_selection_response, multilevel_selection_decomposition
from .ddm import ddm_analytic_accuracy, ddm_mean_decision_time
from .popgen import (
    hardy_weinberg_genotype_freqs,
    selection_update,
    mutation_update,
    fixation_probability,
    watterson_theta,
    heterozygosity_decay,
    inbreeding_coefficient,
    equilibrium_heterozygosity_infinite_alleles,
    island_model_update,
    mutation_selection_balance_recessive,
    mutation_selection_balance_dominant,
)
from .ld import ld_coefficients, r_squared, ld_decay_r2
from .ld import (
    haldane_d_to_c,
    haldane_c_to_d,
    kosambi_d_to_c,
    kosambi_c_to_d,
    expected_r2_from_Ne_c,
)
from .coalescent import (
    expected_time_to_mrca,
    expected_total_branch_length,
    expected_pairwise_diversity,
    tajima_constants,
    tajimas_D,
    expected_segregating_sites,
)
from .quantgen import (
    narrow_sense_heritability,
    breeders_equation_response,
    lande_equation_response,
    realized_heritability,
)
from .dynamics import logistic_map, lotka_volterra_step
from .egt import replicator_step, replicator_derivative
from .epidemiology import sir_step, basic_reproduction_number, seir_step, herd_immunity_threshold, sis_step, effective_reproduction_number
from .fst import fst_from_heterozygosity, fst_from_allele_freqs
from .effective_size import harmonic_mean_effective_size, effective_size_sex_ratio, effective_size_from_family_size_variance

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


