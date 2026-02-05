"""Ecology analysis subpackage.

Provides community ecology, ordination, indicator species,
functional ecology, and macroecology analysis tools.
"""

from __future__ import annotations

from . import community, functional, indicators, macroecology, ordination

# Re-export community functions
from .community import (
    alpha_beta_gamma_diversity,
    beta_diversity,
    calculate_biodiversity_indices,
    calculate_diversity,
    calculate_evenness,
    calculate_single_diversity,
    chao1_estimator,
    community_metrics,
    community_similarity_matrix,
    dominance_diversity_curve,
    generate_ecology_report,
    nestedness_temperature_calculator,
    pielou_evenness,
    rank_abundance_curve,
    rarefaction_curve,
    shannon_diversity,
    simpson_diversity,
    species_accumulation_curve,
    species_area_relationship,
    species_richness,
    species_richness_simple,
)

# Re-export functional ecology functions
from .functional import (
    community_weighted_mean,
    functional_beta_diversity,
    functional_dispersion,
    functional_divergence,
    functional_diversity_suite,
    functional_evenness,
    functional_redundancy,
    functional_richness,
    raos_quadratic_entropy,
    trait_distance_matrix,
)

# Re-export indicator species functions
from .indicators import (
    anosim,
    cluster_communities,
    indval,
    multivariate_dispersion,
    permanova,
    simper,
)

# Re-export macroecology functions
from .macroecology import (
    compare_sad_models,
    distance_decay,
    endemism_index,
    fit_broken_stick,
    fit_geometric_series,
    fit_lognormal,
    fit_logseries,
    metabolic_scaling,
    occupancy_frequency,
    species_area_logarithmic,
    species_area_power,
    taylors_power_law,
)

# Re-export ordination functions
from .ordination import (
    cca,
    distance_matrix,
    nmds,
    pcoa,
    procrustes,
)

__all__ = [
    # Submodules
    "community",
    "ordination",
    "indicators",
    "functional",
    "macroecology",
    # Community ecology
    "shannon_diversity",
    "simpson_diversity",
    "species_richness",
    "species_richness_simple",
    "pielou_evenness",
    "chao1_estimator",
    "beta_diversity",
    "alpha_beta_gamma_diversity",
    "community_similarity_matrix",
    "rarefaction_curve",
    "species_accumulation_curve",
    "rank_abundance_curve",
    "dominance_diversity_curve",
    "nestedness_temperature_calculator",
    "species_area_relationship",
    "calculate_biodiversity_indices",
    "calculate_diversity",
    "calculate_single_diversity",
    "calculate_evenness",
    "community_metrics",
    "generate_ecology_report",
    # Ordination
    "distance_matrix",
    "pcoa",
    "nmds",
    "cca",
    "procrustes",
    # Indicator species
    "indval",
    "anosim",
    "permanova",
    "cluster_communities",
    "simper",
    "multivariate_dispersion",
    # Functional ecology
    "functional_richness",
    "functional_evenness",
    "functional_divergence",
    "functional_dispersion",
    "raos_quadratic_entropy",
    "community_weighted_mean",
    "functional_redundancy",
    "trait_distance_matrix",
    "functional_beta_diversity",
    "functional_diversity_suite",
    # Macroecology
    "fit_logseries",
    "fit_lognormal",
    "fit_broken_stick",
    "fit_geometric_series",
    "compare_sad_models",
    "species_area_power",
    "species_area_logarithmic",
    "distance_decay",
    "occupancy_frequency",
    "metabolic_scaling",
    "endemism_index",
    "taylors_power_law",
]
