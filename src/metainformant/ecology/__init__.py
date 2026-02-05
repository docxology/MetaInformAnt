"""Ecology and community analysis module for METAINFORMANT.

This module provides tools for ecological analysis, including community
structure, diversity metrics, ordination, indicator species, functional
ecology, macroecology, and ecological visualization.

Usage:
    from metainformant.ecology import shannon_diversity, simpson_diversity
    from metainformant.ecology import pcoa, distance_matrix
    from metainformant.ecology import indval, anosim
    from metainformant.ecology import functional_richness, raos_quadratic_entropy
    from metainformant.ecology import fit_logseries, compare_sad_models
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import visualization

# Re-export all analysis functions at package level for convenient access
from .analysis import (
    # Submodules
    community,
    ordination,
    indicators,
    functional,
    macroecology,
    # Community ecology
    shannon_diversity,
    simpson_diversity,
    species_richness,
    species_richness_simple,
    pielou_evenness,
    chao1_estimator,
    beta_diversity,
    alpha_beta_gamma_diversity,
    community_similarity_matrix,
    rarefaction_curve,
    species_accumulation_curve,
    rank_abundance_curve,
    dominance_diversity_curve,
    nestedness_temperature_calculator,
    species_area_relationship,
    calculate_biodiversity_indices,
    calculate_diversity,
    calculate_single_diversity,
    calculate_evenness,
    community_metrics,
    generate_ecology_report,
    # Ordination
    distance_matrix,
    pcoa,
    nmds,
    cca,
    procrustes,
    # Indicator species
    indval,
    anosim,
    permanova,
    cluster_communities,
    simper,
    multivariate_dispersion,
    # Functional ecology
    functional_richness,
    functional_evenness,
    functional_divergence,
    functional_dispersion,
    raos_quadratic_entropy,
    community_weighted_mean,
    functional_redundancy,
    trait_distance_matrix,
    functional_beta_diversity,
    functional_diversity_suite,
    # Macroecology
    fit_logseries,
    fit_lognormal,
    fit_broken_stick,
    fit_geometric_series,
    compare_sad_models,
    species_area_power,
    species_area_logarithmic,
    distance_decay,
    occupancy_frequency,
    metabolic_scaling,
    endemism_index,
    taylors_power_law,
)

# Re-export visualization functions
from .visualization.visualization import (
    plot_diversity_indices_comparison,
    plot_community_composition,
    plot_rank_abundance_curve_comparison,
    plot_species_abundance_distribution,
    plot_beta_diversity_ordination,
    plot_biodiversity_rarefaction,
    plot_ecological_distance_heatmap,
    plot_ecological_network,
    plot_diversity_accumulation_curve,
    create_interactive_ecology_dashboard,
)

__all__ = [
    # Subpackages
    "analysis",
    "visualization",
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
    # Visualization
    "plot_diversity_indices_comparison",
    "plot_community_composition",
    "plot_rank_abundance_curve_comparison",
    "plot_species_abundance_distribution",
    "plot_beta_diversity_ordination",
    "plot_biodiversity_rarefaction",
    "plot_ecological_distance_heatmap",
    "plot_ecological_network",
    "plot_diversity_accumulation_curve",
    "create_interactive_ecology_dashboard",
]
