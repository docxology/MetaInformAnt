"""Population genetics visualization utilities.

This module provides plotting and visualization functions for population
genetics data including F_ST heatmaps, selection plots, and demographic visualizations.

This module re-exports all public symbols from :mod:`visualization_core` and
:mod:`visualization_stats` for backward compatibility.
"""

from __future__ import annotations

from metainformant.dna.population.visualization_core import (
    create_population_summary_plot,
    plot_allele_frequency_spectrum,
    plot_bootstrap_distribution,
    plot_demographic_history,
    plot_fst_matrix,
    plot_ld_decay,
    plot_mutation_spectrum,
    plot_population_diversity,
    plot_population_structure,
    plot_selection_statistics,
    plot_tajima_d_distribution,
)
from metainformant.dna.population.visualization_stats import (
    plot_demographic_comparison,
    plot_diversity_comparison,
    plot_fst_comparison,
    plot_hardy_weinberg_test,
    plot_heterozygosity_distribution,
    plot_kinship_matrix,
    plot_linkage_disequilibrium_decay,
    plot_neutrality_test_suite,
    plot_neutrality_test_summary,
    plot_outlier_detection,
    plot_pca_results,
    plot_permutation_test,
    plot_pi_vs_theta,
    plot_site_frequency_spectrum,
    plot_statistic_correlation_matrix,
    plot_statistic_distribution,
    plot_summary_statistics_grid,
    plot_tajimas_d_comparison,
)

__all__ = [
    # Core visualization functions
    "plot_fst_matrix",
    "plot_tajima_d_distribution",
    "plot_selection_statistics",
    "plot_population_diversity",
    "plot_ld_decay",
    "plot_population_structure",
    "plot_demographic_history",
    "create_population_summary_plot",
    "plot_mutation_spectrum",
    "plot_allele_frequency_spectrum",
    "plot_bootstrap_distribution",
    # Statistical comparison functions
    "plot_demographic_comparison",
    "plot_diversity_comparison",
    "plot_fst_comparison",
    "plot_hardy_weinberg_test",
    "plot_heterozygosity_distribution",
    "plot_kinship_matrix",
    "plot_linkage_disequilibrium_decay",
    "plot_neutrality_test_suite",
    "plot_neutrality_test_summary",
    "plot_outlier_detection",
    "plot_pca_results",
    "plot_permutation_test",
    "plot_pi_vs_theta",
    "plot_site_frequency_spectrum",
    "plot_statistic_correlation_matrix",
    "plot_statistic_distribution",
    "plot_summary_statistics_grid",
    "plot_tajimas_d_comparison",
]
