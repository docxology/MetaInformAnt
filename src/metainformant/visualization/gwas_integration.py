"""GWAS visualization integration.

This module provides unified access to GWAS visualization functions
from the gwas module, integrating them into the central visualization package.
"""

from __future__ import annotations

# Re-export GWAS visualization functions
try:
    from ..gwas.visualization import manhattan_plot, qq_plot, regional_plot
    from ..gwas.visualization_genome import (
        circular_manhattan_plot,
        chromosome_ideogram,
        genome_wide_ld_heatmap,
        manhattan_plot as manhattan_plot_genome,
    )
    from ..gwas.visualization_statistical import (
        lambda_gc_plot,
        power_plot,
        qq_plot_stratified,
        volcano_plot as volcano_plot_gwas,
    )
    from ..gwas.visualization_regional import (
        gene_annotation_plot,
        recombination_rate_plot,
        regional_ld_plot,
        regional_plot as regional_plot_detailed,
    )
    from ..gwas.visualization_population import (
        admixture_plot,
        kinship_heatmap,
        pca_plot as pca_plot_gwas,
        pca_scree_plot as pca_scree_plot_gwas,
        population_tree,
    )
    from ..gwas.visualization_variants import (
        allelic_series_plot,
        hwe_deviation_plot,
        maf_distribution,
        missingness_plot,
        transition_transversion_plot,
        variant_density_plot,
    )
    from ..gwas.visualization_effects import (
        effect_direction_plot,
        effect_size_forest_plot,
        functional_enrichment_plot,
    )
    from ..gwas.visualization_comparison import (
        miami_plot,
        multi_trait_manhattan,
    )

    GWAS_VISUALIZATION_AVAILABLE = True
except ImportError:
    GWAS_VISUALIZATION_AVAILABLE = False
    # Define placeholder functions
    def _not_available(*args, **kwargs):
        raise ImportError("GWAS visualization functions not available. Install gwas module dependencies.")
    
    manhattan_plot = _not_available
    qq_plot = _not_available
    regional_plot = _not_available
    circular_manhattan_plot = _not_available
    chromosome_ideogram = _not_available
    genome_wide_ld_heatmap = _not_available
    manhattan_plot_genome = _not_available
    lambda_gc_plot = _not_available
    power_plot = _not_available
    qq_plot_stratified = _not_available
    volcano_plot_gwas = _not_available
    gene_annotation_plot = _not_available
    recombination_rate_plot = _not_available
    regional_ld_plot = _not_available
    regional_plot_detailed = _not_available
    admixture_plot = _not_available
    kinship_heatmap = _not_available
    pca_plot_gwas = _not_available
    pca_scree_plot_gwas = _not_available
    population_tree = _not_available
    allelic_series_plot = _not_available
    hwe_deviation_plot = _not_available
    maf_distribution = _not_available
    missingness_plot = _not_available
    transition_transversion_plot = _not_available
    variant_density_plot = _not_available
    effect_direction_plot = _not_available
    effect_size_forest_plot = _not_available
    functional_enrichment_plot = _not_available
    miami_plot = _not_available
    multi_trait_manhattan = _not_available

__all__ = [
    # Basic GWAS plots
    "manhattan_plot",
    "qq_plot",
    "regional_plot",
    # Genome-wide visualizations
    "circular_manhattan_plot",
    "chromosome_ideogram",
    "genome_wide_ld_heatmap",
    "manhattan_plot_genome",
    # Statistical diagnostics
    "lambda_gc_plot",
    "power_plot",
    "qq_plot_stratified",
    "volcano_plot_gwas",
    # Regional plots
    "gene_annotation_plot",
    "recombination_rate_plot",
    "regional_ld_plot",
    "regional_plot_detailed",
    # Population structure
    "admixture_plot",
    "kinship_heatmap",
    "pca_plot_gwas",
    "pca_scree_plot_gwas",
    "population_tree",
    # Variant properties
    "allelic_series_plot",
    "hwe_deviation_plot",
    "maf_distribution",
    "missingness_plot",
    "transition_transversion_plot",
    "variant_density_plot",
    # Effect sizes
    "effect_direction_plot",
    "effect_size_forest_plot",
    "functional_enrichment_plot",
    # Comparison plots
    "miami_plot",
    "multi_trait_manhattan",
    # Availability flag
    "GWAS_VISUALIZATION_AVAILABLE",
]

