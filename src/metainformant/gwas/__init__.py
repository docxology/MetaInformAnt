"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

This module provides comprehensive GWAS analysis capabilities, including
quality control, population structure analysis, association testing,
multiple testing correction, and publication-quality visualization.
"""

from __future__ import annotations

# Import subpackages (making them available at package level)
from . import data  # Configuration and data acquisition
from . import visualization  # Visualization tools
from . import analysis  # Core analysis logic (depends on data)
from . import workflow  # Workflow orchestration (depends on all)
from . import finemapping  # Statistical fine-mapping (credible sets, coloc, SuSiE)
from . import heritability as heritability_submodule  # Heritability estimation (LDSC, GREML)

# Import modules from subpackages for backward compatibility and ease of use
from .analysis import (
    annotation,
    association,
    calling,
    correction,
    heritability,
    ld_pruning,
    mixed_model,
    quality,
    structure,
    summary_stats,
)
from .data import (
    config,
    download,
    genome,
    metadata,
    sra_download,
)
from .workflow import workflow as workflow_module

# Visualization modules
# Note: 'visualization' module was renamed to 'general' inside the subpackage
from .visualization import general as visualization_module
from .visualization import (
    visualization_comparison,
    visualization_composite,
    visualization_effects,
    visualization_finemapping,
    visualization_genome,
    visualization_geography,
    visualization_interactive,
    visualization_ld,
    visualization_phenotype,
    visualization_population,
    visualization_regional,
    visualization_statistical,
    visualization_suite,
    visualization_variants,
)

# Direct imports of commonly used functions
from .analysis.annotation import annotate_variants_with_genes, classify_variant_location
from .analysis.association import association_test_linear, association_test_logistic
from .analysis.calling import check_bcftools_available
from .analysis.correction import bonferroni_correction, fdr_correction, genomic_control
from .analysis.ld_pruning import ld_prune
from .analysis.mixed_model import association_test_mixed, run_mixed_model_gwas
from .analysis.quality import apply_qc_filters, parse_vcf_full, extract_variant_regions
from .analysis.structure import compute_kinship_matrix, compute_pca, estimate_population_structure
from .analysis.summary_stats import write_summary_statistics, write_significant_hits, create_results_summary
from .analysis.heritability import estimate_heritability, partition_heritability_by_chromosome, heritability_bar_chart

from .data.config import load_gwas_config
from .data.metadata import load_sample_metadata, validate_metadata, get_population_labels, get_geographic_coordinates
from .data.download import download_reference_genome, download_variant_data
from .data.genome import normalize_chromosome_name, AMEL_HAV3_CHROMOSOMES, AMEL_HAV3_CHROM_SIZES

# Visualization aliases
from .visualization.general import qq_plot, manhattan_plot, kinship_heatmap
from .visualization.visualization_variants import variant_density_plot
from .visualization.visualization_population import (
    pca_scree_plot,
    admixture_plot,
    kinship_dendrogram,
    kinship_clustermap,
)
from .visualization.visualization_comparison import multi_trait_manhattan
from .visualization.visualization_effects import functional_enrichment_plot
from .visualization.visualization_statistical import power_plot
from .visualization.visualization_regional import (
    recombination_rate_plot,
    effect_direction_plot,
    regional_plot,
    regional_ld_plot,
)
from .visualization.visualization_genome import genome_wide_ld_heatmap
from .visualization.visualization_ld import compute_ld_decay, ld_decay_plot, ld_heatmap_region
from .visualization.visualization_phenotype import (
    phenotype_distribution,
    phenotype_correlation_matrix,
    genotype_phenotype_boxplot,
    phenotype_pca_correlation,
)
from .visualization.visualization_composite import gwas_summary_panel, population_structure_panel
from .visualization.visualization_geography import sample_map, allele_frequency_map, population_count_map
from .visualization.visualization_finemapping import credible_set_plot, compute_credible_set, pip_vs_ld_plot
from .visualization.visualization_interactive import interactive_manhattan, interactive_pca, interactive_volcano
from .visualization.visualization_population import pca_multi_panel, pca_3d
from .visualization.config import PlotStyle, THEMES, get_style, apply_style

from .workflow.workflow import GWASWorkflowConfig, execute_gwas_workflow, run_gwas

# Fine-mapping submodule imports
from .finemapping.credible_sets import (
    compute_credible_set as finemapping_credible_set,
    susie_regression,
    compute_bayes_factors,
    colocalization,
    conditional_analysis,
    annotate_credible_set,
)
from .finemapping.colocalization import (
    multi_trait_coloc,
    eqtl_coloc,
    compute_clpp,
    regional_coloc,
)

# Heritability submodule imports
from .heritability.estimation import (
    estimate_h2_ldsc,
    partitioned_h2,
    genetic_correlation,
    haseman_elston_regression,
    greml_simple,
    compute_liability_h2,
)

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "data",
    "visualization",
    "workflow",
    # Core GWAS analysis
    "apply_qc_filters",
    "parse_vcf_full",
    "quality",
    "extract_variant_regions",
    "structure",
    "compute_pca",
    "compute_kinship_matrix",
    "estimate_population_structure",
    "association",
    "correction",
    # New analysis modules
    "annotation",
    "annotate_variants_with_genes",
    "classify_variant_location",
    "ld_pruning",
    "ld_prune",
    "mixed_model",
    "association_test_mixed",
    "run_mixed_model_gwas",
    "summary_stats",
    "write_summary_statistics",
    "write_significant_hits",
    "create_results_summary",
    # Genome mapping
    "genome",
    "normalize_chromosome_name",
    "AMEL_HAV3_CHROMOSOMES",
    "AMEL_HAV3_CHROM_SIZES",
    # Data acquisition
    "download",
    "download_variant_data",
    "sra_download",
    "calling",
    "config",
    # Visualization suite
    "visualization_module",
    "visualization_comparison",
    "visualization_effects",
    "visualization_genome",
    "visualization_population",
    "visualization_regional",
    "visualization_statistical",
    "visualization_suite",
    "visualization_variants",
    "variant_density_plot",
    "pca_scree_plot",
    "multi_trait_manhattan",
    "functional_enrichment_plot",
    "power_plot",
    "recombination_rate_plot",
    "qq_plot",
    "kinship_heatmap",
    "manhattan_plot",
    "effect_direction_plot",
    "regional_plot",
    "regional_ld_plot",
    "genome_wide_ld_heatmap",
    "admixture_plot",
    "kinship_dendrogram",
    "kinship_clustermap",
    "pca_multi_panel",
    "pca_3d",
    # New visualization modules
    "visualization_composite",
    "visualization_finemapping",
    "visualization_geography",
    "visualization_interactive",
    "visualization_ld",
    "visualization_phenotype",
    # LD decay
    "compute_ld_decay",
    "ld_decay_plot",
    "ld_heatmap_region",
    # Phenotype plots
    "phenotype_distribution",
    "phenotype_correlation_matrix",
    "genotype_phenotype_boxplot",
    "phenotype_pca_correlation",
    # Composite figures
    "gwas_summary_panel",
    "population_structure_panel",
    # Geography
    "sample_map",
    "allele_frequency_map",
    "population_count_map",
    # Fine-mapping
    "credible_set_plot",
    "compute_credible_set",
    "pip_vs_ld_plot",
    # Interactive
    "interactive_manhattan",
    "interactive_pca",
    "interactive_volcano",
    # Viz config
    "PlotStyle",
    "THEMES",
    "get_style",
    "apply_style",
    # Heritability
    "heritability",
    "estimate_heritability",
    "partition_heritability_by_chromosome",
    "heritability_bar_chart",
    # Metadata
    "metadata",
    "load_sample_metadata",
    "validate_metadata",
    "get_population_labels",
    "get_geographic_coordinates",
    # Workflow
    "execute_gwas_workflow",
    "run_gwas",
    "GWASWorkflowConfig",
    # Fine-mapping submodule
    "finemapping",
    "finemapping_credible_set",
    "susie_regression",
    "compute_bayes_factors",
    "colocalization",
    "conditional_analysis",
    "annotate_credible_set",
    "multi_trait_coloc",
    "eqtl_coloc",
    "compute_clpp",
    "regional_coloc",
    # Heritability submodule (LDSC, GREML, etc.)
    "heritability_submodule",
    "estimate_h2_ldsc",
    "partitioned_h2",
    "genetic_correlation",
    "haseman_elston_regression",
    "greml_simple",
    "compute_liability_h2",
]
