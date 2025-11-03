"""GWAS (Genome-Wide Association Studies) functionality.

This module provides a comprehensive end-to-end GWAS pipeline including
variant data acquisition, quality control, population structure analysis,
association testing, and visualization.
"""

from .association import association_test_linear, association_test_logistic, run_gwas
from .calling import call_variants_bcftools, call_variants_gatk, merge_vcf_files
from .config import GWASWorkflowConfig, load_gwas_config
from .correction import bonferroni_correction, fdr_correction, genomic_control
from .download import download_reference_genome, download_variant_data, extract_variant_regions
from .quality import apply_qc_filters, parse_vcf_full
from .sra_download import (
    check_sra_tools_available,
    download_sra_project,
    download_sra_run,
    search_sra_for_organism,
)
from .structure import compute_pca, compute_kinship_matrix, estimate_population_structure
from .visualization import manhattan_plot, qq_plot, regional_plot

# Import all visualization modules
from .visualization_genome import (
    manhattan_plot as manhattan_plot_genome,
    circular_manhattan_plot,
    chromosome_ideogram,
    genome_wide_ld_heatmap,
)
from .visualization_statistical import (
    qq_plot as qq_plot_statistical,
    qq_plot_stratified,
    lambda_gc_plot,
    volcano_plot,
    power_plot,
)
from .visualization_regional import (
    regional_plot as regional_plot_detailed,
    regional_ld_plot,
    gene_annotation_plot,
    recombination_rate_plot,
)
from .visualization_population import (
    pca_plot,
    pca_scree_plot,
    kinship_heatmap,
    admixture_plot,
    population_tree,
)
from .visualization_variants import (
    maf_distribution,
    variant_density_plot,
    hwe_deviation_plot,
    missingness_plot,
    transition_transversion_plot,
)
from .visualization_effects import (
    effect_size_forest_plot,
    effect_direction_plot,
    functional_enrichment_plot,
    allelic_series_plot,
)
from .visualization_comparison import (
    miami_plot,
    multi_trait_manhattan,
    cross_cohort_forest,
    concordance_plot,
)
from .visualization_comprehensive import generate_all_plots

try:
    from .workflow import execute_gwas_workflow

    _HAS_WORKFLOW = True
except Exception:
    _HAS_WORKFLOW = False

__all__ = [
    # Configuration and workflow
    "GWASWorkflowConfig",
    "load_gwas_config",
    # Data acquisition
    "download_reference_genome",
    "download_variant_data",
    "extract_variant_regions",
    "check_sra_tools_available",
    "download_sra_run",
    "download_sra_project",
    "search_sra_for_organism",
    # Variant calling
    "call_variants_bcftools",
    "call_variants_gatk",
    "merge_vcf_files",
    # Quality control
    "parse_vcf_full",
    "apply_qc_filters",
    # Population structure
    "compute_pca",
    "compute_kinship_matrix",
    "estimate_population_structure",
    # Association testing
    "association_test_linear",
    "association_test_logistic",
    "run_gwas",
    # Multiple testing correction
    "bonferroni_correction",
    "fdr_correction",
    "genomic_control",
    # Basic visualization (legacy)
    "manhattan_plot",
    "qq_plot",
    "regional_plot",
    # Genome-wide visualization
    "manhattan_plot_genome",
    "circular_manhattan_plot",
    "chromosome_ideogram",
    "genome_wide_ld_heatmap",
    # Statistical visualization
    "qq_plot_statistical",
    "qq_plot_stratified",
    "lambda_gc_plot",
    "volcano_plot",
    "power_plot",
    # Regional visualization
    "regional_plot_detailed",
    "regional_ld_plot",
    "gene_annotation_plot",
    "recombination_rate_plot",
    # Population visualization
    "pca_plot",
    "pca_scree_plot",
    "kinship_heatmap",
    "admixture_plot",
    "population_tree",
    # Variant property visualization
    "maf_distribution",
    "variant_density_plot",
    "hwe_deviation_plot",
    "missingness_plot",
    "transition_transversion_plot",
    # Effect size visualization
    "effect_size_forest_plot",
    "effect_direction_plot",
    "functional_enrichment_plot",
    "allelic_series_plot",
    # Comparison visualization
    "miami_plot",
    "multi_trait_manhattan",
    "cross_cohort_forest",
    "concordance_plot",
    # Comprehensive visualization
    "generate_all_plots",
]

if _HAS_WORKFLOW:
    __all__.append("execute_gwas_workflow")

