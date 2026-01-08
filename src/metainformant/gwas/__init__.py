"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

This module provides comprehensive GWAS analysis capabilities, including
quality control, population structure analysis, association testing,
multiple testing correction, and publication-quality visualization.
"""

from __future__ import annotations

# Import subpackages (making them available at package level)
from . import data         # Configuration and data acquisition
from . import visualization # Visualization tools
from . import analysis     # Core analysis logic (depends on data)
from . import workflow     # Workflow orchestration (depends on all)

# Import modules from subpackages for backward compatibility and ease of use
from .analysis import (
    association,
    calling,
    correction,
    quality,
    structure,
)
from .data import (
    config,
    download,
    sra_download,
)
from .workflow import workflow as workflow_module

# Visualization modules
# Note: 'visualization' module was renamed to 'general' inside the subpackage
from .visualization import general as visualization_module
from .visualization import (
    visualization_comparison,
    visualization_effects,
    visualization_genome,
    visualization_population,
    visualization_regional,
    visualization_statistical,
    visualization_suite,
    visualization_variants,
)

# Direct imports of commonly used functions
from .analysis.association import association_test_linear, association_test_logistic
from .analysis.calling import check_bcftools_available
from .analysis.correction import bonferroni_correction, fdr_correction, genomic_control
from .analysis.quality import apply_qc_filters, parse_vcf_full, extract_variant_regions
from .analysis.structure import compute_kinship_matrix, compute_pca, estimate_population_structure

from .data.config import load_gwas_config
from .data.download import download_reference_genome, download_variant_data

# Visualization aliases
from .visualization.general import qq_plot, manhattan_plot, kinship_heatmap
from .visualization.visualization_variants import variant_density_plot
from .visualization.visualization_population import pca_scree_plot, admixture_plot
from .visualization.visualization_comparison import multi_trait_manhattan
from .visualization.visualization_effects import functional_enrichment_plot
from .visualization.visualization_statistical import power_plot
from .visualization.visualization_regional import recombination_rate_plot, effect_direction_plot, regional_plot, regional_ld_plot
from .visualization.visualization_genome import genome_wide_ld_heatmap

from .workflow.workflow import GWASWorkflowConfig, execute_gwas_workflow, run_gwas

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
    "estimate_population_structure",
    "association",
    "correction",
    
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
    
    # Workflow
    "execute_gwas_workflow",
    "run_gwas",
    "GWASWorkflowConfig",
]
