"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

This module provides comprehensive GWAS analysis capabilities, including
quality control, population structure analysis, association testing,
multiple testing correction, and publication-quality visualization.
"""

from __future__ import annotations

# Import all GWAS analysis submodules
from . import (
    association,
    calling,
    config,
    correction,
    download,
    quality,
    sra_download,
    structure,
    visualization,
    visualization_comparison,
    visualization_effects,
    visualization_genome,
    visualization_population,
    visualization_regional,
    visualization_statistical,
    visualization_suite,
    visualization_variants,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import visualization_comparison
except ImportError:
    visualization_comparison = None

try:
    from . import visualization_effects
except ImportError:
    visualization_effects = None

try:
    from . import visualization_genome
except ImportError:
    visualization_genome = None

try:
    from . import visualization_population
except ImportError:
    visualization_population = None

try:
    from . import visualization_regional
except ImportError:
    visualization_regional = None

try:
    from . import visualization_statistical
except ImportError:
    visualization_statistical = None

try:
    from . import visualization_suite
except ImportError:
    visualization_suite = None

try:
    from . import visualization_variants
except ImportError:
    visualization_variants = None

# Direct imports of commonly used functions
from .association import association_test_linear, association_test_logistic
from .calling import check_bcftools_available
from .config import load_gwas_config
from .correction import bonferroni_correction, fdr_correction, genomic_control
from .download import download_reference_genome, download_variant_data
from .quality import apply_qc_filters, parse_vcf_full, extract_variant_regions
from .structure import compute_kinship_matrix, compute_pca, estimate_population_structure
from .visualization import qq_plot, manhattan_plot, kinship_heatmap
from .visualization_variants import variant_density_plot
from .visualization_population import pca_scree_plot
from .visualization_comparison import multi_trait_manhattan
from .visualization_effects import functional_enrichment_plot
from .visualization_statistical import power_plot
from .visualization_regional import recombination_rate_plot, effect_direction_plot, regional_plot, regional_ld_plot
from .visualization_genome import genome_wide_ld_heatmap
from .visualization_population import admixture_plot
from .workflow import GWASWorkflowConfig, execute_gwas_workflow, run_gwas

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
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
    "workflow",

    # Data acquisition
    "download",
    "download_variant_data",
    "sra_download",
    "calling",

    # Visualization suite
    "visualization",
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

    # Configuration
    "config",
]



