"""GWAS (Genome-Wide Association Studies) functionality.

This module provides a comprehensive end-to-end GWAS pipeline including
variant data acquisition, quality control, population structure analysis,
association testing, and visualization.
"""

from .association import association_test_linear, association_test_logistic, run_gwas
from .config import GWASWorkflowConfig, load_gwas_config
from .correction import bonferroni_correction, fdr_correction, genomic_control
from .download import download_reference_genome, download_variant_data, extract_variant_regions
from .quality import apply_qc_filters, parse_vcf_full
from .structure import compute_pca, compute_kinship_matrix, estimate_population_structure
from .visualization import manhattan_plot, qq_plot, regional_plot

try:
    from .workflow import execute_gwas_workflow

    _HAS_WORKFLOW = True
except Exception:
    _HAS_WORKFLOW = False

__all__ = [
    "GWASWorkflowConfig",
    "load_gwas_config",
    "download_reference_genome",
    "download_variant_data",
    "extract_variant_regions",
    "parse_vcf_full",
    "apply_qc_filters",
    "compute_pca",
    "compute_kinship_matrix",
    "estimate_population_structure",
    "association_test_linear",
    "association_test_logistic",
    "run_gwas",
    "bonferroni_correction",
    "fdr_correction",
    "genomic_control",
    "manhattan_plot",
    "qq_plot",
    "regional_plot",
]

if _HAS_WORKFLOW:
    __all__.append("execute_gwas_workflow")

