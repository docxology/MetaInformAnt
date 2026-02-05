"""RNA analysis modules for expression analysis, QC, and validation.

This subpackage provides tools for analyzing RNA-seq data including:
- Expression analysis (normalization, differential expression)
- Quality control metrics and outlier detection
- Sample validation and pipeline status checking
- RNA-protein integration analysis
"""

from __future__ import annotations

from .protein_integration import (
    calculate_translation_efficiency,
    predict_protein_abundance_from_rna,
    ribosome_profiling_integration,
)
from .validation import (
    get_sample_pipeline_status,
    save_validation_report,
    validate_all_samples,
    validate_sample_end_to_end,
)

# Conditional imports for new modules (may not exist yet)
try:
    from .expression import (
        adjust_pvalues,
        compute_sample_distances,
        differential_expression,
        estimate_size_factors,
        filter_low_expression,
        get_highly_variable_genes,
        normalize_counts,
        pca_analysis,
        prepare_ma_data,
        prepare_volcano_data,
    )

    _HAS_EXPRESSION = True
except ImportError:
    _HAS_EXPRESSION = False

try:
    from .qc import (
        classify_expression_level,
        compute_correlation_matrix,
        compute_gene_metrics,
        compute_sample_metrics,
        compute_saturation_curve,
        detect_batch_effects,
        detect_gc_bias,
        detect_length_bias,
        detect_outlier_samples,
        estimate_library_complexity,
        generate_qc_report,
    )

    _HAS_QC = True
except ImportError:
    _HAS_QC = False

try:
    from .cross_species import (
        build_ortholog_map,
        compare_expression_across_species,
        compute_expression_conservation,
        compute_expression_divergence_matrix,
        cross_species_pca,
        identify_divergent_genes,
        map_expression_to_orthologs,
        phylogenetic_expression_profile,
    )

    _HAS_CROSS_SPECIES = True
except ImportError:
    _HAS_CROSS_SPECIES = False

__all__ = [
    # Validation
    "get_sample_pipeline_status",
    "validate_sample_end_to_end",
    "validate_all_samples",
    "save_validation_report",
    # Protein integration
    "calculate_translation_efficiency",
    "predict_protein_abundance_from_rna",
    "ribosome_profiling_integration",
]

# Add expression exports if available
if _HAS_EXPRESSION:
    __all__.extend(
        [
            "normalize_counts",
            "estimate_size_factors",
            "differential_expression",
            "adjust_pvalues",
            "pca_analysis",
            "compute_sample_distances",
            "filter_low_expression",
            "get_highly_variable_genes",
            "prepare_volcano_data",
            "prepare_ma_data",
        ]
    )

# Add QC exports if available
if _HAS_QC:
    __all__.extend(
        [
            "compute_sample_metrics",
            "detect_outlier_samples",
            "compute_gene_metrics",
            "classify_expression_level",
            "estimate_library_complexity",
            "compute_saturation_curve",
            "detect_batch_effects",
            "compute_correlation_matrix",
            "detect_gc_bias",
            "detect_length_bias",
            "generate_qc_report",
        ]
    )

# Add cross-species exports if available
if _HAS_CROSS_SPECIES:
    __all__.extend(
        [
            "build_ortholog_map",
            "map_expression_to_orthologs",
            "compare_expression_across_species",
            "compute_expression_conservation",
            "identify_divergent_genes",
            "compute_expression_divergence_matrix",
            "phylogenetic_expression_profile",
            "cross_species_pca",
        ]
    )
