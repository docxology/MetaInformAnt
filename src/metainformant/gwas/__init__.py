"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, data, finemapping, heritability, visualization, workflow

# --- Analysis convenience exports ---
from .analysis.association import association_test_linear, association_test_logistic
from .analysis.benchmarking import (
    ComputeTimeEstimate,
    StepTiming,
    benchmark_subset_run,
    extrapolate_full_genome_time,
    scaling_model,
)
from .analysis.correction import bonferroni_correction, fdr_correction
from .analysis.eqtl import run_eqtl_analysis
from .analysis.heritability import estimate_heritability
from .analysis.mixed_model import association_test_mixed, run_mixed_model_gwas
from .analysis.quality import parse_vcf_full
from .analysis.structure import compute_kinship_matrix, compute_pca

# --- Data convenience exports ---
from .data.download import (
    download_annotation,
    download_reference_genome,
    download_sra_run,
    download_variant_database,
    download_variant_data,
    search_sra_for_organism,
)
from .data.expression import ExpressionLoader
from .data.metadata import (
    get_geographic_coordinates,
    get_population_labels,
    load_sample_metadata,
    validate_metadata,
)

# --- Finemapping convenience exports ---
from .finemapping.colocalization import (
    eqtl_coloc,
    multi_trait_coloc,
    compute_clpp,
    regional_coloc,
)
from .finemapping.eqtl import (
    cis_eqtl_scan,
    trans_eqtl_scan,
    conditional_eqtl,
    eqtl_effect_sizes,
    eqtl_summary_stats,
)

# --- Heritability convenience exports ---
from .heritability.estimation import (
    estimate_h2_ldsc,
    genetic_correlation,
    greml_simple,
    haseman_elston_regression,
    partitioned_h2,
)

# --- Workflow convenience exports ---
from .workflow.workflow_execution import (
    execute_gwas_workflow,
    run_gwas,
    run_multi_trait_gwas,
)

__all__ = [
    # Sub-packages
    "analysis",
    "data",
    "finemapping",
    "heritability",
    "visualization",
    "workflow",
    # Analysis
    "association_test_linear",
    "association_test_logistic",
    "association_test_mixed",
    "run_mixed_model_gwas",
    "bonferroni_correction",
    "fdr_correction",
    "estimate_heritability",
    "parse_vcf_full",
    "compute_kinship_matrix",
    "compute_pca",
    "run_eqtl_analysis",
    # Benchmarking
    "benchmark_subset_run",
    "extrapolate_full_genome_time",
    "scaling_model",
    "ComputeTimeEstimate",
    "StepTiming",
    # Data
    "download_reference_genome",
    "download_variant_database",
    "download_sra_run",
    "download_annotation",
    "download_variant_data",
    "search_sra_for_organism",
    "ExpressionLoader",
    "get_geographic_coordinates",
    "get_population_labels",
    "load_sample_metadata",
    "validate_metadata",
    # Finemapping
    "eqtl_coloc",
    "multi_trait_coloc",
    "compute_clpp",
    "regional_coloc",
    "cis_eqtl_scan",
    "trans_eqtl_scan",
    "conditional_eqtl",
    "eqtl_effect_sizes",
    "eqtl_summary_stats",
    # Heritability
    "estimate_h2_ldsc",
    "partitioned_h2",
    "genetic_correlation",
    "greml_simple",
    "haseman_elston_regression",
    # Workflow
    "execute_gwas_workflow",
    "run_gwas",
    "run_multi_trait_gwas",
]
