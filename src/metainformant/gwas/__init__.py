"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

Provides a complete, end-to-end GWAS pipeline for population genomics:

* **QC** – VCF parsing, MAF/missing/HWE filtering, haplodiploidy detection
  (:mod:`~metainformant.gwas.analysis.quality`)
* **Structure** – VanRaden/IBS/Astle/Yang kinship and SVD-based PCA
  (:mod:`~metainformant.gwas.analysis.structure`)
* **Association** – linear, logistic, and EMMA mixed-model (Kang et al. 2008)
  (:mod:`~metainformant.gwas.analysis.association`,
  :mod:`~metainformant.gwas.analysis.mixed_model`)
* **Correction** – Bonferroni, Benjamini-Hochberg FDR, genomic control, q-values
  (:mod:`~metainformant.gwas.analysis.correction`)
* **Heritability** – GREML/REML h² estimation and per-chromosome partitioning
  (:mod:`~metainformant.gwas.heritability`)
* **Fine-mapping** – credible sets, eQTL co-localisation
  (:mod:`~metainformant.gwas.finemapping`)
* **Visualization** – Manhattan, Q-Q, LD heatmap, phenotype, geographic maps
  (:mod:`~metainformant.gwas.visualization`)
* **Workflow** – config-driven end-to-end orchestration, multi-trait GWAS
  (:mod:`~metainformant.gwas.workflow`)

Quick start::

    from metainformant.gwas import run_gwas
    results = run_gwas(
        vcf_path="data/variants.vcf",
        phenotype_path="data/phenotypes.tsv",
        config={"model": "mixed", "significance_threshold": 5e-8},
        output_dir="output/gwas/my_study",
    )
"""
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

# --- Visualization convenience exports ---
from .visualization.config import THEMES, PlotStyle, apply_style, get_style, style_from_config
from .visualization.genomic.ld import compute_ld_decay, ld_decay_plot, ld_heatmap_region
from .visualization.interactive.composite import gwas_summary_panel, population_structure_panel
from .visualization.interactive.finemapping import compute_credible_set, credible_set_plot
from .visualization.interactive.interactive import interactive_manhattan
from .visualization.interactive.phenotype import (
    genotype_phenotype_boxplot,
    phenotype_correlation_matrix,
    phenotype_distribution,
    phenotype_pca_correlation,
    top_hits_genotype_phenotype,
)
from .visualization.population.geography import population_count_map, sample_map

# --- Data/metadata additional exports ---
from .data.metadata import merge_metadata_with_phenotypes

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
    "merge_metadata_with_phenotypes",
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
    # Visualization
    "THEMES",
    "PlotStyle",
    "apply_style",
    "get_style",
    "style_from_config",
    "compute_ld_decay",
    "ld_decay_plot",
    "ld_heatmap_region",
    "gwas_summary_panel",
    "population_structure_panel",
    "compute_credible_set",
    "credible_set_plot",
    "interactive_manhattan",
    "genotype_phenotype_boxplot",
    "phenotype_correlation_matrix",
    "phenotype_distribution",
    "phenotype_pca_correlation",
    "top_hits_genotype_phenotype",
    "population_count_map",
    "sample_map",
    # Workflow
    "execute_gwas_workflow",
    "run_gwas",
    "run_multi_trait_gwas",
]
