"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

The package exposes the historical convenience API lazily. Importing
``metainformant.gwas`` keeps subpackages and common functions discoverable
without eagerly importing the full workflow and visualization stack.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import TYPE_CHECKING, Any

_SUBPACKAGES = {
    "analysis",
    "data",
    "finemapping",
    "heritability",
    "reporting",
    "simulation",
    "validation",
    "visualization",
    "workflow",
}

_EXPORTS = {
    # Analysis
    "association_test_linear": "metainformant.gwas.analysis.association",
    "association_test_logistic": "metainformant.gwas.analysis.association",
    "association_test_mixed": "metainformant.gwas.analysis.mixed_model",
    "run_mixed_model_gwas": "metainformant.gwas.analysis.mixed_model",
    "bonferroni_correction": "metainformant.gwas.analysis.correction",
    "fdr_correction": "metainformant.gwas.analysis.correction",
    "estimate_heritability": "metainformant.gwas.analysis.heritability",
    "parse_vcf_full": "metainformant.gwas.analysis.quality",
    "compute_kinship_matrix": "metainformant.gwas.analysis.structure",
    "compute_pca": "metainformant.gwas.analysis.structure",
    "run_eqtl_analysis": "metainformant.gwas.analysis.eqtl",
    # Benchmarking
    "benchmark_subset_run": "metainformant.gwas.analysis.benchmarking",
    "extrapolate_full_genome_time": "metainformant.gwas.analysis.benchmarking",
    "scaling_model": "metainformant.gwas.analysis.benchmarking",
    "ComputeTimeEstimate": "metainformant.gwas.analysis.benchmarking",
    "StepTiming": "metainformant.gwas.analysis.benchmarking",
    # Data
    "download_reference_genome": "metainformant.gwas.data.download",
    "download_variant_database": "metainformant.gwas.data.download",
    "download_sra_run": "metainformant.gwas.data.download",
    "download_annotation": "metainformant.gwas.data.download",
    "download_variant_data": "metainformant.gwas.data.download",
    "search_sra_for_organism": "metainformant.gwas.data.download",
    "ExpressionLoader": "metainformant.gwas.data.expression",
    "get_geographic_coordinates": "metainformant.gwas.data.metadata",
    "get_population_labels": "metainformant.gwas.data.metadata",
    "load_sample_metadata": "metainformant.gwas.data.metadata",
    "merge_metadata_with_phenotypes": "metainformant.gwas.data.metadata",
    "validate_metadata": "metainformant.gwas.data.metadata",
    # Finemapping
    "eqtl_coloc": "metainformant.gwas.finemapping.colocalization",
    "multi_trait_coloc": "metainformant.gwas.finemapping.colocalization",
    "compute_clpp": "metainformant.gwas.finemapping.colocalization",
    "regional_coloc": "metainformant.gwas.finemapping.colocalization",
    "cis_eqtl_scan": "metainformant.gwas.finemapping.eqtl",
    "trans_eqtl_scan": "metainformant.gwas.finemapping.eqtl",
    "conditional_eqtl": "metainformant.gwas.finemapping.eqtl",
    "eqtl_effect_sizes": "metainformant.gwas.finemapping.eqtl",
    "eqtl_summary_stats": "metainformant.gwas.finemapping.eqtl",
    # Heritability
    "estimate_h2_ldsc": "metainformant.gwas.heritability.estimation",
    "partitioned_h2": "metainformant.gwas.heritability.estimation",
    "genetic_correlation": "metainformant.gwas.heritability.estimation",
    "greml_simple": "metainformant.gwas.heritability.estimation",
    "haseman_elston_regression": "metainformant.gwas.heritability.estimation",
    # Visualization
    "THEMES": "metainformant.gwas.visualization.config",
    "PlotStyle": "metainformant.gwas.visualization.config",
    "apply_style": "metainformant.gwas.visualization.config",
    "get_style": "metainformant.gwas.visualization.config",
    "style_from_config": "metainformant.gwas.visualization.config",
    "compute_ld_decay": "metainformant.gwas.visualization.genomic.ld",
    "ld_decay_plot": "metainformant.gwas.visualization.genomic.ld",
    "ld_heatmap_region": "metainformant.gwas.visualization.genomic.ld",
    "gwas_summary_panel": "metainformant.gwas.visualization.interactive.composite",
    "population_structure_panel": "metainformant.gwas.visualization.interactive.composite",
    "compute_credible_set": "metainformant.gwas.visualization.interactive.finemapping",
    "credible_set_plot": "metainformant.gwas.visualization.interactive.finemapping",
    "interactive_manhattan": "metainformant.gwas.visualization.interactive.interactive",
    "genotype_phenotype_boxplot": "metainformant.gwas.visualization.interactive.phenotype",
    "phenotype_correlation_matrix": "metainformant.gwas.visualization.interactive.phenotype",
    "phenotype_distribution": "metainformant.gwas.visualization.interactive.phenotype",
    "phenotype_pca_correlation": "metainformant.gwas.visualization.interactive.phenotype",
    "top_hits_genotype_phenotype": "metainformant.gwas.visualization.interactive.phenotype",
    "population_count_map": "metainformant.gwas.visualization.population.geography",
    "sample_map": "metainformant.gwas.visualization.population.geography",
    # Workflow
    "execute_gwas_workflow": "metainformant.gwas.workflow.workflow_execution",
    "run_gwas": "metainformant.gwas.workflow.workflow_execution",
    "run_multi_trait_gwas": "metainformant.gwas.workflow.workflow_execution",
}

__all__ = [
    "analysis",
    "data",
    "finemapping",
    "heritability",
    "reporting",
    "simulation",
    "validation",
    "visualization",
    "workflow",
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
    "benchmark_subset_run",
    "extrapolate_full_genome_time",
    "scaling_model",
    "ComputeTimeEstimate",
    "StepTiming",
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
    "eqtl_coloc",
    "multi_trait_coloc",
    "compute_clpp",
    "regional_coloc",
    "cis_eqtl_scan",
    "trans_eqtl_scan",
    "conditional_eqtl",
    "eqtl_effect_sizes",
    "eqtl_summary_stats",
    "estimate_h2_ldsc",
    "partitioned_h2",
    "genetic_correlation",
    "greml_simple",
    "haseman_elston_regression",
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
    "execute_gwas_workflow",
    "run_gwas",
    "run_multi_trait_gwas",
]


def __getattr__(name: str) -> Any:
    """Lazily load GWAS subpackages and convenience exports."""
    if name in _SUBPACKAGES:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    if name in _EXPORTS:
        module = importlib.import_module(_EXPORTS[name])
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Return public lazy exports for interactive help."""
    return sorted(set(globals()) | set(__all__))


if TYPE_CHECKING:
    from . import analysis, data, finemapping, heritability, reporting, simulation, validation, visualization, workflow
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
    from .data.download import (
        download_annotation,
        download_reference_genome,
        download_sra_run,
        download_variant_data,
        download_variant_database,
        search_sra_for_organism,
    )
    from .data.expression import ExpressionLoader
    from .data.metadata import (
        get_geographic_coordinates,
        get_population_labels,
        load_sample_metadata,
        merge_metadata_with_phenotypes,
        validate_metadata,
    )
    from .finemapping.colocalization import compute_clpp, eqtl_coloc, multi_trait_coloc, regional_coloc
    from .finemapping.eqtl import (
        cis_eqtl_scan,
        conditional_eqtl,
        eqtl_effect_sizes,
        eqtl_summary_stats,
        trans_eqtl_scan,
    )
    from .heritability.estimation import (
        estimate_h2_ldsc,
        genetic_correlation,
        greml_simple,
        haseman_elston_regression,
        partitioned_h2,
    )
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
    from .workflow.workflow_execution import execute_gwas_workflow, run_gwas, run_multi_trait_gwas
