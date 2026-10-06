"""RNA analysis modules for expression analysis, QC, and validation.

This subpackage provides tools for analyzing RNA-seq data including:
- Expression analysis (normalization, differential expression)
- Quality control metrics and outlier detection
- Sample validation and pipeline status checking
- RNA-protein integration analysis
- Cross-species profile conservation, TPM distribution summaries, and
  per-tissue completeness tables (descriptive only)"""

from __future__ import annotations

from metainformant.rna.analysis import native_artifact_validation

from . import (
    across_species_orchestrator,
    atlas_plots,
    conservation_profiles,
    cross_species,
    counting_statistics,
    expression_io,
    ortholog_diagnostics,
    ortholog_profiles,
    expression,
    expression_analysis,
    inferential_comparative,
    ortholog_mapping,
    protein_integration,
    phylogenetic_comparative,
    qc,
    qc_filtering,
    qc_metrics,
    statistics_contract,
    statistics_io,
    tissue_specificity,
    tissue_specificity_evolution,
    validation,
    within_species_orchestrator,
)

__all__ = [
    "across_species_orchestrator",
    "atlas_plots",
    "conservation_profiles",
    "cross_species",
    "counting_statistics",
    "native_artifact_validation",
    "expression_io",
    "ortholog_diagnostics",
    "ortholog_profiles",
    "expression",
    "expression_analysis",
    "inferential_comparative",
    "ortholog_mapping",
    "phylogenetic_comparative",
    "protein_integration",
    "qc",
    "qc_filtering",
    "qc_metrics",
    "statistics_contract",
    "statistics_io",
    "tissue_specificity",
    "tissue_specificity_evolution",
    "validation",
    "within_species_orchestrator",
]

