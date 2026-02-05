"""Pharmacogenomics module for METAINFORMANT.

This module provides comprehensive pharmacogenomics analysis capabilities, including
star allele calling, diplotype determination, metabolizer phenotype prediction,
CPIC/PharmGKB guideline annotation, ACMG variant classification, drug-gene
interaction analysis, clinical report generation, and publication-quality visualization.

Subpackages:
    annotations: CPIC, PharmGKB, and FDA drug label annotation parsing
    alleles: Star allele calling, diplotype determination, phenotype prediction
    clinical: Pathogenicity classification, drug interactions, clinical reporting
    visualization: Metabolizer status, allele frequency, and drug response plots
"""

from __future__ import annotations

# Type checking imports
from typing import TYPE_CHECKING

# Import subpackages
from . import alleles, annotations, clinical, interaction, metabolism, visualization
from .alleles.diplotype import (
    Diplotype,
    calculate_activity_score,
    determine_diplotype,
    phased_diplotype,
    resolve_ambiguous_diplotypes,
)
from .alleles.phenotype import (
    MetabolizerPhenotype,
    classify_phenotype,
    get_phenotype_thresholds,
    population_phenotype_frequencies,
    predict_metabolizer_status,
)

# Data structure imports
# Allele imports
from .alleles.star_allele import (
    StarAllele,
    call_star_alleles,
    detect_novel_alleles,
    handle_cyp2d6_cnv,
    load_allele_definitions,
    match_allele_definition,
)

# Annotation imports
from .annotations.cpic import (
    get_dosing_recommendation,
    list_actionable_genes,
    load_cpic_guidelines,
    lookup_drug_gene,
    parse_cpic_allele_definitions,
)
from .annotations.drug_labels import (
    classify_label_type,
    extract_biomarker_info,
    parse_drug_label,
    search_labels_by_gene,
)
from .annotations.pharmgkb import (
    get_evidence_level,
    get_variant_annotations,
    parse_clinical_annotations,
    query_pharmgkb_annotations,
    search_drug_pathways,
)
from .clinical.drug_interaction import (
    DrugRecommendation,
    InteractionSeverity,
    analyze_drug_gene_interactions,
    calculate_interaction_severity,
    check_contraindications,
    polypharmacy_analysis,
    suggest_alternatives,
)

# Clinical imports
from .clinical.pathogenicity import (
    ACMGClassification,
    ACMGCriteria,
    aggregate_evidence,
    apply_acmg_criteria,
    check_gnomad_frequency,
    classify_variant_acmg,
    query_clinvar,
)
from .clinical.reporting import (
    add_disclaimer,
    export_report,
    format_recommendation,
    generate_clinical_report,
    generate_summary_table,
)

# Interaction imports
from .interaction.drug_interactions import (
    cyp_inhibition_prediction,
    default_interaction_database,
    polypharmacy_risk,
    predict_drug_interaction,
)

# Metabolism imports
from .metabolism.metabolizer_status import (
    classify_metabolizer,
    compute_activity_score,
    default_allele_function_table,
    dose_adjustment,
)
from .metabolism.metabolizer_status import predict_metabolizer_status as predict_metabolizer

# Visualization imports
from .visualization.plots import (
    plot_acmg_criteria,
    plot_activity_score_distribution,
    plot_allele_frequencies,
    plot_drug_response_heatmap,
    plot_metabolizer_status,
    plot_population_comparison,
)

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "annotations",
    "alleles",
    "clinical",
    "visualization",
    # Data structures
    "StarAllele",
    "Diplotype",
    "MetabolizerPhenotype",
    "ACMGClassification",
    "ACMGCriteria",
    "DrugRecommendation",
    "InteractionSeverity",
    # Annotations - CPIC
    "load_cpic_guidelines",
    "lookup_drug_gene",
    "get_dosing_recommendation",
    "list_actionable_genes",
    "parse_cpic_allele_definitions",
    # Annotations - PharmGKB
    "query_pharmgkb_annotations",
    "parse_clinical_annotations",
    "get_evidence_level",
    "search_drug_pathways",
    "get_variant_annotations",
    # Annotations - Drug labels
    "parse_drug_label",
    "extract_biomarker_info",
    "classify_label_type",
    "search_labels_by_gene",
    # Alleles - Star allele
    "call_star_alleles",
    "match_allele_definition",
    "load_allele_definitions",
    "detect_novel_alleles",
    "handle_cyp2d6_cnv",
    # Alleles - Diplotype
    "determine_diplotype",
    "calculate_activity_score",
    "resolve_ambiguous_diplotypes",
    "phased_diplotype",
    # Alleles - Phenotype
    "predict_metabolizer_status",
    "get_phenotype_thresholds",
    "classify_phenotype",
    "population_phenotype_frequencies",
    # Clinical - Pathogenicity
    "classify_variant_acmg",
    "apply_acmg_criteria",
    "query_clinvar",
    "aggregate_evidence",
    "check_gnomad_frequency",
    # Clinical - Drug interaction
    "analyze_drug_gene_interactions",
    "check_contraindications",
    "calculate_interaction_severity",
    "polypharmacy_analysis",
    "suggest_alternatives",
    # Clinical - Reporting
    "generate_clinical_report",
    "format_recommendation",
    "generate_summary_table",
    "export_report",
    "add_disclaimer",
    # Visualization
    "plot_metabolizer_status",
    "plot_allele_frequencies",
    "plot_activity_score_distribution",
    "plot_drug_response_heatmap",
    "plot_population_comparison",
    "plot_acmg_criteria",
    # Interaction
    "interaction",
    "predict_drug_interaction",
    "polypharmacy_risk",
    "cyp_inhibition_prediction",
    "default_interaction_database",
    # Metabolism
    "metabolism",
    "predict_metabolizer",
    "compute_activity_score",
    "classify_metabolizer",
    "dose_adjustment",
    "default_allele_function_table",
]
