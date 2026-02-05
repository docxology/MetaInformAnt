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

# Import subpackages
from . import annotations
from . import alleles
from . import clinical
from . import interaction
from . import metabolism
from . import visualization

# Annotation imports
from .annotations.cpic import (
    load_cpic_guidelines,
    lookup_drug_gene,
    get_dosing_recommendation,
    list_actionable_genes,
    parse_cpic_allele_definitions,
)
from .annotations.pharmgkb import (
    query_pharmgkb_annotations,
    parse_clinical_annotations,
    get_evidence_level,
    search_drug_pathways,
    get_variant_annotations,
)
from .annotations.drug_labels import (
    parse_drug_label,
    extract_biomarker_info,
    classify_label_type,
    search_labels_by_gene,
)

# Allele imports
from .alleles.star_allele import (
    call_star_alleles,
    match_allele_definition,
    load_allele_definitions,
    detect_novel_alleles,
    handle_cyp2d6_cnv,
)
from .alleles.diplotype import (
    determine_diplotype,
    calculate_activity_score,
    resolve_ambiguous_diplotypes,
    phased_diplotype,
)
from .alleles.phenotype import (
    predict_metabolizer_status,
    get_phenotype_thresholds,
    classify_phenotype,
    population_phenotype_frequencies,
)

# Clinical imports
from .clinical.pathogenicity import (
    classify_variant_acmg,
    apply_acmg_criteria,
    query_clinvar,
    aggregate_evidence,
    check_gnomad_frequency,
)
from .clinical.drug_interaction import (
    analyze_drug_gene_interactions,
    check_contraindications,
    calculate_interaction_severity,
    polypharmacy_analysis,
    suggest_alternatives,
)
from .clinical.reporting import (
    generate_clinical_report,
    format_recommendation,
    generate_summary_table,
    export_report,
    add_disclaimer,
)

# Visualization imports
from .visualization.plots import (
    plot_metabolizer_status,
    plot_allele_frequencies,
    plot_activity_score_distribution,
    plot_drug_response_heatmap,
    plot_population_comparison,
    plot_acmg_criteria,
)

# Interaction imports
from .interaction.drug_interactions import (
    predict_drug_interaction,
    polypharmacy_risk,
    cyp_inhibition_prediction,
    default_interaction_database,
)

# Metabolism imports
from .metabolism.metabolizer_status import (
    predict_metabolizer_status as predict_metabolizer,
    compute_activity_score,
    classify_metabolizer,
    dose_adjustment,
    default_allele_function_table,
)

# Data structure imports
from .alleles.star_allele import StarAllele
from .alleles.diplotype import Diplotype
from .alleles.phenotype import MetabolizerPhenotype
from .clinical.pathogenicity import ACMGClassification, ACMGCriteria
from .clinical.drug_interaction import DrugRecommendation, InteractionSeverity

# Type checking imports
from typing import TYPE_CHECKING

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
