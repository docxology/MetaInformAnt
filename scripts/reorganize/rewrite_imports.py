#!/usr/bin/env python3
"""Rewrite __init__.py files and test imports to use canonical submodule paths.

Part 1: Replaces bloated __init__.py files with clean submodule-only imports.
Part 2: Rewrites test files to use canonical deep import paths instead of
         shortcut imports through __init__.py re-exports.

Usage:
    python scripts/reorganize/rewrite_imports.py              # dry-run (default)
    python scripts/reorganize/rewrite_imports.py --apply      # actually write files
    python scripts/reorganize/rewrite_imports.py --init-only  # only rewrite __init__.py
    python scripts/reorganize/rewrite_imports.py --tests-only # only rewrite test imports
"""
from __future__ import annotations

import argparse
import ast
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
SRC = ROOT / "src" / "metainformant"
TESTS = ROOT / "tests"

# ============================================================================
# Part 1: Clean __init__.py definitions
# ============================================================================

CLEAN_INITS: dict[str, str] = {
    "life_events": '''\
"""Life events and trajectory analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, core, models, survival, visualization, workflow

__all__ = ["analysis", "core", "models", "survival", "visualization", "workflow"]
''',
    "math": '''\
"""Mathematical biology and theoretical modeling module for METAINFORMANT."""
from __future__ import annotations

from . import bayesian, core, decision_theory, epidemiology, evolutionary_dynamics, population_genetics, quantitative_genetics

__all__ = ["bayesian", "core", "decision_theory", "epidemiology", "evolutionary_dynamics", "population_genetics", "quantitative_genetics"]
''',
    "information": '''\
"""Information theory analysis module for METAINFORMANT."""
from __future__ import annotations

from . import integration, metrics, network_info, workflow

__all__ = ["integration", "metrics", "network_info", "workflow"]
''',
    "networks": '''\
"""Network analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, config, interaction, regulatory, workflow

__all__ = ["analysis", "config", "interaction", "regulatory", "workflow"]
''',
    "dna": '''\
"""DNA sequence analysis and genomics module for METAINFORMANT."""
from __future__ import annotations

from . import alignment, annotation, expression, external, integration, io, phylogeny, population, sequence, variation

__all__ = ["alignment", "annotation", "expression", "external", "integration", "io", "phylogeny", "population", "sequence", "variation"]
''',
    "epigenome": '''\
"""Epigenome analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, assays, chromatin_state, peak_calling, visualization, workflow

__all__ = ["analysis", "assays", "chromatin_state", "peak_calling", "visualization", "workflow"]
''',
    "pharmacogenomics": '''\
"""Pharmacogenomics module for METAINFORMANT."""
from __future__ import annotations

from . import alleles, annotations, clinical, interaction, metabolism, visualization

__all__ = ["alleles", "annotations", "clinical", "interaction", "metabolism", "visualization"]
''',
    "longread": '''\
"""Long-read sequencing analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, assembly, io, methylation, phasing, quality, utils, visualization, workflow

__all__ = ["analysis", "assembly", "io", "methylation", "phasing", "quality", "utils", "visualization", "workflow"]
''',
    "spatial": '''\
"""Spatial transcriptomics module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, communication, deconvolution, integration, io, visualization

__all__ = ["analysis", "communication", "deconvolution", "integration", "io", "visualization"]
''',
    "rna": '''\
"""RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT."""
from __future__ import annotations

from . import amalgkit, analysis, core, deconvolution, engine, retrieval, splicing

__all__ = ["amalgkit", "analysis", "core", "deconvolution", "engine", "retrieval", "splicing"]
''',
    "gwas": '''\
"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, data, finemapping, heritability, visualization, workflow

__all__ = ["analysis", "data", "finemapping", "heritability", "visualization", "workflow"]
''',
    "visualization": '''\
"""Visualization and plotting utilities module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, config, dashboards, genomics, interactive_dashboards, plots

__all__ = ["analysis", "config", "dashboards", "genomics", "interactive_dashboards", "plots"]
''',
    "core": '''\
"""Core utilities for METAINFORMANT bioinformatics toolkit."""
from __future__ import annotations

# isort: skip_file
from . import data, engine, execution, io, ui, utils

# Optional dependencies
try:
    from .data import db
except ImportError:
    db = None

__all__ = ["data", "db", "engine", "execution", "io", "ui", "utils"]
''',
    # Additional modules that also have bloated re-exports
    "ecology": '''\
"""Ecology and community analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, phylogenetic, visualization

__all__ = ["analysis", "phylogenetic", "visualization"]
''',
    "simulation": '''\
"""Simulation and synthetic data generation module for METAINFORMANT."""
from __future__ import annotations

from . import benchmark, models, visualization, workflow

__all__ = ["benchmark", "models", "visualization", "workflow"]
''',
    "menu": '''\
"""Interactive menu and discovery system module for METAINFORMANT."""
from __future__ import annotations

from . import core, ui

__all__ = ["core", "ui"]
''',
    "quality": '''\
"""Quality control and metrics module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, io, reporting

__all__ = ["analysis", "io", "reporting"]
''',
    "singlecell": '''\
"""Single-cell analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, data, visualization

__all__ = ["analysis", "data", "visualization"]
''',
    "phenotype": '''\
"""Phenotype analysis module for METAINFORMANT."""
from __future__ import annotations

from . import analysis, behavior, chemical, data, electronic, gwas_integration, morphological, sonic, workflow

__all__ = ["analysis", "behavior", "chemical", "data", "electronic", "gwas_integration", "morphological", "sonic", "workflow"]
''',
    "multiomics": '''\
"""Multi-omics integration module for METAINFORMANT."""
from __future__ import annotations

from . import analysis

__all__ = ["analysis"]
''',
    "ontology": '''\
"""Ontology and gene set analysis module for METAINFORMANT."""
from __future__ import annotations

from . import core, pathway_enrichment, query, types, visualization

__all__ = ["core", "pathway_enrichment", "query", "types", "visualization"]
''',
    "protein": '''\
"""Protein analysis module for METAINFORMANT."""
from __future__ import annotations

from . import database, function, orchestration, sequence, structure

__all__ = ["database", "function", "orchestration", "sequence", "structure"]
''',
    "ml": '''\
"""Machine learning module for METAINFORMANT."""
from __future__ import annotations

from . import automl, evaluation, features, interpretability, models

__all__ = ["automl", "evaluation", "features", "interpretability", "models"]
''',
    "structural_variants": '''\
"""Structural variant detection and analysis module for METAINFORMANT."""
from __future__ import annotations

from . import annotation, detection, filtering, population, visualization

__all__ = ["annotation", "detection", "filtering", "population", "visualization"]
''',
    "metagenomics": '''\
"""Metagenomics analysis module for METAINFORMANT."""
from __future__ import annotations

from . import amplicon, comparative, diversity, functional, shotgun, visualization

__all__ = ["amplicon", "comparative", "diversity", "functional", "shotgun", "visualization"]
''',
}


# ============================================================================
# Part 2: Canonical import mapping (symbol -> canonical module path)
# ============================================================================
# Format: { "module_name": { "symbol_name": "canonical.module.path" } }
# The canonical path is used as: from <canonical_path> import <symbol>

CANONICAL_MAP: dict[str, dict[str, str]] = {
    "core": {
        "ProgressState": "metainformant.core.ui.tui",
        "SampleStage": "metainformant.core.engine.workflow_manager",
        "SampleState": "metainformant.core.engine.workflow_manager",
        "TerminalInterface": "metainformant.core.ui.tui",
        "Timer": "metainformant.core.utils.timing",
        "WorkflowManager": "metainformant.core.engine.workflow_manager",
        "atomic_replace": "metainformant.core.io.atomic",
        "atomic_write": "metainformant.core.io.atomic",
        "cache": "metainformant.core.io.cache",
        "compute_checksums_batch": "metainformant.core.io.checksums",
        "compute_md5": "metainformant.core.io.checksums",
        "compute_sha256": "metainformant.core.io.checksums",
        "cpu_count": "metainformant.core.execution.parallel",
        "create_sample_config": "metainformant.core.execution.workflow",
        "data": "metainformant.core.io.cache",
        "discover_configs": "metainformant.core.execution.discovery",
        "discover_functions": "metainformant.core.execution.discovery",
        "download": "metainformant.core.io.download",
        "download_and_process_data": "metainformant.core.execution.workflow",
        "download_file": "metainformant.core.io.io",
        "dump_json": "metainformant.core.io.io",
        "errors": "metainformant.core.execution.workflow",
        "load_json": "metainformant.core.io.io",
        "load_yaml": "metainformant.core.io.io",
        "parallel_batch": "metainformant.core.execution.parallel",
        "rate_limiter": "metainformant.core.utils.timing",
        "run_config_based_workflow": "metainformant.core.execution.workflow",
        "safe_write_bytes": "metainformant.core.io.atomic",
        "safe_write_text": "metainformant.core.io.atomic",
        "temp_directory": "metainformant.core.io.atomic",
        "text": "metainformant.core.utils.config",
        "thread_map": "metainformant.core.execution.parallel",
        "timed": "metainformant.core.utils.timing",
        "timeout_after": "metainformant.core.utils.timing",
        "validate_config_file": "metainformant.core.execution.workflow",
        "validate_not_empty": "metainformant.core.data.validation",
        "validate_not_none": "metainformant.core.data.validation",
        "validate_path_exists": "metainformant.core.data.validation",
        "validate_range": "metainformant.core.data.validation",
        "validate_schema": "metainformant.core.data.validation",
        "validate_type": "metainformant.core.data.validation",
        "verify_checksum": "metainformant.core.io.checksums",
        "verify_checksum_file": "metainformant.core.io.checksums",
        "write_checksum_file": "metainformant.core.io.checksums",
    },
    "dna": {
        "FastqRecord": "metainformant.dna.io.fastq",
        "codon": "metainformant.dna.annotation.functional",
        "consensus": "metainformant.dna.alignment.msa",
        "kmer": "metainformant.dna.phylogeny.tree_construction",
        "mutations": "metainformant.dna.variation.mutations",
        "sequence": "metainformant.dna.external.entrez",
        "sequences": "metainformant.dna.alignment.msa",
        "variants": "metainformant.dna.variation.variants",
    },
    "ecology": {
        "alpha_beta_gamma_diversity": "metainformant.ecology.analysis.community",
        "anosim": "metainformant.ecology.analysis.indicators",
        "beta_diversity": "metainformant.ecology.analysis.community",
        "build_simple_tree": "metainformant.ecology.phylogenetic.diversity",
        "calculate_biodiversity_indices": "metainformant.ecology.analysis.community",
        "calculate_diversity": "metainformant.ecology.analysis.community",
        "calculate_evenness": "metainformant.ecology.analysis.community",
        "calculate_single_diversity": "metainformant.ecology.analysis.community",
        "cca": "metainformant.ecology.analysis.ordination",
        "chao1_estimator": "metainformant.ecology.analysis.community",
        "cluster_communities": "metainformant.ecology.analysis.indicators",
        "community_metrics": "metainformant.ecology.analysis.community",
        "community_similarity_matrix": "metainformant.ecology.analysis.community",
        "community_weighted_mean": "metainformant.ecology.analysis.functional",
        "compare_sad_models": "metainformant.ecology.analysis.macroecology",
        "compute_unifrac": "metainformant.ecology.phylogenetic.diversity",
        "create_interactive_ecology_dashboard": "metainformant.ecology.visualization.visualization",
        "distance_decay": "metainformant.ecology.analysis.macroecology",
        "distance_matrix": "metainformant.ecology.analysis.ordination",
        "dominance_diversity_curve": "metainformant.ecology.analysis.community",
        "endemism_index": "metainformant.ecology.analysis.macroecology",
        "faiths_pd": "metainformant.ecology.phylogenetic.diversity",
        "fit_broken_stick": "metainformant.ecology.analysis.macroecology",
        "fit_geometric_series": "metainformant.ecology.analysis.macroecology",
        "fit_lognormal": "metainformant.ecology.analysis.macroecology",
        "fit_logseries": "metainformant.ecology.analysis.macroecology",
        "functional_beta_diversity": "metainformant.ecology.analysis.functional",
        "functional_dispersion": "metainformant.ecology.analysis.functional",
        "functional_divergence": "metainformant.ecology.analysis.functional",
        "functional_diversity_suite": "metainformant.ecology.analysis.functional",
        "functional_evenness": "metainformant.ecology.analysis.functional",
        "functional_redundancy": "metainformant.ecology.analysis.functional",
        "functional_richness": "metainformant.ecology.analysis.functional",
        "generate_ecology_report": "metainformant.ecology.analysis.community",
        "indval": "metainformant.ecology.analysis.indicators",
        "metabolic_scaling": "metainformant.ecology.analysis.macroecology",
        "multivariate_dispersion": "metainformant.ecology.analysis.indicators",
        "nestedness_temperature_calculator": "metainformant.ecology.analysis.community",
        "nmds": "metainformant.ecology.analysis.ordination",
        "nri_nti": "metainformant.ecology.phylogenetic.diversity",
        "occupancy_frequency": "metainformant.ecology.analysis.macroecology",
        "pcoa": "metainformant.ecology.analysis.ordination",
        "permanova": "metainformant.ecology.analysis.indicators",
        "phylogenetic_beta_diversity": "metainformant.ecology.phylogenetic.diversity",
        "phylogenetic_signal": "metainformant.ecology.phylogenetic.diversity",
        "pielou_evenness": "metainformant.ecology.analysis.community",
        "plot_beta_diversity_ordination": "metainformant.ecology.visualization.visualization",
        "plot_biodiversity_rarefaction": "metainformant.ecology.visualization.visualization",
        "plot_community_composition": "metainformant.ecology.visualization.visualization",
        "plot_diversity_accumulation_curve": "metainformant.ecology.visualization.visualization",
        "plot_diversity_indices_comparison": "metainformant.ecology.visualization.visualization",
        "plot_ecological_distance_heatmap": "metainformant.ecology.visualization.visualization",
        "plot_ecological_network": "metainformant.ecology.visualization.visualization",
        "plot_rank_abundance_curve_comparison": "metainformant.ecology.visualization.visualization",
        "plot_species_abundance_distribution": "metainformant.ecology.visualization.visualization",
        "procrustes": "metainformant.ecology.analysis.ordination",
        "rank_abundance_curve": "metainformant.ecology.analysis.community",
        "raos_quadratic_entropy": "metainformant.ecology.analysis.functional",
        "rarefaction_curve": "metainformant.ecology.analysis.community",
        "shannon_diversity": "metainformant.ecology.analysis.community",
        "simper": "metainformant.ecology.analysis.indicators",
        "simpson_diversity": "metainformant.ecology.analysis.community",
        "species_accumulation_curve": "metainformant.ecology.analysis.community",
        "species_area_logarithmic": "metainformant.ecology.analysis.macroecology",
        "species_area_power": "metainformant.ecology.analysis.macroecology",
        "species_area_relationship": "metainformant.ecology.analysis.community",
        "species_richness": "metainformant.ecology.analysis.community",
        "species_richness_simple": "metainformant.ecology.analysis.community",
        "taylors_power_law": "metainformant.ecology.analysis.macroecology",
        "trait_distance_matrix": "metainformant.ecology.analysis.functional",
    },
    "epigenome": {
        "assign_states": "metainformant.epigenome.chromatin_state.state_learning",
        "call_peaks_broad": "metainformant.epigenome.peak_calling.peak_detection",
        "call_peaks_simple": "metainformant.epigenome.peak_calling.peak_detection",
        "compare_chromatin_states": "metainformant.epigenome.chromatin_state.state_learning",
        "compute_beta_values": "metainformant.epigenome.assays.methylation",
        "compute_frip": "metainformant.epigenome.peak_calling.peak_detection",
        "compute_state_enrichment": "metainformant.epigenome.chromatin_state.state_learning",
        "differential_peaks": "metainformant.epigenome.peak_calling.peak_detection",
        "filter_peaks": "metainformant.epigenome.peak_calling.peak_detection",
        "interpret_states": "metainformant.epigenome.chromatin_state.state_learning",
        "learn_chromatin_states": "metainformant.epigenome.chromatin_state.state_learning",
        "load_cpg_table": "metainformant.epigenome.assays.methylation",
        "merge_peaks": "metainformant.epigenome.peak_calling.peak_detection",
        "read_bedgraph": "metainformant.epigenome.analysis.tracks",
        "segment_genome": "metainformant.epigenome.chromatin_state.state_learning",
        "summarize_beta_by_chromosome": "metainformant.epigenome.assays.methylation",
    },
    "gwas": {
        "GWASWorkflowConfig": "metainformant.gwas.workflow.workflow_config",
        "PlotStyle": "metainformant.gwas.visualization.config",
        "admixture_plot": "metainformant.gwas.visualization.population.population_admixture",
        "allele_frequency_map": "metainformant.gwas.visualization.population.geography",
        "annotate_credible_set": "metainformant.gwas.finemapping.credible_sets",
        "annotate_variants_with_genes": "metainformant.gwas.analysis.annotation",
        "annotation": "metainformant.gwas.analysis.annotation",
        "apply_qc_filters": "metainformant.gwas.analysis.quality",
        "apply_style": "metainformant.gwas.visualization.config",
        "association_test_linear": "metainformant.gwas.analysis.association",
        "association_test_logistic": "metainformant.gwas.analysis.association",
        "association_test_mixed": "metainformant.gwas.analysis.mixed_model",
        "bonferroni_correction": "metainformant.gwas.analysis.correction",
        "calling": "metainformant.gwas.data.config",
        "check_bcftools_available": "metainformant.gwas.analysis.calling",
        "classify_variant_location": "metainformant.gwas.analysis.annotation",
        "colocalization": "metainformant.gwas.finemapping.credible_sets",
        "compute_bayes_factors": "metainformant.gwas.finemapping.credible_sets",
        "compute_clpp": "metainformant.gwas.finemapping.colocalization",
        "compute_credible_set": "metainformant.gwas.finemapping.credible_sets",
        "compute_kinship_matrix": "metainformant.gwas.analysis.structure",
        "compute_ld_decay": "metainformant.gwas.visualization.genomic.ld",
        "compute_liability_h2": "metainformant.gwas.heritability.estimation",
        "compute_pca": "metainformant.gwas.analysis.structure",
        "conditional_analysis": "metainformant.gwas.finemapping.credible_sets",
        "config": "metainformant.gwas.analysis.structure",
        "correction": "metainformant.gwas.data.config",
        "create_results_summary": "metainformant.gwas.analysis.summary_stats",
        "credible_set_plot": "metainformant.gwas.visualization.interactive.finemapping",
        "data": "metainformant.gwas.analysis.quality",
        "download_reference_genome": "metainformant.gwas.data.download",
        "download_variant_data": "metainformant.gwas.data.download",
        "effect_direction_plot": "metainformant.gwas.visualization.genomic.regional",
        "eqtl_coloc": "metainformant.gwas.finemapping.colocalization",
        "estimate_h2_ldsc": "metainformant.gwas.heritability.estimation",
        "estimate_heritability": "metainformant.gwas.analysis.heritability",
        "estimate_population_structure": "metainformant.gwas.analysis.structure",
        "execute_gwas_workflow": "metainformant.gwas.workflow.workflow_execution",
        "extract_variant_regions": "metainformant.gwas.analysis.quality",
        "fdr_correction": "metainformant.gwas.analysis.correction",
        "functional_enrichment_plot": "metainformant.gwas.visualization.general",
        "genetic_correlation": "metainformant.gwas.heritability.estimation",
        "genome_wide_ld_heatmap": "metainformant.gwas.visualization.genomic.genome",
        "genomic_control": "metainformant.gwas.analysis.correction",
        "genotype_phenotype_boxplot": "metainformant.gwas.visualization.interactive.phenotype",
        "get_geographic_coordinates": "metainformant.gwas.data.metadata",
        "get_population_labels": "metainformant.gwas.data.metadata",
        "get_style": "metainformant.gwas.visualization.config",
        "greml_simple": "metainformant.gwas.heritability.estimation",
        "gwas_summary_panel": "metainformant.gwas.visualization.interactive.composite",
        "haseman_elston_regression": "metainformant.gwas.heritability.estimation",
        "heritability": "metainformant.gwas.analysis.mixed_model",
        "heritability_bar_chart": "metainformant.gwas.analysis.heritability",
        "interactive_manhattan": "metainformant.gwas.visualization.interactive.interactive",
        "interactive_pca": "metainformant.gwas.visualization.interactive.interactive",
        "interactive_volcano": "metainformant.gwas.visualization.interactive.interactive",
        "kinship_clustermap": "metainformant.gwas.visualization.population.population_admixture",
        "kinship_dendrogram": "metainformant.gwas.visualization.population.population_admixture",
        "kinship_heatmap": "metainformant.gwas.visualization.general",
        "ld_decay_plot": "metainformant.gwas.visualization.genomic.ld",
        "ld_heatmap_region": "metainformant.gwas.visualization.genomic.ld",
        "ld_prune": "metainformant.gwas.analysis.ld_pruning",
        "load_gwas_config": "metainformant.gwas.data.config",
        "load_sample_metadata": "metainformant.gwas.data.metadata",
        "manhattan_plot": "metainformant.gwas.visualization.general",
        "metadata": "metainformant.gwas.analysis.quality",
        "multi_trait_coloc": "metainformant.gwas.finemapping.colocalization",
        "multi_trait_manhattan": "metainformant.gwas.visualization.statistical.comparison",
        "normalize_chromosome_name": "metainformant.gwas.data.genome",
        "parse_vcf_full": "metainformant.gwas.analysis.quality",
        "partition_heritability_by_chromosome": "metainformant.gwas.analysis.heritability",
        "partitioned_h2": "metainformant.gwas.heritability.estimation",
        "pca_3d": "metainformant.gwas.visualization.population.population_pca",
        "pca_multi_panel": "metainformant.gwas.visualization.population.population_pca",
        "pca_scree_plot": "metainformant.gwas.visualization.population.population_pca",
        "phenotype_correlation_matrix": "metainformant.gwas.visualization.interactive.phenotype",
        "phenotype_distribution": "metainformant.gwas.visualization.interactive.phenotype",
        "phenotype_pca_correlation": "metainformant.gwas.visualization.interactive.phenotype",
        "pip_vs_ld_plot": "metainformant.gwas.visualization.interactive.finemapping",
        "population_count_map": "metainformant.gwas.visualization.population.geography",
        "population_structure_panel": "metainformant.gwas.visualization.interactive.composite",
        "power_plot": "metainformant.gwas.visualization.statistical.statistical",
        "qq_plot": "metainformant.gwas.visualization.general",
        "recombination_rate_plot": "metainformant.gwas.visualization.genomic.regional",
        "regional_coloc": "metainformant.gwas.finemapping.colocalization",
        "regional_ld_plot": "metainformant.gwas.visualization.genomic.regional",
        "regional_plot": "metainformant.gwas.visualization.general",
        "run_gwas": "metainformant.gwas.workflow.workflow_execution",
        "run_mixed_model_gwas": "metainformant.gwas.analysis.mixed_model",
        "sample_map": "metainformant.gwas.visualization.population.geography",
        "susie_regression": "metainformant.gwas.finemapping.credible_sets",
        "validate_metadata": "metainformant.gwas.data.metadata",
        "variant_density_plot": "metainformant.gwas.visualization.genomic.variants",
        "write_significant_hits": "metainformant.gwas.analysis.summary_stats",
        "write_summary_statistics": "metainformant.gwas.analysis.summary_stats",
    },
    "information": {
        "analyze_sequence_information": "metainformant.information.metrics.analysis.analysis",
        "annotation_specificity": "metainformant.information.metrics.advanced.semantic",
        "batch_entropy_analysis": "metainformant.information.workflow.workflows",
        "bias_correction": "metainformant.information.metrics.core.estimation",
        "binding_information": "metainformant.information.metrics.analysis.advanced_analysis",
        "channel_capacity": "metainformant.information.metrics.advanced.channel",
        "channel_mutual_information": "metainformant.information.metrics.advanced.channel",
        "co_information": "metainformant.information.metrics.advanced.decomposition",
        "compare_datasets": "metainformant.information.workflow.workflows",
        "compare_sequences_information": "metainformant.information.metrics.analysis.analysis",
        "conditional_entropy": "metainformant.information.metrics.core.syntactic",
        "conditional_entropy_continuous": "metainformant.information.metrics.core.continuous",
        "conditional_mutual_information": "metainformant.information.metrics.core.syntactic",
        "copula_entropy": "metainformant.information.metrics.core.continuous",
        "cross_entropy": "metainformant.information.metrics.core.syntactic",
        "differential_entropy": "metainformant.information.metrics.core.continuous",
        "dual_total_correlation": "metainformant.information.metrics.advanced.decomposition",
        "entropy_bootstrap_confidence": "metainformant.information.metrics.core.estimation",
        "entropy_confidence_interval": "metainformant.information.metrics.advanced.hypothesis",
        "entropy_estimation": "metainformant.information.metrics.core.continuous",
        "entropy_estimator": "metainformant.information.metrics.core.estimation",
        "entropy_power_inequality": "metainformant.information.metrics.advanced.information_projection",
        "entropy_rate_estimator": "metainformant.information.metrics.core.estimation",
        "entropy_rate_test": "metainformant.information.metrics.advanced.hypothesis",
        "exponential_family_entropy": "metainformant.information.metrics.advanced.fisher_rao",
        "fisher_information": "metainformant.information.metrics.analysis.advanced_analysis",
        "fisher_information_matrix": "metainformant.information.metrics.analysis.advanced_analysis",
        "fisher_rao_distance": "metainformant.information.metrics.advanced.fisher_rao",
        "granger_causality": "metainformant.information.network_info.information_flow",
        "hellinger_distance": "metainformant.information.metrics.advanced.fisher_rao",
        "independence_test": "metainformant.information.metrics.advanced.hypothesis",
        "information_bottleneck": "metainformant.information.metrics.advanced.channel",
        "information_coefficient": "metainformant.information.metrics.core.syntactic",
        "information_community_detection": "metainformant.information.integration.networks",
        "information_content": "metainformant.information.metrics.advanced.semantic",
        "information_content_from_annotations": "metainformant.information.metrics.advanced.semantic",
        "information_dimension": "metainformant.information.metrics.advanced.information_projection",
        "information_flow": "metainformant.information.integration.networks",
        "information_flow_network": "metainformant.information.metrics.core.continuous",
        "information_graph_distance": "metainformant.information.integration.networks",
        "information_profile": "metainformant.information.metrics.analysis.analysis",
        "information_projection": "metainformant.information.metrics.advanced.information_projection",
        "information_report": "metainformant.information.workflow.workflows",
        "information_signature": "metainformant.information.metrics.analysis.analysis",
        "information_significance_filter": "metainformant.information.metrics.advanced.hypothesis",
        "information_workflow": "metainformant.information.workflow.workflows",
        "interaction_information": "metainformant.information.metrics.analysis.advanced_analysis",
        "jensen_shannon_divergence": "metainformant.information.metrics.core.syntactic",
        "joint_entropy": "metainformant.information.metrics.core.syntactic",
        "kl_divergence": "metainformant.information.metrics.advanced.information_projection",
        "kl_divergence_continuous": "metainformant.information.metrics.core.continuous",
        "kl_divergence_estimator": "metainformant.information.metrics.core.estimation",
        "lautum_information": "metainformant.information.metrics.analysis.advanced_analysis",
        "metrics": "metainformant.information.integration.integration",
        "mi_permutation_test": "metainformant.information.metrics.advanced.hypothesis",
        "mutual_information": "metainformant.information.metrics.core.syntactic",
        "mutual_information_continuous": "metainformant.information.metrics.core.continuous",
        "mutual_information_estimator": "metainformant.information.metrics.core.estimation",
        "mutual_information_network": "metainformant.information.network_info.information_flow",
        "natural_gradient": "metainformant.information.metrics.advanced.fisher_rao",
        "network_entropy": "metainformant.information.integration.networks",
        "network_information_centrality": "metainformant.information.integration.networks",
        "network_motif_information": "metainformant.information.integration.networks",
        "noisy_channel_capacity": "metainformant.information.metrics.advanced.channel",
        "normalized_mutual_information": "metainformant.information.metrics.core.syntactic",
        "o_information": "metainformant.information.metrics.advanced.decomposition",
        "ontology_complexity": "metainformant.information.metrics.advanced.semantic",
        "panzeri_treves_bias_correction": "metainformant.information.metrics.core.estimation",
        "partial_information_decomposition": "metainformant.information.metrics.advanced.decomposition",
        "rate_distortion": "metainformant.information.metrics.advanced.channel",
        "redundant_information": "metainformant.information.metrics.advanced.decomposition",
        "relative_information_gain": "metainformant.information.metrics.analysis.advanced_analysis",
        "renyi_entropy": "metainformant.information.metrics.core.syntactic",
        "semantic_distance": "metainformant.information.metrics.advanced.semantic",
        "semantic_entropy": "metainformant.information.metrics.advanced.semantic",
        "semantic_similarity": "metainformant.information.metrics.advanced.semantic",
        "semantic_similarity_matrix": "metainformant.information.metrics.advanced.semantic",
        "shannon_entropy": "metainformant.information.metrics.advanced.channel",
        "shannon_entropy_from_counts": "metainformant.information.metrics.core.syntactic",
        "statistical_divergence": "metainformant.information.metrics.advanced.information_projection",
        "synergistic_information": "metainformant.information.metrics.advanced.decomposition",
        "term_redundancy": "metainformant.information.metrics.advanced.semantic",
        "term_specificity": "metainformant.information.metrics.advanced.semantic",
        "total_correlation": "metainformant.information.metrics.core.syntactic",
        "transfer_entropy": "metainformant.information.metrics.core.syntactic",
        "transfer_entropy_continuous": "metainformant.information.metrics.core.continuous",
        "tsallis_entropy": "metainformant.information.metrics.core.syntactic",
        "unique_information": "metainformant.information.metrics.advanced.decomposition",
        "variation_of_information": "metainformant.information.metrics.analysis.advanced_analysis",
    },
    "life_events": {
        "EnsemblePredictor": "metainformant.life_events.models.predictor",
        "Event": "metainformant.life_events.core.events",
        "EventDatabase": "metainformant.life_events.core.events",
        "EventSequence": "metainformant.life_events.core.events",
        "EventSequencePredictor": "metainformant.life_events.models.predictor",
        "GRUSequenceModel": "metainformant.life_events.models.sequence_models",
        "LSTMSequenceModel": "metainformant.life_events.models.sequence_models",
        "LifeEventsWorkflowConfig": "metainformant.life_events.core.config",
        "MultiTaskPredictor": "metainformant.life_events.models.statistical_models",
        "SurvivalPredictor": "metainformant.life_events.models.statistical_models",
        "add_temporal_noise": "metainformant.life_events.core.utils",
        "analyze_life_course": "metainformant.life_events.workflow.workflow",
        "attention_weights": "metainformant.life_events.analysis.interpretability",
        "biological_embedding": "metainformant.life_events.models.embeddings",
        "compare_populations": "metainformant.life_events.workflow.workflow",
        "competing_risks": "metainformant.life_events.survival.time_to_event",
        "config": "metainformant.life_events.core.config",
        "convert_sequences_to_tokens": "metainformant.life_events.core.utils",
        "cox_ph_model": "metainformant.life_events.survival.time_to_event",
        "domain_specific_embeddings": "metainformant.life_events.models.embeddings",
        "embeddings": "metainformant.life_events.analysis.interpretability",
        "event_importance": "metainformant.life_events.analysis.interpretability",
        "events": "metainformant.life_events.analysis.interpretability",
        "feature_attribution": "metainformant.life_events.analysis.interpretability",
        "generate_cohort_sequences": "metainformant.life_events.core.utils",
        "generate_event_chain": "metainformant.life_events.core.utils",
        "generate_synthetic_life_events": "metainformant.life_events.core.utils",
        "get_event_statistics": "metainformant.life_events.core.utils",
        "intervention_analysis": "metainformant.life_events.workflow.workflow",
        "kaplan_meier_estimator": "metainformant.life_events.survival.time_to_event",
        "learn_event_embeddings": "metainformant.life_events.analysis.interpretability",
        "load_life_events_config": "metainformant.life_events.core.config",
        "load_sequences_from_json": "metainformant.life_events.core.utils",
        "plot_attention_heatmap": "metainformant.life_events.visualization.statistical",
        "plot_domain_distribution": "metainformant.life_events.visualization.timeline",
        "plot_domain_timeline": "metainformant.life_events.visualization.timeline",
        "plot_embedding_clusters": "metainformant.life_events.visualization.statistical",
        "plot_event_cooccurrence": "metainformant.life_events.visualization.statistical",
        "plot_event_embeddings": "metainformant.life_events.visualization.statistical",
        "plot_event_frequency_heatmap": "metainformant.life_events.visualization.statistical",
        "plot_event_timeline": "metainformant.life_events.visualization.timeline",
        "plot_intervention_effects": "metainformant.life_events.visualization.statistical",
        "plot_outcome_distribution": "metainformant.life_events.visualization.statistical",
        "plot_population_comparison": "metainformant.life_events.visualization.statistical",
        "plot_prediction_accuracy": "metainformant.life_events.visualization.statistical",
        "plot_prediction_importance": "metainformant.life_events.visualization.statistical",
        "plot_sequence_length_distribution": "metainformant.life_events.visualization.statistical",
        "plot_sequence_similarity": "metainformant.life_events.visualization.statistical",
        "plot_temporal_density": "metainformant.life_events.visualization.timeline",
        "plot_temporal_patterns": "metainformant.life_events.visualization.timeline",
        "plot_transition_network": "metainformant.life_events.visualization.network",
        "recurrent_events": "metainformant.life_events.survival.time_to_event",
        "sequence_embeddings": "metainformant.life_events.core.utils",
        "temporal_patterns": "metainformant.life_events.analysis.interpretability",
        "time_varying_covariates": "metainformant.life_events.survival.time_to_event",
        "validate_sequence": "metainformant.life_events.core.utils",
    },
    "longread": {
        "BatchResult": "metainformant.longread.utils.batch",
        "LongReadOrchestrator": "metainformant.longread.workflow.orchestrator_core",
        "PipelineResult": "metainformant.longread.workflow.orchestrator_core",
        "PipelineStep": "metainformant.longread.workflow.orchestrator_core",
        "QCReport": "metainformant.longread.workflow.reporting",
        "RunSummary": "metainformant.longread.utils.summary",
        "aggregate_methylation": "metainformant.longread.analysis.modified_bases",
        "allele_specific_analysis": "metainformant.longread.phasing.haplotyping",
        "batch_compute_metrics": "metainformant.longread.utils.batch",
        "batch_filter_reads": "metainformant.longread.utils.batch",
        "build_haplotype_blocks": "metainformant.longread.analysis.phasing",
        "build_run_summary": "metainformant.longread.utils.summary",
        "calculate_alignment_stats": "metainformant.longread.io.bam",
        "calculate_consensus_quality": "metainformant.longread.assembly.consensus",
        "calculate_n50": "metainformant.longread.quality.metrics",
        "calculate_nx": "metainformant.longread.quality.metrics",
        "calculate_phase_block_stats": "metainformant.longread.analysis.phasing",
        "calculate_throughput": "metainformant.longread.quality.metrics",
        "call_5mc": "metainformant.longread.analysis.modified_bases",
        "call_6ma": "metainformant.longread.analysis.modified_bases",
        "call_methylation_from_signal": "metainformant.longread.methylation.calling",
        "compare_run_summaries": "metainformant.longread.utils.summary",
        "compute_methylation_stats": "metainformant.longread.methylation.calling",
        "compute_overlap_graph": "metainformant.longread.assembly.overlap",
        "compute_switch_errors": "metainformant.longread.phasing.haplotyping",
        "convert_pod5_to_fast5": "metainformant.longread.io.formats",
        "correct_with_short_reads": "metainformant.longread.assembly.hybrid",
        "detect_adapters": "metainformant.longread.quality.filtering",
        "detect_dmrs": "metainformant.longread.methylation.calling",
        "detect_insertions": "metainformant.longread.analysis.structural",
        "detect_inversions": "metainformant.longread.analysis.structural",
        "detect_methylation": "metainformant.longread.analysis.modified_bases",
        "detect_sv_from_long_reads": "metainformant.longread.analysis.structural",
        "differential_methylation": "metainformant.longread.analysis.modified_bases",
        "estimate_accuracy": "metainformant.longread.quality.metrics",
        "export_run_summary": "metainformant.longread.utils.summary",
        "extract_basecalls": "metainformant.longread.io.fast5",
        "extract_methylation_tags": "metainformant.longread.io.bam",
        "extract_signal": "metainformant.longread.io.fast5",
        "fast5_to_fastq": "metainformant.longread.io.formats",
        "filter_by_length": "metainformant.longread.quality.filtering",
        "filter_by_quality": "metainformant.longread.quality.filtering",
        "filter_contained_reads": "metainformant.longread.assembly.overlap",
        "find_overlaps": "metainformant.longread.assembly.overlap",
        "generate_assembly_summary": "metainformant.longread.utils.summary",
        "generate_consensus": "metainformant.longread.assembly.consensus",
        "generate_methylation_summary": "metainformant.longread.utils.summary",
        "generate_qc_summary": "metainformant.longread.utils.summary",
        "generate_sv_summary": "metainformant.longread.utils.summary",
        "get_read_metadata": "metainformant.longread.io.fast5",
        "get_supplementary_alignments": "metainformant.longread.io.bam",
        "haplotag_reads": "metainformant.longread.phasing.haplotyping",
        "hybrid_assemble": "metainformant.longread.assembly.hybrid",
        "methylation_pattern_analysis": "metainformant.longread.methylation.calling",
        "minimizer_sketch": "metainformant.longread.assembly.overlap",
        "multiple_sequence_alignment": "metainformant.longread.assembly.consensus",
        "phase_reads": "metainformant.longread.analysis.phasing",
        "phase_structural_variants": "metainformant.longread.analysis.structural",
        "plot_alignment_view": "metainformant.longread.visualization.plots",
        "plot_dotplot": "metainformant.longread.visualization.plots",
        "plot_methylation_track": "metainformant.longread.visualization.plots",
        "plot_phasing_blocks": "metainformant.longread.visualization.plots",
        "plot_quality_vs_length": "metainformant.longread.visualization.plots",
        "plot_read_length_histogram": "metainformant.longread.visualization.plots",
        "polish_consensus": "metainformant.longread.assembly.consensus",
        "process_batch": "metainformant.longread.utils.batch",
        "quality": "metainformant.longread.analysis.phasing",
        "quality_score_distribution": "metainformant.longread.quality.metrics",
        "read_fast5": "metainformant.longread.io.fast5",
        "read_length_stats": "metainformant.longread.quality.metrics",
        "read_long_read_bam": "metainformant.longread.io.bam",
        "scaffold_with_long_reads": "metainformant.longread.assembly.hybrid",
        "split_chimeric_reads": "metainformant.longread.quality.filtering",
        "tag_reads_by_haplotype": "metainformant.longread.analysis.phasing",
        "trim_adapters": "metainformant.longread.quality.filtering",
        "write_paf": "metainformant.longread.io.formats",
    },
    "math": {
        "abc_rejection": "metainformant.math.bayesian.inference",
        "basic_reproduction_number": "metainformant.math.epidemiology.models",
        "breeders_equation_response": "metainformant.math.population_genetics.selection",
        "compute_bayes_factor": "metainformant.math.bayesian.inference",
        "compute_dic": "metainformant.math.bayesian.inference",
        "compute_waic": "metainformant.math.bayesian.inference",
        "conjugate_beta_binomial": "metainformant.math.bayesian.inference",
        "conjugate_normal": "metainformant.math.bayesian.inference",
        "correlation": "metainformant.math.core.utilities",
        "correlation_coefficient": "metainformant.math.core.utilities",
        "covariance": "metainformant.math.core.utilities",
        "ddm_analytic_accuracy": "metainformant.math.decision_theory.ddm",
        "ddm_mean_decision_time": "metainformant.math.decision_theory.ddm",
        "delta_mean_trait": "metainformant.math.quantitative_genetics.price",
        "effective_reproduction_number": "metainformant.math.epidemiology.models",
        "effective_size_from_family_size_variance": "metainformant.math.population_genetics.statistics",
        "effective_size_sex_ratio": "metainformant.math.population_genetics.effective_size",
        "equilibrium_heterozygosity_infinite_alleles": "metainformant.math.population_genetics.statistics",
        "expectation": "metainformant.math.quantitative_genetics.price",
        "expected_pairwise_diversity": "metainformant.math.population_genetics.coalescent",
        "expected_r2_from_Ne_c": "metainformant.math.population_genetics.statistics",
        "expected_segregating_sites": "metainformant.math.population_genetics.coalescent",
        "expected_time_to_mrca": "metainformant.math.population_genetics.coalescent",
        "expected_total_branch_length": "metainformant.math.population_genetics.coalescent",
        "fisher_exact_test": "metainformant.math.core.utilities",
        "fixation_probability": "metainformant.math.population_genetics.statistics",
        "fst_from_allele_freqs": "metainformant.math.population_genetics.fst",
        "fst_from_heterozygosity": "metainformant.math.population_genetics.fst",
        "haldane_c_to_d": "metainformant.math.population_genetics.ld",
        "haldane_d_to_c": "metainformant.math.population_genetics.ld",
        "hardy_weinberg_genotype_freqs": "metainformant.math.population_genetics.core",
        "harmonic_mean_effective_size": "metainformant.math.population_genetics.effective_size",
        "herd_immunity_threshold": "metainformant.math.epidemiology.models",
        "heterozygosity_decay": "metainformant.math.population_genetics.core",
        "inbreeding_coefficient": "metainformant.math.population_genetics.core",
        "island_model_update": "metainformant.math.decision_theory.ddm",
        "jensen_shannon_divergence": "metainformant.math.core.utilities",
        "kin_selection_response": "metainformant.math.population_genetics.selection",
        "kosambi_c_to_d": "metainformant.math.population_genetics.ld",
        "kosambi_d_to_c": "metainformant.math.population_genetics.ld",
        "kurtosis": "metainformant.math.population_genetics.statistics",
        "lande_equation_response": "metainformant.math.quantitative_genetics.core",
        "ld_coefficients": "metainformant.math.population_genetics.ld",
        "ld_decay_r2": "metainformant.math.population_genetics.ld",
        "linear_regression": "metainformant.math.core.utilities",
        "logistic_map": "metainformant.math.evolutionary_dynamics.core",
        "lotka_volterra_step": "metainformant.math.evolutionary_dynamics.core",
        "metropolis_hastings": "metainformant.math.bayesian.inference",
        "mutation_selection_balance_dominant": "metainformant.math.population_genetics.core",
        "mutation_selection_balance_recessive": "metainformant.math.population_genetics.core",
        "mutation_update": "metainformant.math.population_genetics.selection",
        "narrow_sense_heritability": "metainformant.math.quantitative_genetics.core",
        "price_equation": "metainformant.math.quantitative_genetics.price",
        "r_squared": "metainformant.math.core.utilities",
        "realized_heritability": "metainformant.math.quantitative_genetics.core",
        "relative_fitness": "metainformant.math.population_genetics.selection",
        "replicator_derivative": "metainformant.math.evolutionary_dynamics.core",
        "replicator_step": "metainformant.math.evolutionary_dynamics.core",
        "seir_step": "metainformant.math.epidemiology.models",
        "selection_differential": "metainformant.math.population_genetics.selection",
        "selection_gradient": "metainformant.math.population_genetics.selection",
        "selection_intensity": "metainformant.math.population_genetics.selection",
        "selection_update": "metainformant.math.population_genetics.selection",
        "shannon_entropy": "metainformant.math.core.utilities",
        "sir_step": "metainformant.math.epidemiology.models",
        "sis_step": "metainformant.math.epidemiology.models",
        "skewness": "metainformant.math.population_genetics.statistics",
        "standard_deviation": "metainformant.math.population_genetics.statistics",
        "tajima_constants": "metainformant.math.population_genetics.coalescent",
        "tajimas_D": "metainformant.math.population_genetics.coalescent",
        "variance": "metainformant.math.population_genetics.fst",
        "watterson_theta": "metainformant.math.population_genetics.coalescent",
        "weighted_correlation": "metainformant.math.quantitative_genetics.price",
        "weighted_covariance": "metainformant.math.quantitative_genetics.price",
        "weighted_variance": "metainformant.math.quantitative_genetics.price",
    },
    "menu": {
        "Menu": "metainformant.menu.ui.navigation",
        "MenuHistory": "metainformant.menu.ui.navigation",
        "MenuItem": "metainformant.menu.ui.navigation",
        "MenuSystem": "metainformant.menu.ui.navigation",
        "ScriptInfo": "metainformant.menu.core.discovery",
        "categorize_script": "metainformant.menu.core.discovery",
        "clear_screen": "metainformant.menu.ui.display",
        "discover_scripts": "metainformant.menu.core.discovery",
        "execute_bash_script": "metainformant.menu.core.executor",
        "execute_python_script": "metainformant.menu.core.executor",
        "execute_script": "metainformant.menu.core.executor",
        "extract_script_metadata": "metainformant.menu.core.discovery",
        "format_breadcrumb": "metainformant.menu.ui.display",
        "format_menu": "metainformant.menu.ui.display",
        "generate_menu_from_scripts": "metainformant.menu.core.discovery",
        "get_choice": "metainformant.menu.ui.display",
        "get_current_menu": "metainformant.menu.ui.navigation",
        "go_back": "metainformant.menu.ui.navigation",
        "navigate_to_submenu": "metainformant.menu.ui.navigation",
        "prompt_for_args": "metainformant.menu.core.executor",
        "show_menu": "metainformant.menu.ui.display",
        "validate_script_executable": "metainformant.menu.core.executor",
    },
    "metagenomics": {
        "alpha_diversity": "metainformant.metagenomics.diversity.metrics",
        "annotate_genes": "metainformant.metagenomics.functional.annotation",
        "assemble_contigs": "metainformant.metagenomics.shotgun.assembly",
        "assess_bin_quality": "metainformant.metagenomics.shotgun.binning",
        "beta_diversity": "metainformant.metagenomics.diversity.metrics",
        "bin_contigs": "metainformant.metagenomics.shotgun.binning",
        "biomarker_discovery": "metainformant.metagenomics.comparative.differential_abundance",
        "build_kmer_index": "metainformant.metagenomics.shotgun.profiling",
        "build_taxonomy_tree": "metainformant.metagenomics.amplicon.taxonomy",
        "calculate_assembly_stats": "metainformant.metagenomics.shotgun.assembly",
        "calculate_confidence": "metainformant.metagenomics.amplicon.taxonomy",
        "calculate_identity": "metainformant.metagenomics.amplicon.otu_clustering",
        "calculate_pathway_completeness": "metainformant.metagenomics.functional.pathways",
        "calculate_relative_abundance": "metainformant.metagenomics.shotgun.profiling",
        "calculate_tetranucleotide_freq": "metainformant.metagenomics.shotgun.binning",
        "classify_gene_families": "metainformant.metagenomics.functional.annotation",
        "classify_taxonomy": "metainformant.metagenomics.amplicon.taxonomy",
        "clr_transform": "metainformant.metagenomics.comparative.differential_abundance",
        "cluster_otus": "metainformant.metagenomics.amplicon.otu_clustering",
        "compare_pathway_profiles": "metainformant.metagenomics.functional.pathways",
        "denoise_sequences": "metainformant.metagenomics.amplicon.asv_denoising",
        "differential_abundance": "metainformant.metagenomics.comparative.differential_abundance",
        "effect_size_analysis": "metainformant.metagenomics.comparative.differential_abundance",
        "estimate_error_rates": "metainformant.metagenomics.amplicon.asv_denoising",
        "filter_chimeras": "metainformant.metagenomics.amplicon.otu_clustering",
        "indicator_species": "metainformant.metagenomics.comparative.differential_abundance",
        "merge_paired_reads": "metainformant.metagenomics.amplicon.asv_denoising",
        "ordination": "metainformant.metagenomics.diversity.metrics",
        "permanova": "metainformant.metagenomics.diversity.metrics",
        "plot_alpha_diversity": "metainformant.metagenomics.visualization.plots",
        "plot_heatmap": "metainformant.metagenomics.visualization.plots",
        "plot_krona_chart": "metainformant.metagenomics.visualization.plots",
        "plot_ordination": "metainformant.metagenomics.visualization.plots",
        "plot_rarefaction_curves": "metainformant.metagenomics.visualization.plots",
        "plot_stacked_bar": "metainformant.metagenomics.visualization.plots",
        "predict_orfs": "metainformant.metagenomics.functional.annotation",
        "profile_community": "metainformant.metagenomics.shotgun.profiling",
        "rarefaction_curve": "metainformant.metagenomics.diversity.metrics",
        "rarefy": "metainformant.metagenomics.diversity.metrics",
        "reconstruct_pathways": "metainformant.metagenomics.functional.pathways",
        "refine_bins": "metainformant.metagenomics.shotgun.binning",
        "scaffold_contigs": "metainformant.metagenomics.shotgun.assembly",
    },
    "ml": {
        "auto_preprocess": "metainformant.ml.automl.optimization",
        "bayesian_optimization": "metainformant.ml.automl.optimization",
        "boruta_selection": "metainformant.ml.interpretability.feature_selection",
        "compute_attention_weights": "metainformant.ml.interpretability.explainers",
        "compute_lime_explanation": "metainformant.ml.interpretability.explainers",
        "compute_permutation_importance": "metainformant.ml.interpretability.explainers",
        "compute_shap_values_kernel": "metainformant.ml.interpretability.explainers",
        "feature_interaction": "metainformant.ml.interpretability.explainers",
        "features": "metainformant.ml.features.dimensionality",
        "grid_search": "metainformant.ml.automl.optimization",
        "model_selection": "metainformant.ml.automl.optimization",
        "mutual_information_selection": "metainformant.ml.interpretability.feature_selection",
        "partial_dependence": "metainformant.ml.interpretability.explainers",
        "random_search": "metainformant.ml.automl.optimization",
        "recursive_elimination": "metainformant.ml.interpretability.feature_selection",
        "stability_selection": "metainformant.ml.interpretability.feature_selection",
    },
    "networks": {
        "BiologicalNetwork": "metainformant.networks.analysis.graph_core",
        "CommunityDetectionConfig": "metainformant.networks.config",
        "GRNConfig": "metainformant.networks.config",
        "GeneRegulatoryNetwork": "metainformant.networks.interaction.regulatory_core",
        "NetworkConfig": "metainformant.networks.config",
        "NetworkWorkflow": "metainformant.networks.workflow",
        "NetworkWorkflowConfig": "metainformant.networks.config",
        "PPIConfig": "metainformant.networks.config",
        "PathwayEnrichmentConfig": "metainformant.networks.config",
        "PathwayNetwork": "metainformant.networks.analysis.pathway",
        "ProteinNetwork": "metainformant.networks.interaction.ppi",
        "analysis": "metainformant.networks.analysis.pathway",
        "build_pwm": "metainformant.networks.regulatory.motif_analysis",
        "centrality_measures": "metainformant.networks.analysis.graph_algorithms",
        "community": "metainformant.networks.config",
        "community_metrics": "metainformant.networks.analysis.community",
        "compute_network_motifs": "metainformant.networks.regulatory.grn_inference",
        "create_network": "metainformant.networks.analysis.graph_core",
        "detect_communities": "metainformant.networks.analysis.community",
        "evaluate_communities": "metainformant.networks.analysis.community",
        "export_network": "metainformant.networks.analysis.graph_algorithms",
        "filter_network": "metainformant.networks.analysis.graph_algorithms",
        "find_tf_binding_motifs": "metainformant.networks.regulatory.motif_analysis",
        "get_connected_components": "metainformant.networks.analysis.graph_algorithms",
        "get_network_summary": "metainformant.networks.analysis.graph_core",
        "graph": "metainformant.networks.analysis.graph_algorithms",
        "greedy_modularity_communities": "metainformant.networks.analysis.community",
        "import_network": "metainformant.networks.analysis.graph_algorithms",
        "infer_grn": "metainformant.networks.interaction.regulatory_analysis",
        "infer_grn_correlation": "metainformant.networks.regulatory.grn_inference",
        "infer_grn_mutual_info": "metainformant.networks.regulatory.grn_inference",
        "infer_grn_regression": "metainformant.networks.regulatory.grn_inference",
        "label_propagation_communities": "metainformant.networks.analysis.community",
        "load_network": "metainformant.networks.analysis.graph_core",
        "load_pathway_database": "metainformant.networks.analysis.pathway",
        "load_ppi_network": "metainformant.networks.interaction.ppi",
        "louvain_communities": "metainformant.networks.analysis.community",
        "network_enrichment_analysis": "metainformant.networks.analysis.pathway",
        "network_intersection": "metainformant.networks.analysis.graph_algorithms",
        "network_metrics": "metainformant.networks.analysis.graph_algorithms",
        "network_similarity": "metainformant.networks.analysis.graph_algorithms",
        "network_union": "metainformant.networks.analysis.graph_algorithms",
        "pathway": "metainformant.networks.config",
        "pathway_enrichment": "metainformant.networks.analysis.pathway",
        "pathway_enrichment_analysis": "metainformant.networks.analysis.pathway",
        "pathway_regulation_analysis": "metainformant.networks.interaction.regulatory_analysis",
        "ppi": "metainformant.networks.config",
        "ppi_network_analysis": "metainformant.networks.interaction.ppi",
        "predict_interactions": "metainformant.networks.interaction.ppi",
        "regulatory_motifs": "metainformant.networks.interaction.regulatory_analysis",
        "save_network": "metainformant.networks.analysis.graph_core",
        "scan_sequence_for_motifs": "metainformant.networks.regulatory.motif_analysis",
        "score_motif_match": "metainformant.networks.regulatory.motif_analysis",
        "score_regulators": "metainformant.networks.regulatory.grn_inference",
        "shortest_paths": "metainformant.networks.analysis.graph_algorithms",
        "validate_grn": "metainformant.networks.regulatory.grn_inference",
    },
    "ontology": {
        "compare_enrichments": "metainformant.ontology.pathway_enrichment.enrichment",
        "compute_enrichment_score": "metainformant.ontology.pathway_enrichment.enrichment",
        "gsea": "metainformant.ontology.pathway_enrichment.enrichment",
        "over_representation_analysis": "metainformant.ontology.pathway_enrichment.enrichment",
        "pathway_network": "metainformant.ontology.pathway_enrichment.enrichment",
    },
    "pharmacogenomics": {
        "ACMGClassification": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "ACMGCriteria": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "Diplotype": "metainformant.pharmacogenomics.alleles.diplotype",
        "DrugRecommendation": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "InteractionSeverity": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "MetabolizerPhenotype": "metainformant.pharmacogenomics.alleles.phenotype",
        "StarAllele": "metainformant.pharmacogenomics.alleles.star_allele",
        "add_disclaimer": "metainformant.pharmacogenomics.clinical.reporting",
        "aggregate_evidence": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "alleles": "metainformant.pharmacogenomics.alleles.phenotype",
        "analyze_drug_gene_interactions": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "apply_acmg_criteria": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "calculate_activity_score": "metainformant.pharmacogenomics.alleles.diplotype",
        "calculate_interaction_severity": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "call_star_alleles": "metainformant.pharmacogenomics.alleles.star_allele",
        "check_contraindications": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "check_gnomad_frequency": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "classify_label_type": "metainformant.pharmacogenomics.annotations.drug_labels",
        "classify_metabolizer": "metainformant.pharmacogenomics.metabolism.metabolizer_status",
        "classify_phenotype": "metainformant.pharmacogenomics.alleles.phenotype",
        "classify_variant_acmg": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "compute_activity_score": "metainformant.pharmacogenomics.metabolism.metabolizer_status",
        "cyp_inhibition_prediction": "metainformant.pharmacogenomics.interaction.drug_interactions",
        "default_allele_function_table": "metainformant.pharmacogenomics.metabolism.metabolizer_status",
        "default_interaction_database": "metainformant.pharmacogenomics.interaction.drug_interactions",
        "detect_novel_alleles": "metainformant.pharmacogenomics.alleles.star_allele",
        "determine_diplotype": "metainformant.pharmacogenomics.alleles.diplotype",
        "dose_adjustment": "metainformant.pharmacogenomics.metabolism.metabolizer_status",
        "export_report": "metainformant.pharmacogenomics.clinical.reporting",
        "extract_biomarker_info": "metainformant.pharmacogenomics.annotations.drug_labels",
        "format_recommendation": "metainformant.pharmacogenomics.clinical.reporting",
        "generate_clinical_report": "metainformant.pharmacogenomics.clinical.reporting",
        "generate_summary_table": "metainformant.pharmacogenomics.clinical.reporting",
        "get_dosing_recommendation": "metainformant.pharmacogenomics.annotations.cpic",
        "get_evidence_level": "metainformant.pharmacogenomics.annotations.pharmgkb",
        "get_phenotype_thresholds": "metainformant.pharmacogenomics.alleles.phenotype",
        "get_variant_annotations": "metainformant.pharmacogenomics.annotations.pharmgkb",
        "handle_cyp2d6_cnv": "metainformant.pharmacogenomics.alleles.star_allele",
        "list_actionable_genes": "metainformant.pharmacogenomics.annotations.cpic",
        "load_allele_definitions": "metainformant.pharmacogenomics.alleles.star_allele",
        "load_cpic_guidelines": "metainformant.pharmacogenomics.annotations.cpic",
        "lookup_drug_gene": "metainformant.pharmacogenomics.annotations.cpic",
        "match_allele_definition": "metainformant.pharmacogenomics.alleles.star_allele",
        "parse_clinical_annotations": "metainformant.pharmacogenomics.annotations.pharmgkb",
        "parse_cpic_allele_definitions": "metainformant.pharmacogenomics.annotations.cpic",
        "parse_drug_label": "metainformant.pharmacogenomics.annotations.drug_labels",
        "phased_diplotype": "metainformant.pharmacogenomics.alleles.diplotype",
        "plot_acmg_criteria": "metainformant.pharmacogenomics.visualization.plots",
        "plot_activity_score_distribution": "metainformant.pharmacogenomics.visualization.plots",
        "plot_allele_frequencies": "metainformant.pharmacogenomics.visualization.plots",
        "plot_drug_response_heatmap": "metainformant.pharmacogenomics.visualization.plots",
        "plot_metabolizer_status": "metainformant.pharmacogenomics.visualization.plots",
        "plot_population_comparison": "metainformant.pharmacogenomics.visualization.plots",
        "polypharmacy_analysis": "metainformant.pharmacogenomics.clinical.drug_interaction",
        "polypharmacy_risk": "metainformant.pharmacogenomics.interaction.drug_interactions",
        "population_phenotype_frequencies": "metainformant.pharmacogenomics.alleles.phenotype",
        "predict_drug_interaction": "metainformant.pharmacogenomics.interaction.drug_interactions",
        "predict_metabolizer_status": "metainformant.pharmacogenomics.alleles.phenotype",
        "query_clinvar": "metainformant.pharmacogenomics.clinical.pathogenicity",
        "query_pharmgkb_annotations": "metainformant.pharmacogenomics.annotations.pharmgkb",
        "resolve_ambiguous_diplotypes": "metainformant.pharmacogenomics.alleles.diplotype",
        "search_drug_pathways": "metainformant.pharmacogenomics.annotations.pharmgkb",
        "search_labels_by_gene": "metainformant.pharmacogenomics.annotations.drug_labels",
        "suggest_alternatives": "metainformant.pharmacogenomics.clinical.drug_interaction",
    },
    "phenotype": {
        "AcousticSignal": "metainformant.phenotype.sonic.signal",
        "BehaviorSequence": "metainformant.phenotype.behavior.sequence",
        "ChemicalProfile": "metainformant.phenotype.chemical.profile",
        "Compound": "metainformant.phenotype.chemical.compound",
        "Ethogram": "metainformant.phenotype.behavior.ethogram",
        "Measurement": "metainformant.phenotype.morphological.measurement",
        "MorphometricProfile": "metainformant.phenotype.morphological.profile",
        "PhenotypePipeline": "metainformant.phenotype.workflow.pipeline",
        "PipelineConfig": "metainformant.phenotype.workflow.pipeline",
        "PipelineResult": "metainformant.phenotype.workflow.pipeline",
        "TrackingPoint": "metainformant.phenotype.electronic.tracking",
        "Trajectory": "metainformant.phenotype.electronic.tracking",
        "behavior": "metainformant.phenotype.data.antwiki",
        "categorize_phenotypes": "metainformant.phenotype.gwas_integration.phewas",
        "data": "metainformant.phenotype.data.antwiki",
        "genetic_risk_score": "metainformant.phenotype.gwas_integration.phewas",
        "phenotype_correlation_matrix": "metainformant.phenotype.gwas_integration.phewas",
        "phenotype_heritability_screen": "metainformant.phenotype.gwas_integration.phewas",
        "run_phewas": "metainformant.phenotype.gwas_integration.phewas",
    },
    "protein": {
        "contacts": "metainformant.protein.structure.pdb",
        "domains": "metainformant.protein.database.interpro",
        "sequence": "metainformant.protein.orchestration",
        "sequences": "metainformant.protein.orchestration",
        "structure": "metainformant.protein.orchestration",
    },
    "quality": {
        "aggregate_sample_qc": "metainformant.quality.reporting.multiqc_integration",
        "check_qc_thresholds": "metainformant.quality.reporting.multiqc_integration",
        "default_qc_thresholds": "metainformant.quality.reporting.multiqc_integration",
        "generate_qc_report": "metainformant.quality.reporting.multiqc_integration",
        "metrics": "metainformant.quality.analysis.metrics",
        "qc_trend_analysis": "metainformant.quality.reporting.multiqc_integration",
    },
    "rna": {
        "AmalgkitRunLayout": "metainformant.rna.core.configs",
        "AmalgkitWorkflowConfig": "metainformant.rna.engine.workflow_core",
        "ProgressTracker": "metainformant.rna.engine.progress_tracker",
        "RNAPipelineConfig": "metainformant.rna.core.configs",
        "batch_deconvolve": "metainformant.rna.deconvolution.bulk_deconvolution",
        "build_isoform_graph": "metainformant.rna.splicing.isoforms",
        "build_signature_matrix": "metainformant.rna.deconvolution.bulk_deconvolution",
        "classify_splicing_events": "metainformant.rna.splicing.splice_analysis",
        "compare_isoform_usage": "metainformant.rna.splicing.isoforms",
        "compute_isoform_diversity": "metainformant.rna.splicing.isoforms",
        "compute_psi": "metainformant.rna.splicing.splice_analysis",
        "compute_sample_metrics": "metainformant.rna.analysis.qc_metrics",
        "compute_splice_site_strength": "metainformant.rna.splicing.splice_sites",
        "configs": "metainformant.rna.engine.orchestration",
        "deconvolve_nnls": "metainformant.rna.deconvolution.bulk_deconvolution",
        "deconvolve_svr": "metainformant.rna.deconvolution.bulk_deconvolution",
        "deps": "metainformant.rna.core.environment",
        "detect_splice_junctions": "metainformant.rna.splicing.splice_sites",
        "differential_expression": "metainformant.rna.analysis.expression_analysis",
        "differential_splicing": "metainformant.rna.splicing.splice_analysis",
        "enumerate_isoforms": "metainformant.rna.splicing.isoforms",
        "filter_low_expression": "metainformant.rna.analysis.expression_core",
        "find_novel_junctions": "metainformant.rna.splicing.splice_analysis",
        "generate_qc_report": "metainformant.rna.analysis.qc_filtering",
        "normalize_counts": "metainformant.rna.analysis.expression_core",
        "pca_analysis": "metainformant.rna.analysis.expression_analysis",
        "quantify_isoforms": "metainformant.rna.splicing.isoforms",
        "select_marker_genes": "metainformant.rna.deconvolution.bulk_deconvolution",
        "summarize_curate_tables": "metainformant.rna.engine.pipeline",
        "validate_deconvolution": "metainformant.rna.deconvolution.bulk_deconvolution",
        "validate_environment": "metainformant.rna.core.environment",
        "validation": "metainformant.rna.analysis.validation",
    },
    "simulation": {
        "AMINO_ACIDS": "metainformant.simulation.models.sequences",
        "Agent": "metainformant.simulation.models.agents",
        "DNA_BASES": "metainformant.simulation.models.sequences",
        "Ecosystem": "metainformant.simulation.models.agents",
        "GridAgent": "metainformant.simulation.models.agents",
        "GridWorld": "metainformant.simulation.models.agents",
        "RNA_BASES": "metainformant.simulation.models.sequences",
        "SimulationConfig": "metainformant.simulation.workflow.workflow",
        "add_agent": "metainformant.simulation.models.agents",
        "add_technical_noise": "metainformant.simulation.models.rna",
        "analyze_sequence_divergence": "metainformant.simulation.models.sequences",
        "animate_population_dynamics": "metainformant.simulation.visualization.visualization",
        "animate_sequence_evolution": "metainformant.simulation.visualization.visualization",
        "benchmark_suite": "metainformant.simulation.benchmark.generators",
        "calculate_biodiversity_metrics": "metainformant.simulation.models.agents",
        "calculate_sequence_similarity": "metainformant.simulation.models.sequences",
        "calibrate_simulation_parameters": "metainformant.simulation.workflow.workflow",
        "count_agents_by_type": "metainformant.simulation.models.agents",
        "create_ecosystem": "metainformant.simulation.models.agents",
        "create_interactive_simulation_dashboard": "metainformant.simulation.visualization.visualization",
        "create_simulation_config": "metainformant.simulation.workflow.workflow",
        "evaluate_benchmark": "metainformant.simulation.benchmark.generators",
        "evolve_sequence": "metainformant.simulation.models.sequences",
        "generate_benchmark_dataset": "metainformant.simulation.benchmark.generators",
        "generate_coding_sequence": "metainformant.simulation.models.sequences",
        "generate_genotype_matrix": "metainformant.simulation.models.popgen",
        "generate_linkage_disequilibrium_data": "metainformant.simulation.models.popgen",
        "generate_population_sequences": "metainformant.simulation.models.popgen",
        "generate_random_dna": "metainformant.simulation.models.sequences",
        "generate_random_protein": "metainformant.simulation.models.sequences",
        "generate_sequence_family": "metainformant.simulation.models.sequences",
        "generate_site_frequency_spectrum": "metainformant.simulation.models.popgen",
        "generate_synthetic_expression": "metainformant.simulation.benchmark.generators",
        "generate_synthetic_variants": "metainformant.simulation.benchmark.generators",
        "generate_two_populations": "metainformant.simulation.models.popgen",
        "get_population_dynamics": "metainformant.simulation.models.agents",
        "mutate_sequence": "metainformant.simulation.models.sequences",
        "plot_agent_based_model_results": "metainformant.simulation.visualization.visualization",
        "plot_evolutionary_simulation_summary": "metainformant.simulation.visualization.visualization",
        "plot_population_dynamics_simulation": "metainformant.simulation.visualization.visualization",
        "plot_rnaseq_simulation_results": "metainformant.simulation.visualization.visualization",
        "plot_sequence_evolution": "metainformant.simulation.visualization.visualization",
        "plot_simulation_parameter_sensitivity": "metainformant.simulation.visualization.visualization",
        "plot_simulation_validation_comparison": "metainformant.simulation.visualization.visualization",
        "remove_agent": "metainformant.simulation.models.agents",
        "reverse_transcribe_protein_to_dna": "metainformant.simulation.models.sequences",
        "run_benchmark_simulation": "metainformant.simulation.workflow.workflow",
        "run_simulation": "metainformant.simulation.models.agents",
        "run_simulation_workflow": "metainformant.simulation.workflow.workflow",
        "sequences": "metainformant.simulation.models.popgen",
        "simulate_admixture": "metainformant.simulation.models.popgen",
        "simulate_bottleneck_population": "metainformant.simulation.models.popgen",
        "simulate_bulk_rnaseq": "metainformant.simulation.models.rna",
        "simulate_competition": "metainformant.simulation.models.agents",
        "simulate_counts_negative_binomial": "metainformant.simulation.models.rna",
        "simulate_differential_expression": "metainformant.simulation.models.rna",
        "simulate_gene_duplication": "metainformant.simulation.models.sequences",
        "simulate_population_expansion": "metainformant.simulation.models.popgen",
        "simulate_predator_prey": "metainformant.simulation.models.agents",
        "simulate_rnaseq_counts": "metainformant.simulation.models.rna",
        "simulate_selection": "metainformant.simulation.models.popgen",
        "simulate_single_cell_rnaseq": "metainformant.simulation.models.rna",
        "simulate_spatial_expression": "metainformant.simulation.models.rna",
        "simulate_time_series_expression": "metainformant.simulation.models.rna",
        "simulation_step": "metainformant.simulation.models.agents",
        "translate_dna_to_protein": "metainformant.simulation.models.sequences",
        "validate_simulation_output": "metainformant.simulation.workflow.workflow",
    },
    "singlecell": {
        "data": "metainformant.singlecell.analysis.nonlinear_methods",
    },
    "spatial": {
        "CellBoundary": "metainformant.spatial.io.xenium",
        "CellMetadata": "metainformant.spatial.io.merfish",
        "DeconvolutionResult": "metainformant.spatial.analysis.deconvolution",
        "GearyCResult": "metainformant.spatial.analysis.autocorrelation",
        "GetisOrdResult": "metainformant.spatial.analysis.autocorrelation",
        "ImputationResult": "metainformant.spatial.integration.scrna_mapping",
        "InteractionResult": "metainformant.spatial.analysis.neighborhood",
        "LocalMoransResult": "metainformant.spatial.analysis.autocorrelation",
        "MERFISHDataset": "metainformant.spatial.io.merfish",
        "MappingResult": "metainformant.spatial.integration.scrna_mapping",
        "MoransIResult": "metainformant.spatial.analysis.autocorrelation",
        "NeighborhoodEnrichmentResult": "metainformant.spatial.analysis.neighborhood",
        "NicheResult": "metainformant.spatial.analysis.neighborhood",
        "RipleyKResult": "metainformant.spatial.analysis.neighborhood",
        "SpatialClusterResult": "metainformant.spatial.analysis.clustering",
        "SpatialDataset": "metainformant.spatial.io.visium",
        "TissuePosition": "metainformant.spatial.io.visium",
        "TranscriptSpot": "metainformant.spatial.io.merfish",
        "VariogramResult": "metainformant.spatial.analysis.autocorrelation",
        "XeniumDataset": "metainformant.spatial.io.xenium",
        "XeniumTranscript": "metainformant.spatial.io.xenium",
        "aggregate_to_cells": "metainformant.spatial.io.merfish",
        "anchor_based_transfer": "metainformant.spatial.integration.scrna_mapping",
        "build_communication_network": "metainformant.spatial.communication.cell_communication",
        "build_spatial_graph": "metainformant.spatial.analysis.clustering",
        "communication_pattern_analysis": "metainformant.spatial.communication.cell_communication",
        "compute_interaction_matrix": "metainformant.spatial.analysis.neighborhood",
        "compute_ligand_receptor_interactions": "metainformant.spatial.communication.cell_communication",
        "correlation_mapping": "metainformant.spatial.integration.scrna_mapping",
        "create_reference_profiles": "metainformant.spatial.analysis.deconvolution",
        "create_spatial_dataset": "metainformant.spatial.io.visium",
        "deconvolve_spots": "metainformant.spatial.analysis.deconvolution",
        "default_lr_database": "metainformant.spatial.communication.cell_communication",
        "enrichment_score": "metainformant.spatial.analysis.deconvolution",
        "estimate_cell_fractions": "metainformant.spatial.analysis.deconvolution",
        "filter_tissue_spots": "metainformant.spatial.io.visium",
        "gearys_c": "metainformant.spatial.analysis.autocorrelation",
        "getis_ord_g": "metainformant.spatial.analysis.autocorrelation",
        "impute_spatial_genes": "metainformant.spatial.integration.scrna_mapping",
        "leiden_clustering": "metainformant.spatial.analysis.clustering",
        "ligand_receptor_spatial": "metainformant.spatial.analysis.neighborhood",
        "load_cell_boundaries": "metainformant.spatial.io.xenium",
        "load_merfish": "metainformant.spatial.io.merfish",
        "load_transcript_spots": "metainformant.spatial.io.merfish",
        "load_visium": "metainformant.spatial.io.visium",
        "load_xenium": "metainformant.spatial.io.xenium",
        "local_morans_i": "metainformant.spatial.analysis.autocorrelation",
        "louvain_clustering": "metainformant.spatial.analysis.clustering",
        "map_scrna_to_spatial": "metainformant.spatial.integration.scrna_mapping",
        "morans_i": "metainformant.spatial.analysis.autocorrelation",
        "neighborhood_enrichment": "metainformant.spatial.analysis.neighborhood",
        "niche_detection": "metainformant.spatial.analysis.neighborhood",
        "niche_identification": "metainformant.spatial.deconvolution.spatial_deconvolution",
        "nnls_deconvolution": "metainformant.spatial.analysis.deconvolution",
        "parse_cell_metadata": "metainformant.spatial.io.merfish",
        "plot_cell_type_map": "metainformant.spatial.visualization.plots",
        "plot_deconvolution_pie": "metainformant.spatial.visualization.plots",
        "plot_gene_expression_map": "metainformant.spatial.visualization.plots",
        "plot_interaction_heatmap": "metainformant.spatial.visualization.plots",
        "plot_neighborhood_graph": "metainformant.spatial.visualization.plots",
        "plot_spatial_autocorrelation": "metainformant.spatial.visualization.plots",
        "plot_spatial_scatter": "metainformant.spatial.visualization.plots",
        "plot_tissue_overlay": "metainformant.spatial.visualization.plots",
        "read_cell_features": "metainformant.spatial.io.xenium",
        "read_spatial_image": "metainformant.spatial.io.visium",
        "read_tissue_positions": "metainformant.spatial.io.visium",
        "read_transcripts": "metainformant.spatial.io.xenium",
        "ripley_k": "metainformant.spatial.analysis.neighborhood",
        "spatial_cell_type_mapping": "metainformant.spatial.deconvolution.spatial_deconvolution",
        "spatial_cluster": "metainformant.spatial.analysis.clustering",
        "spatial_domains": "metainformant.spatial.analysis.clustering",
        "spatial_interaction_score": "metainformant.spatial.communication.cell_communication",
        "spatial_variogram": "metainformant.spatial.analysis.autocorrelation",
        "spatial_weights_matrix": "metainformant.spatial.analysis.autocorrelation",
        "validate_deconvolution": "metainformant.spatial.deconvolution.spatial_deconvolution",
    },
    "structural_variants": {
        "annotate_gene_overlap": "metainformant.structural_variants.annotation.overlap",
        "annotate_regulatory_overlap": "metainformant.structural_variants.annotation.overlap",
        "apply_blacklist": "metainformant.structural_variants.filtering.quality_filter",
        "assess_dosage_sensitivity": "metainformant.structural_variants.annotation.functional_impact",
        "calculate_breakpoint_confidence": "metainformant.structural_variants.detection.breakpoints",
        "calculate_log2_ratio": "metainformant.structural_variants.detection.cnv",
        "calculate_overlap_fraction": "metainformant.structural_variants.annotation.overlap",
        "calculate_reciprocal_overlap": "metainformant.structural_variants.filtering.merge",
        "call_cnv_states": "metainformant.structural_variants.detection.cnv",
        "call_structural_variants": "metainformant.structural_variants.detection.sv_calling",
        "classify_sv_type": "metainformant.structural_variants.detection.sv_calling",
        "cluster_breakpoints": "metainformant.structural_variants.detection.breakpoints",
        "deduplicate_variants": "metainformant.structural_variants.filtering.merge",
        "detect_cnv_from_depth": "metainformant.structural_variants.detection.cnv",
        "detect_discordant_pairs": "metainformant.structural_variants.detection.sv_calling",
        "detect_microhomology": "metainformant.structural_variants.detection.breakpoints",
        "detect_split_reads": "metainformant.structural_variants.detection.sv_calling",
        "filter_by_frequency": "metainformant.structural_variants.filtering.quality_filter",
        "filter_by_quality": "metainformant.structural_variants.filtering.quality_filter",
        "filter_by_size": "metainformant.structural_variants.filtering.quality_filter",
        "find_nearest_gene": "metainformant.structural_variants.annotation.overlap",
        "genotype_sv": "metainformant.structural_variants.detection.sv_calling",
        "genotype_sv_population": "metainformant.structural_variants.population.sv_population",
        "merge_adjacent_segments": "metainformant.structural_variants.detection.cnv",
        "merge_callsets": "metainformant.structural_variants.filtering.merge",
        "merge_sv_callsets": "metainformant.structural_variants.population.sv_population",
        "plot_breakpoint_detail": "metainformant.structural_variants.visualization.plots",
        "plot_circos": "metainformant.structural_variants.visualization.plots",
        "plot_cnv_profile": "metainformant.structural_variants.visualization.plots",
        "plot_coverage_track": "metainformant.structural_variants.visualization.plots",
        "plot_sv_size_distribution": "metainformant.structural_variants.visualization.plots",
        "plot_sv_type_summary": "metainformant.structural_variants.visualization.plots",
        "predict_functional_impact": "metainformant.structural_variants.annotation.functional_impact",
        "predict_tad_disruption": "metainformant.structural_variants.annotation.functional_impact",
        "refine_breakpoints": "metainformant.structural_variants.detection.breakpoints",
        "score_pathogenicity": "metainformant.structural_variants.annotation.functional_impact",
        "segment_coverage": "metainformant.structural_variants.detection.cnv",
        "survivor_merge": "metainformant.structural_variants.filtering.merge",
        "sv_allele_frequency": "metainformant.structural_variants.population.sv_population",
        "sv_association_test": "metainformant.structural_variants.population.sv_population",
        "sv_ld_analysis": "metainformant.structural_variants.population.sv_population",
        "sv_population_structure": "metainformant.structural_variants.population.sv_population",
    },
    "visualization": {
        "alternating_pair": "metainformant.visualization.config.palettes",
        "animate_clustering": "metainformant.visualization.plots.animations",
        "animate_evolution": "metainformant.visualization.plots.animations",
        "animate_network": "metainformant.visualization.plots.animations",
        "animate_time_series": "metainformant.visualization.plots.animations",
        "animate_trajectory": "metainformant.visualization.plots.animations",
        "apply_theme": "metainformant.visualization.config.themes",
        "area_plot": "metainformant.visualization.plots.basic",
        "bar_plot": "metainformant.visualization.plots.basic",
        "categorical": "metainformant.visualization.config.palettes",
        "chromosome_palette": "metainformant.visualization.config.palettes",
        "circular_tree_plot": "metainformant.visualization.genomics.trees",
        "config": "metainformant.visualization.interactive_dashboards.dashboards",
        "correlation_heatmap": "metainformant.visualization.analysis.statistical",
        "create_dashboard": "metainformant.visualization.interactive_dashboards.dashboards",
        "create_genome_browser_track": "metainformant.visualization.interactive_dashboards.dashboards",
        "create_interactive_heatmap": "metainformant.visualization.interactive_dashboards.dashboards",
        "create_interactive_scatter": "metainformant.visualization.interactive_dashboards.dashboards",
        "create_interactive_volcano": "metainformant.visualization.interactive_dashboards.dashboards",
        "export_to_html": "metainformant.visualization.interactive_dashboards.dashboards",
        "expression_gradient": "metainformant.visualization.config.palettes",
        "expression_heatmap": "metainformant.visualization.plots.general",
        "genomic_overview": "metainformant.visualization.dashboards.composite",
        "get_theme": "metainformant.visualization.config.themes",
        "heatmap": "metainformant.visualization.plots.basic",
        "heatmap_cmap": "metainformant.visualization.config.palettes",
        "interactive_heatmap": "metainformant.visualization.dashboards.interactive",
        "interactive_manhattan": "metainformant.visualization.dashboards.interactive",
        "interactive_scatter": "metainformant.visualization.dashboards.interactive",
        "interactive_volcano": "metainformant.visualization.dashboards.interactive",
        "lineplot": "metainformant.visualization.plots.basic",
        "list_themes": "metainformant.visualization.config.themes",
        "manhattan_plot": "metainformant.visualization.genomics.genomics",
        "multi_panel": "metainformant.visualization.dashboards.composite",
        "pca_plot": "metainformant.visualization.plots.general",
        "pie_chart": "metainformant.visualization.plots.basic",
        "plot_3d_scatter": "metainformant.visualization.plots.multidim",
        "plot_alluvial_diagram": "metainformant.visualization.plots.specialized",
        "plot_chord_diagram": "metainformant.visualization.plots.specialized",
        "plot_circular_barplot": "metainformant.visualization.plots.specialized",
        "plot_network_circular_layout": "metainformant.visualization.plots.specialized",
        "plot_pairwise_relationships": "metainformant.visualization.plots.multidim",
        "plot_parallel_coordinates": "metainformant.visualization.plots.multidim",
        "plot_phylo_tree": "metainformant.visualization.genomics.trees",
        "plot_radar_chart": "metainformant.visualization.plots.multidim",
        "plot_sankey_diagram": "metainformant.visualization.plots.specialized",
        "plot_upset_plot": "metainformant.visualization.plots.specialized",
        "plot_venn_diagram": "metainformant.visualization.plots.specialized",
        "qc_summary": "metainformant.visualization.dashboards.composite",
        "qq_plot": "metainformant.visualization.analysis.statistical",
        "register_theme": "metainformant.visualization.config.themes",
        "reset_theme": "metainformant.visualization.config.themes",
        "scatter_plot": "metainformant.visualization.plots.basic",
        "significance_color": "metainformant.visualization.config.palettes",
        "significance_palette": "metainformant.visualization.config.palettes",
        "step_plot": "metainformant.visualization.plots.basic",
        "theme": "metainformant.visualization.config.themes",
        "unrooted_tree_plot": "metainformant.visualization.genomics.trees",
        "volcano_plot": "metainformant.visualization.genomics.genomics",
    },
}


# ============================================================================
# Part 1 implementation: rewrite __init__.py files
# ============================================================================


def rewrite_init_files(*, apply: bool = False) -> int:
    """Rewrite all __init__.py files with clean submodule-only imports.

    Returns the number of files that would be (or were) changed.
    """
    changed = 0

    for module_name, clean_content in sorted(CLEAN_INITS.items()):
        init_path = SRC / module_name / "__init__.py"

        if not init_path.exists():
            print(f"  [SKIP] {module_name}/__init__.py does not exist")
            continue

        current = init_path.read_text(encoding="utf-8")

        if current.strip() == clean_content.strip():
            print(f"  [OK]   {module_name}/__init__.py already clean")
            continue

        changed += 1
        rel = init_path.relative_to(ROOT)

        if apply:
            init_path.write_text(clean_content, encoding="utf-8")
            print(f"  [WRITE] {rel}")
        else:
            print(f"  [WOULD WRITE] {rel}")

    return changed


# ============================================================================
# Part 2 implementation: rewrite test imports
# ============================================================================

# Regex to match: from metainformant.<module> import <names>
# This handles single-line and multi-line (parenthesized) imports.
_IMPORT_RE = re.compile(
    r"^(\s*)from\s+metainformant\.(\w+)\s+import\s+"
    r"(?:\(([^)]+)\)|(.+))$",
    re.MULTILINE,
)

# Regex to match multi-line import blocks: from metainformant.<mod> import (
_MULTILINE_START_RE = re.compile(
    r"^(\s*)from\s+metainformant\.(\w+)\s+import\s+\(\s*$"
)
_MULTILINE_END_RE = re.compile(r"^(\s*)\)\s*$")


def parse_import_block(lines: list[str], start_idx: int) -> tuple[str, str, list[str], int] | None:
    """Parse an import block starting at start_idx.

    Returns (indent, module_name, symbols, end_idx) or None.
    Handles both single-line and multi-line parenthesized imports.
    """
    line = lines[start_idx]

    # Check for: from metainformant.<mod> import (
    m = _MULTILINE_START_RE.match(line)
    if m:
        indent = m.group(1)
        module_name = m.group(2)
        symbols = []
        idx = start_idx + 1
        while idx < len(lines):
            end_m = _MULTILINE_END_RE.match(lines[idx])
            if end_m:
                return indent, module_name, symbols, idx
            # Parse symbol names from the line (strip comments, commas, whitespace)
            sym_line = lines[idx].split("#")[0].strip().rstrip(",")
            if sym_line:
                # Handle "name as alias" and plain "name"
                for part in sym_line.split(","):
                    part = part.strip()
                    if part:
                        symbols.append(part)
            idx += 1
        return None  # Unterminated

    # Check for single-line: from metainformant.<mod> import a, b, c
    single_re = re.match(
        r"^(\s*)from\s+metainformant\.(\w+)\s+import\s+(.+)$", line
    )
    if single_re:
        indent = single_re.group(1)
        module_name = single_re.group(2)
        names_str = single_re.group(3).split("#")[0].strip()

        # Skip if it looks like a deeper import (from metainformant.mod.sub import ...)
        # We already checked it's from metainformant.<word> only (no dots in \w+)
        # But also handle parenthesized single-line: from metainformant.mod import (a, b)
        if names_str.startswith("(") and names_str.endswith(")"):
            names_str = names_str[1:-1]

        symbols = []
        for part in names_str.split(","):
            part = part.strip()
            if part:
                symbols.append(part)

        return indent, module_name, symbols, start_idx

    return None


def rewrite_test_imports(*, apply: bool = False) -> int:
    """Rewrite shortcut imports in test files and cross-module source files.

    Returns the number of files that would be (or were) changed.
    """
    changed = 0

    # Collect all .py files to scan
    files_to_scan: list[Path] = []
    files_to_scan.extend(sorted(TESTS.glob("*.py")))
    # Also scan source files for cross-module shortcut imports
    for py_file in sorted(SRC.rglob("*.py")):
        if "__pycache__" in str(py_file):
            continue
        if py_file.name == "__init__.py":
            continue
        files_to_scan.append(py_file)

    for filepath in files_to_scan:
        try:
            content = filepath.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue

        lines = content.split("\n")
        new_lines: list[str] = []
        idx = 0
        file_changed = False
        # Track which module this file belongs to (for skipping intra-module)
        file_module = None
        try:
            rel_to_src = filepath.relative_to(SRC)
            file_module = rel_to_src.parts[0] if rel_to_src.parts else None
        except ValueError:
            pass  # Not in src, probably a test file

        while idx < len(lines):
            line = lines[idx]

            # Try to parse as a shortcut import from metainformant.<module>
            parsed = parse_import_block(lines, idx)

            if parsed is None:
                new_lines.append(line)
                idx += 1
                continue

            indent, module_name, symbols, end_idx = parsed

            # Skip if this module has no canonical map
            if module_name not in CANONICAL_MAP:
                for i in range(idx, end_idx + 1):
                    new_lines.append(lines[i])
                idx = end_idx + 1
                continue

            # Skip intra-module imports in source files
            if file_module and file_module == module_name:
                for i in range(idx, end_idx + 1):
                    new_lines.append(lines[i])
                idx = end_idx + 1
                continue

            symbol_map = CANONICAL_MAP[module_name]

            # Check if ANY symbol in this import needs rewriting
            needs_rewrite = False
            for sym in symbols:
                # Handle "name as alias"
                name = sym.split(" as ")[0].strip() if " as " in sym else sym.strip()
                if name in symbol_map:
                    needs_rewrite = True
                    break

            if not needs_rewrite:
                for i in range(idx, end_idx + 1):
                    new_lines.append(lines[i])
                idx = end_idx + 1
                continue

            # Group symbols by canonical path
            grouped: dict[str, list[str]] = defaultdict(list)
            unchanged: list[str] = []

            for sym in symbols:
                sym = sym.strip()
                if not sym:
                    continue

                # Parse "name as alias"
                if " as " in sym:
                    name, alias = sym.split(" as ", 1)
                    name = name.strip()
                    alias = alias.strip()
                    sym_display = f"{name} as {alias}"
                else:
                    name = sym
                    sym_display = sym

                if name in symbol_map:
                    canonical = symbol_map[name]
                    grouped[canonical].append(sym_display)
                else:
                    unchanged.append(sym_display)

            # Generate new import lines
            new_import_lines = []

            # Keep unchanged symbols as original import
            if unchanged:
                if len(unchanged) <= 3:
                    names_str = ", ".join(unchanged)
                    new_import_lines.append(
                        f"{indent}from metainformant.{module_name} import {names_str}"
                    )
                else:
                    new_import_lines.append(
                        f"{indent}from metainformant.{module_name} import ("
                    )
                    for u in unchanged:
                        new_import_lines.append(f"{indent}    {u},")
                    new_import_lines.append(f"{indent})")

            # Add canonical imports, grouped by path
            for canonical_path in sorted(grouped.keys()):
                syms = sorted(grouped[canonical_path])
                if len(syms) <= 3:
                    names_str = ", ".join(syms)
                    new_import_lines.append(
                        f"{indent}from {canonical_path} import {names_str}"
                    )
                else:
                    new_import_lines.append(
                        f"{indent}from {canonical_path} import ("
                    )
                    for s in syms:
                        new_import_lines.append(f"{indent}    {s},")
                    new_import_lines.append(f"{indent})")

            new_lines.extend(new_import_lines)
            file_changed = True
            idx = end_idx + 1

        if file_changed:
            changed += 1
            rel = filepath.relative_to(ROOT)

            if apply:
                filepath.write_text("\n".join(new_lines), encoding="utf-8")
                print(f"  [WRITE] {rel}")
            else:
                print(f"  [WOULD WRITE] {rel}")

    return changed


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite __init__.py and test imports to canonical paths"
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually write changes (default is dry-run)",
    )
    parser.add_argument(
        "--init-only",
        action="store_true",
        help="Only rewrite __init__.py files",
    )
    parser.add_argument(
        "--tests-only",
        action="store_true",
        help="Only rewrite test/source imports",
    )
    args = parser.parse_args()

    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"=== Import Rewriter ({mode}) ===\n")

    total_changed = 0

    if not args.tests_only:
        print("--- Part 1: Rewriting __init__.py files ---")
        n = rewrite_init_files(apply=args.apply)
        total_changed += n
        print(f"  -> {n} __init__.py files {'changed' if args.apply else 'would change'}\n")

    if not args.init_only:
        print("--- Part 2: Rewriting test/source imports ---")
        n = rewrite_test_imports(apply=args.apply)
        total_changed += n
        print(f"  -> {n} files {'changed' if args.apply else 'would change'}\n")

    print(f"=== Total: {total_changed} files {'changed' if args.apply else 'would change'} ===")

    if not args.apply and total_changed > 0:
        print("\nRun with --apply to write changes.")


if __name__ == "__main__":
    main()
