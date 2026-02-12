#!/usr/bin/env python3
"""Rewrite test imports from module-level shortcuts to canonical submodule paths.

This script reads each test file, finds imports from metainformant.<module> that use
shortcut imports (which no longer exist after __init__.py cleanup), and rewrites them
to use canonical submodule paths.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Complete mapping: (module, symbol) -> canonical_module_path
IMPORT_MAP: dict[str, dict[str, str]] = {
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
        "config": "metainformant.core.utils.config",
        "cpu_count": "metainformant.core.execution.parallel",
        "create_sample_config": "metainformant.core.execution.workflow",
        "data": "metainformant.core.io.cache",
        "discover_configs": "metainformant.core.execution.discovery",
        "discover_functions": "metainformant.core.execution.discovery",
        "discovery": "metainformant.core.execution.discovery",
        "download": "metainformant.core.io.download",
        "download_and_process_data": "metainformant.core.execution.workflow",
        "download_file": "metainformant.core.io.io",
        "dump_json": "metainformant.core.io.io",
        "errors": "metainformant.core.execution.workflow",
        "load_json": "metainformant.core.io.io",
        "load_yaml": "metainformant.core.io.io",
        "logging": "metainformant.core.utils.logging",
        "parallel": "metainformant.core.execution.parallel",
        "parallel_batch": "metainformant.core.execution.parallel",
        "paths": "metainformant.core.io.paths",
        "rate_limiter": "metainformant.core.utils.timing",
        "run_config_based_workflow": "metainformant.core.execution.workflow",
        "safe_write_bytes": "metainformant.core.io.atomic",
        "safe_write_text": "metainformant.core.io.atomic",
        "symbols": "metainformant.core.utils.symbols",
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
        "validation": "metainformant.core.data.validation",
        "verify_checksum": "metainformant.core.io.checksums",
        "verify_checksum_file": "metainformant.core.io.checksums",
        "write_checksum_file": "metainformant.core.io.checksums",
        "disk": "metainformant.core.io.disk",
        "hash": "metainformant.core.utils.hash",
        "optional_deps": "metainformant.core.utils.optional_deps",
        "progress": "metainformant.core.utils.progress",
        "workflow": "metainformant.core.execution.workflow",
        "io_module": "metainformant.core.io.io",
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
        # Submodule imports (these are packages, import as submodules)
        "alignment": "metainformant.dna.alignment",
        "msa": "metainformant.dna.alignment.msa",
        "distances": "metainformant.dna.alignment.distances",
        "transcription": "metainformant.dna.expression.transcription",
        "translation": "metainformant.dna.expression.translation",
        "gene_prediction": "metainformant.dna.annotation.gene_prediction",
        "functional": "metainformant.dna.annotation.functional",
        "ncbi": "metainformant.dna.external.ncbi",
        "entrez": "metainformant.dna.external.entrez",
        "genomes": "metainformant.dna.external.genomes",
        "fastq": "metainformant.dna.io.fastq",
        "composition": "metainformant.dna.sequence.composition",
        "motifs": "metainformant.dna.sequence.motifs",
        "restriction": "metainformant.dna.sequence.restriction",
        "calling": "metainformant.dna.variation.calling",
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
    "life_events": {
        "Event": "metainformant.life_events.core.events",
        "EventSequence": "metainformant.life_events.core.events",
        "EventDatabase": "metainformant.life_events.core.events",
        "EventSequencePredictor": "metainformant.life_events.models.models",
        "EnsemblePredictor": "metainformant.life_events.models.models",
        "GRUSequenceModel": "metainformant.life_events.models.models",
        "LSTMSequenceModel": "metainformant.life_events.models.models",
        "MultiTaskPredictor": "metainformant.life_events.models.models",
        "SurvivalPredictor": "metainformant.life_events.models.models",
        "attention_weights": "metainformant.life_events.models.models",
        "LifeEventsWorkflowConfig": "metainformant.life_events.core.config",
        "load_life_events_config": "metainformant.life_events.core.config",
        "analyze_life_course": "metainformant.life_events.workflow.workflow",
        "compare_populations": "metainformant.life_events.workflow.workflow",
        "intervention_analysis": "metainformant.life_events.workflow.workflow",
        "temporal_patterns": "metainformant.life_events.analysis.interpretability",
        "event_importance": "metainformant.life_events.analysis",
        "feature_attribution": "metainformant.life_events.analysis",
        "learn_event_embeddings": "metainformant.life_events.models.embeddings",
        "biological_embedding": "metainformant.life_events.models.embeddings",
        "domain_specific_embeddings": "metainformant.life_events.models.embeddings",
        "add_temporal_noise": "metainformant.life_events.core.utils",
        "convert_sequences_to_tokens": "metainformant.life_events.core.utils",
        "generate_cohort_sequences": "metainformant.life_events.core.utils",
        "generate_event_chain": "metainformant.life_events.core.utils",
        "generate_synthetic_life_events": "metainformant.life_events.core.utils",
        "get_event_statistics": "metainformant.life_events.core.utils",
        "load_sequences_from_json": "metainformant.life_events.core.utils",
        "sequence_embeddings": "metainformant.life_events.core.utils",
        "validate_sequence": "metainformant.life_events.core.utils",
        "plot_event_timeline": "metainformant.life_events.visualization.visualization",
        "plot_attention_heatmap": "metainformant.life_events.visualization.visualization",
        "plot_domain_distribution": "metainformant.life_events.visualization.visualization",
        "plot_domain_timeline": "metainformant.life_events.visualization.visualization",
        "plot_embedding_clusters": "metainformant.life_events.visualization.visualization",
        "plot_event_cooccurrence": "metainformant.life_events.visualization.visualization",
        "plot_event_embeddings": "metainformant.life_events.visualization.visualization",
        "plot_event_frequency_heatmap": "metainformant.life_events.visualization.visualization",
        "plot_intervention_effects": "metainformant.life_events.visualization.visualization",
        "plot_outcome_distribution": "metainformant.life_events.visualization.visualization",
        "plot_population_comparison": "metainformant.life_events.visualization.visualization",
        "plot_prediction_accuracy": "metainformant.life_events.visualization.visualization",
        "plot_prediction_importance": "metainformant.life_events.visualization.visualization",
        "plot_sequence_length_distribution": "metainformant.life_events.visualization.visualization",
        "plot_sequence_similarity": "metainformant.life_events.visualization.visualization",
        "plot_temporal_density": "metainformant.life_events.visualization.visualization",
        "plot_temporal_patterns": "metainformant.life_events.visualization.visualization",
        "plot_transition_network": "metainformant.life_events.visualization.visualization",
        "kaplan_meier_estimator": "metainformant.life_events.survival.time_to_event",
        "cox_ph_model": "metainformant.life_events.survival.time_to_event",
        "competing_risks": "metainformant.life_events.survival.time_to_event",
        "recurrent_events": "metainformant.life_events.survival.time_to_event",
        "time_varying_covariates": "metainformant.life_events.survival.time_to_event",
    },
    "math": {
        "correlation": "metainformant.math.core.utilities",
        "correlation_coefficient": "metainformant.math.core.utilities",
        "covariance": "metainformant.math.core.utilities",
        "fisher_exact_test": "metainformant.math.core.utilities",
        "jensen_shannon_divergence": "metainformant.math.core.utilities",
        "linear_regression": "metainformant.math.core.utilities",
        "r_squared": "metainformant.math.core.utilities",
        "shannon_entropy": "metainformant.math.core.utilities",
        "abc_rejection": "metainformant.math.bayesian.inference",
        "compute_bayes_factor": "metainformant.math.bayesian.inference",
        "compute_dic": "metainformant.math.bayesian.inference",
        "compute_waic": "metainformant.math.bayesian.inference",
        "conjugate_beta_binomial": "metainformant.math.bayesian.inference",
        "conjugate_normal": "metainformant.math.bayesian.inference",
        "metropolis_hastings": "metainformant.math.bayesian.inference",
        "ddm_analytic_accuracy": "metainformant.math.decision_theory.ddm",
        "ddm_mean_decision_time": "metainformant.math.decision_theory.ddm",
        "basic_reproduction_number": "metainformant.math.epidemiology.models",
        "effective_reproduction_number": "metainformant.math.epidemiology.models",
        "herd_immunity_threshold": "metainformant.math.epidemiology.models",
        "seir_step": "metainformant.math.epidemiology.models",
        "sir_step": "metainformant.math.epidemiology.models",
        "sis_step": "metainformant.math.epidemiology.models",
        "logistic_map": "metainformant.math.evolutionary_dynamics.core",
        "lotka_volterra_step": "metainformant.math.evolutionary_dynamics.core",
        "replicator_derivative": "metainformant.math.evolutionary_dynamics.core",
        "replicator_step": "metainformant.math.evolutionary_dynamics.core",
        "expected_pairwise_diversity": "metainformant.math.population_genetics.coalescent",
        "expected_segregating_sites": "metainformant.math.population_genetics.coalescent",
        "expected_time_to_mrca": "metainformant.math.population_genetics.coalescent",
        "expected_total_branch_length": "metainformant.math.population_genetics.coalescent",
        "tajima_constants": "metainformant.math.population_genetics.coalescent",
        "tajimas_D": "metainformant.math.population_genetics.coalescent",
        "watterson_theta": "metainformant.math.population_genetics.coalescent",
        "hardy_weinberg_genotype_freqs": "metainformant.math.population_genetics.core",
        "heterozygosity_decay": "metainformant.math.population_genetics.core",
        "inbreeding_coefficient": "metainformant.math.population_genetics.core",
        "mutation_selection_balance_dominant": "metainformant.math.population_genetics.core",
        "mutation_selection_balance_recessive": "metainformant.math.population_genetics.core",
        "island_model_update": "metainformant.math.decision_theory.ddm",
        "effective_size_sex_ratio": "metainformant.math.population_genetics.effective_size",
        "harmonic_mean_effective_size": "metainformant.math.population_genetics.effective_size",
        "fst_from_allele_freqs": "metainformant.math.population_genetics.fst",
        "fst_from_heterozygosity": "metainformant.math.population_genetics.fst",
        "variance": "metainformant.math.population_genetics.fst",
        "haldane_c_to_d": "metainformant.math.population_genetics.ld",
        "haldane_d_to_c": "metainformant.math.population_genetics.ld",
        "kosambi_c_to_d": "metainformant.math.population_genetics.ld",
        "kosambi_d_to_c": "metainformant.math.population_genetics.ld",
        "ld_coefficients": "metainformant.math.population_genetics.ld",
        "ld_decay_r2": "metainformant.math.population_genetics.ld",
        "breeders_equation_response": "metainformant.math.population_genetics.selection",
        "kin_selection_response": "metainformant.math.population_genetics.selection",
        "mutation_update": "metainformant.math.population_genetics.selection",
        "relative_fitness": "metainformant.math.population_genetics.selection",
        "selection_differential": "metainformant.math.population_genetics.selection",
        "selection_gradient": "metainformant.math.population_genetics.selection",
        "selection_intensity": "metainformant.math.population_genetics.selection",
        "selection_update": "metainformant.math.population_genetics.selection",
        "effective_size_from_family_size_variance": "metainformant.math.population_genetics.statistics",
        "equilibrium_heterozygosity_infinite_alleles": "metainformant.math.population_genetics.statistics",
        "expected_r2_from_Ne_c": "metainformant.math.population_genetics.statistics",
        "fixation_probability": "metainformant.math.population_genetics.statistics",
        "kurtosis": "metainformant.math.population_genetics.statistics",
        "skewness": "metainformant.math.population_genetics.statistics",
        "standard_deviation": "metainformant.math.population_genetics.statistics",
        "lande_equation_response": "metainformant.math.quantitative_genetics.core",
        "narrow_sense_heritability": "metainformant.math.quantitative_genetics.core",
        "realized_heritability": "metainformant.math.quantitative_genetics.core",
        "delta_mean_trait": "metainformant.math.quantitative_genetics.price",
        "expectation": "metainformant.math.quantitative_genetics.price",
        "price_equation": "metainformant.math.quantitative_genetics.price",
        "weighted_correlation": "metainformant.math.quantitative_genetics.price",
        "weighted_covariance": "metainformant.math.quantitative_genetics.price",
        "weighted_variance": "metainformant.math.quantitative_genetics.price",
    },
    "information": {
        "shannon_entropy": "metainformant.information.metrics.core.syntactic",
        "shannon_entropy_from_counts": "metainformant.information.metrics.core.syntactic",
        "joint_entropy": "metainformant.information.metrics.core.syntactic",
        "conditional_entropy": "metainformant.information.metrics.core.syntactic",
        "mutual_information": "metainformant.information.metrics.core.syntactic",
        "conditional_mutual_information": "metainformant.information.metrics.core.syntactic",
        "kl_divergence": "metainformant.information.metrics.core.syntactic",
        "cross_entropy": "metainformant.information.metrics.core.syntactic",
        "jensen_shannon_divergence": "metainformant.information.metrics.core.syntactic",
        "total_correlation": "metainformant.information.metrics.core.syntactic",
        "transfer_entropy": "metainformant.information.metrics.core.syntactic",
        "renyi_entropy": "metainformant.information.metrics.core.syntactic",
        "tsallis_entropy": "metainformant.information.metrics.core.syntactic",
        "normalized_mutual_information": "metainformant.information.metrics.core.syntactic",
        "information_coefficient": "metainformant.information.metrics.core.syntactic",
        "information_content": "metainformant.information.metrics.advanced.semantic",
        "information_content_from_annotations": "metainformant.information.metrics.advanced.semantic",
        "semantic_entropy": "metainformant.information.metrics.advanced.semantic",
        "semantic_similarity": "metainformant.information.metrics.advanced.semantic",
        "semantic_similarity_matrix": "metainformant.information.metrics.advanced.semantic",
        "semantic_distance": "metainformant.information.metrics.advanced.semantic",
        "term_specificity": "metainformant.information.metrics.advanced.semantic",
        "ontology_complexity": "metainformant.information.metrics.advanced.semantic",
        "term_redundancy": "metainformant.information.metrics.advanced.semantic",
        "annotation_specificity": "metainformant.information.metrics.advanced.semantic",
        "differential_entropy": "metainformant.information.metrics.core.continuous",
        "mutual_information_continuous": "metainformant.information.metrics.core.continuous",
        "kl_divergence_continuous": "metainformant.information.metrics.core.continuous",
        "entropy_estimation": "metainformant.information.metrics.core.continuous",
        "copula_entropy": "metainformant.information.metrics.core.continuous",
        "transfer_entropy_continuous": "metainformant.information.metrics.core.continuous",
        "conditional_entropy_continuous": "metainformant.information.metrics.core.continuous",
        "information_flow_network": "metainformant.information.metrics.core.continuous",
        "entropy_estimator": "metainformant.information.metrics.core.estimation",
        "mutual_information_estimator": "metainformant.information.metrics.core.estimation",
        "kl_divergence_estimator": "metainformant.information.metrics.core.estimation",
        "bias_correction": "metainformant.information.metrics.core.estimation",
        "entropy_bootstrap_confidence": "metainformant.information.metrics.core.estimation",
        "panzeri_treves_bias_correction": "metainformant.information.metrics.core.estimation",
        "entropy_rate_estimator": "metainformant.information.metrics.core.estimation",
        "channel_capacity": "metainformant.information.metrics.advanced.channel",
        "channel_mutual_information": "metainformant.information.metrics.advanced.channel",
        "information_bottleneck": "metainformant.information.metrics.advanced.channel",
        "noisy_channel_capacity": "metainformant.information.metrics.advanced.channel",
        "rate_distortion": "metainformant.information.metrics.advanced.channel",
        "co_information": "metainformant.information.metrics.advanced.decomposition",
        "dual_total_correlation": "metainformant.information.metrics.advanced.decomposition",
        "o_information": "metainformant.information.metrics.advanced.decomposition",
        "partial_information_decomposition": "metainformant.information.metrics.advanced.decomposition",
        "redundant_information": "metainformant.information.metrics.advanced.decomposition",
        "synergistic_information": "metainformant.information.metrics.advanced.decomposition",
        "unique_information": "metainformant.information.metrics.advanced.decomposition",
        "entropy_power_inequality": "metainformant.information.metrics.advanced.geometry",
        "exponential_family_entropy": "metainformant.information.metrics.advanced.geometry",
        "fisher_rao_distance": "metainformant.information.metrics.advanced.geometry",
        "hellinger_distance": "metainformant.information.metrics.advanced.geometry",
        "information_dimension": "metainformant.information.metrics.advanced.geometry",
        "information_projection": "metainformant.information.metrics.advanced.geometry",
        "natural_gradient": "metainformant.information.metrics.advanced.geometry",
        "statistical_divergence": "metainformant.information.metrics.advanced.geometry",
        "entropy_confidence_interval": "metainformant.information.metrics.advanced.hypothesis",
        "entropy_rate_test": "metainformant.information.metrics.advanced.hypothesis",
        "independence_test": "metainformant.information.metrics.advanced.hypothesis",
        "information_significance_filter": "metainformant.information.metrics.advanced.hypothesis",
        "mi_permutation_test": "metainformant.information.metrics.advanced.hypothesis",
        "analyze_sequence_information": "metainformant.information.metrics.analysis",
        "binding_information": "metainformant.information.metrics.analysis",
        "compare_sequences_information": "metainformant.information.metrics.analysis",
        "fisher_information": "metainformant.information.metrics.analysis",
        "fisher_information_matrix": "metainformant.information.metrics.analysis",
        "information_profile": "metainformant.information.metrics.analysis",
        "information_signature": "metainformant.information.metrics.analysis",
        "interaction_information": "metainformant.information.metrics.analysis",
        "lautum_information": "metainformant.information.metrics.analysis",
        "relative_information_gain": "metainformant.information.metrics.analysis",
        "variation_of_information": "metainformant.information.metrics.analysis",
        "network_entropy": "metainformant.information.integration.networks",
        "information_flow": "metainformant.information.integration.networks",
        "information_community_detection": "metainformant.information.integration.networks",
        "network_information_centrality": "metainformant.information.integration.networks",
        "network_motif_information": "metainformant.information.integration.networks",
        "information_graph_distance": "metainformant.information.integration.networks",
        "batch_entropy_analysis": "metainformant.information.workflow.workflows",
        "information_workflow": "metainformant.information.workflow.workflows",
        "compare_datasets": "metainformant.information.workflow.workflows",
        "information_report": "metainformant.information.workflow.workflows",
        "granger_causality": "metainformant.information.network_info.information_flow",
        "mutual_information_network": "metainformant.information.network_info.information_flow",
    },
    "visualization": {
        "apply_theme": "metainformant.visualization.config.themes",
        "list_themes": "metainformant.visualization.config.themes",
        "theme": "metainformant.visualization.config.themes",
        "get_theme": "metainformant.visualization.config.themes",
        "reset_theme": "metainformant.visualization.config.themes",
        "register_theme": "metainformant.visualization.config.themes",
        "categorical": "metainformant.visualization.config.palettes",
        "chromosome_palette": "metainformant.visualization.config.palettes",
        "config": "metainformant.visualization.config",
        "genomic_overview": "metainformant.visualization.dashboards.composite",
        "multi_panel": "metainformant.visualization.dashboards.composite",
        "qc_summary": "metainformant.visualization.dashboards.composite",
        "interactive_scatter": "metainformant.visualization.dashboards.interactive",
        "interactive_volcano": "metainformant.visualization.dashboards.interactive",
        "interactive_heatmap": "metainformant.visualization.dashboards.interactive",
        "interactive_manhattan": "metainformant.visualization.dashboards.interactive",
        "plot_phylo_tree": "metainformant.visualization.genomics.trees",
        "circular_tree_plot": "metainformant.visualization.genomics.trees",
        "unrooted_tree_plot": "metainformant.visualization.genomics.trees",
    },
    "gwas": {
        "GWASWorkflowConfig": "metainformant.gwas.workflow.workflow_config",
        "PlotStyle": "metainformant.gwas.visualization.config",
        "admixture_plot": "metainformant.gwas.visualization.population.population_admixture",
    },
    "rna": {
        "RNAPipelineConfig": "metainformant.rna.core.configs",
        "AmalgkitRunLayout": "metainformant.rna.core.configs",
        "ProgressTracker": "metainformant.rna.engine.progress_tracker",
        "AmalgkitWorkflowConfig": "metainformant.rna.engine.workflow",
    },
}


def parse_import_line(line: str) -> tuple[str, list[str]] | None:
    """Parse 'from metainformant.X import a, b, c' into (module, [symbols])."""
    m = re.match(
        r"^from\s+metainformant\.(\w+)\s+import\s+(.+)$",
        line.strip().rstrip("\\").strip(),
    )
    if not m:
        return None
    module = m.group(1)
    symbols_str = m.group(2)
    # Handle parenthesized imports
    symbols_str = symbols_str.strip("()")
    symbols = [s.strip().rstrip(",") for s in symbols_str.split(",") if s.strip()]
    return module, symbols


def rewrite_file(filepath: Path, dry_run: bool = True) -> bool:
    """Rewrite imports in a single file. Returns True if changes were made."""
    content = filepath.read_text()
    lines = content.split("\n")
    new_lines: list[str] = []
    changed = False
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Handle multi-line imports (with parentheses or backslash continuation)
        if stripped.startswith("from metainformant.") and "import" in stripped:
            # Collect full import statement
            full_import = stripped
            if "(" in stripped and ")" not in stripped:
                # Multi-line parenthesized import
                i += 1
                while i < len(lines) and ")" not in lines[i]:
                    full_import += " " + lines[i].strip()
                    i += 1
                if i < len(lines):
                    full_import += " " + lines[i].strip()
            elif stripped.endswith("\\"):
                while i < len(lines) and lines[i].strip().endswith("\\"):
                    i += 1
                    if i < len(lines):
                        full_import += " " + lines[i].strip()

            # Parse the import
            # Clean up parentheses
            clean = full_import.replace("(", "").replace(")", "").replace("\\", "")
            clean = re.sub(r"\s+", " ", clean).strip()
            parsed = parse_import_line(clean)

            if parsed:
                module, symbols = parsed
                if module in IMPORT_MAP:
                    module_map = IMPORT_MAP[module]
                    # Group symbols by their canonical module
                    groups: dict[str, list[str]] = {}
                    unmatched: list[str] = []
                    for sym in symbols:
                        sym = sym.strip()
                        if not sym or sym.startswith("#"):
                            continue
                        if sym in module_map:
                            target = module_map[sym]
                            groups.setdefault(target, []).append(sym)
                        else:
                            unmatched.append(sym)

                    if groups:
                        changed = True
                        # Generate new import lines
                        for target_mod, syms in sorted(groups.items()):
                            syms_str = ", ".join(sorted(syms))
                            new_lines.append(f"from {target_mod} import {syms_str}")
                        # Keep any unmatched symbols as original import
                        if unmatched:
                            syms_str = ", ".join(unmatched)
                            new_lines.append(f"from metainformant.{module} import {syms_str}")
                        i += 1
                        continue

            # No rewrite needed
            new_lines.append(lines[i] if i < len(lines) else "")
            i += 1
            continue

        new_lines.append(line)
        i += 1

    if changed:
        new_content = "\n".join(new_lines)
        if dry_run:
            print(f"  [DRY RUN] Would rewrite {filepath}")
        else:
            filepath.write_text(new_content)
            print(f"  [APPLIED] Rewrote {filepath}")
    return changed


def main() -> None:
    dry_run = "--apply" not in sys.argv
    if dry_run:
        print("DRY RUN MODE - use --apply to write changes\n")

    tests_dir = Path("tests")
    if not tests_dir.exists():
        print("ERROR: Run from project root")
        sys.exit(1)

    # Also check source files for cross-module imports
    src_dir = Path("src/metainformant")

    test_files = sorted(tests_dir.glob("test_*.py"))
    src_files = sorted(src_dir.rglob("*.py"))

    total_changed = 0

    print(f"Processing {len(test_files)} test files...")
    for tf in test_files:
        if rewrite_file(tf, dry_run):
            total_changed += 1

    print(f"\nProcessing {len(src_files)} source files...")
    for sf in src_files:
        # Skip __init__.py files (already cleaned)
        if sf.name == "__init__.py":
            continue
        if rewrite_file(sf, dry_run):
            total_changed += 1

    print(f"\nTotal files {'that would be' if dry_run else ''} changed: {total_changed}")


if __name__ == "__main__":
    main()
