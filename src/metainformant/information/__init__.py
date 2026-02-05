"""Information theory analysis module for METAINFORMANT.

This module provides comprehensive information-theoretic methods for biological data,
including syntactic information (Shannon entropy, mutual information),
semantic information (information content, semantic similarity),
continuous information theory, estimation methods, and analysis workflows.
"""

from __future__ import annotations

# Import subpackages
from . import integration, metrics, network_info, workflow
from .integration import integration as integration_module
from .integration import networks

# Network functions
from .integration.networks import (
    information_community_detection,
    information_flow,
    information_graph_distance,
    network_entropy,
    network_information_centrality,
    network_motif_information,
)

# Import modules from subpackages
from .metrics import (
    advanced,
    analysis,
    channel,
    continuous,
    decomposition,
    estimation,
    geometry,
    hypothesis,
    semantic,
    syntactic,
)

# Advanced measures
from .metrics.advanced import (
    binding_information,
    fisher_information,
    fisher_information_matrix,
    interaction_information,
    lautum_information,
    relative_information_gain,
    variation_of_information,
)

# Analysis functions
from .metrics.analysis import (
    analyze_sequence_information,
    compare_sequences_information,
    information_profile,
    information_signature,
)

# Channel capacity functions
from .metrics.channel import (
    channel_capacity,
    channel_mutual_information,
    information_bottleneck,
    noisy_channel_capacity,
    rate_distortion,
)

# Continuous information functions
from .metrics.continuous import (
    conditional_entropy_continuous,
    copula_entropy,
    differential_entropy,
    entropy_estimation,
    information_flow_network,
    kl_divergence_continuous,
    mutual_information_continuous,
    transfer_entropy_continuous,
)

# Partial Information Decomposition functions
from .metrics.decomposition import (
    co_information,
    dual_total_correlation,
    o_information,
    partial_information_decomposition,
    redundant_information,
    synergistic_information,
    unique_information,
)

# Estimation functions
from .metrics.estimation import (
    bias_correction,
    entropy_bootstrap_confidence,
    entropy_estimator,
    entropy_rate_estimator,
    kl_divergence_estimator,
    mutual_information_estimator,
    panzeri_treves_bias_correction,
)

# Information geometry functions
from .metrics.geometry import (
    entropy_power_inequality,
    exponential_family_entropy,
    fisher_rao_distance,
    hellinger_distance,
    information_dimension,
    information_projection,
    natural_gradient,
    statistical_divergence,
)

# Hypothesis testing functions
from .metrics.hypothesis import (
    entropy_confidence_interval,
    entropy_rate_test,
    independence_test,
    information_significance_filter,
    mi_permutation_test,
)

# Semantic information functions
from .metrics.semantic import (
    annotation_specificity,
    information_content,
    information_content_from_annotations,
    ontology_complexity,
    semantic_distance,
    semantic_entropy,
    semantic_similarity,
    semantic_similarity_matrix,
    term_redundancy,
    term_specificity,
)

# Syntactic information functions
from .metrics.syntactic import (
    conditional_entropy,
    conditional_mutual_information,
    cross_entropy,
    information_coefficient,
    jensen_shannon_divergence,
    joint_entropy,
    kl_divergence,
    mutual_information,
    normalized_mutual_information,
    renyi_entropy,
    shannon_entropy,
    shannon_entropy_from_counts,
    total_correlation,
    transfer_entropy,
    tsallis_entropy,
)

# Network information flow
from .network_info.information_flow import (
    granger_causality,
)
from .network_info.information_flow import information_flow_network as build_information_flow_network
from .network_info.information_flow import (
    mutual_information_network,
)
from .network_info.information_flow import network_entropy as network_von_neumann_entropy
from .network_info.information_flow import transfer_entropy as transfer_entropy_network
from .workflow import workflows as workflow_module

# Workflow functions
from .workflow.workflows import (
    batch_entropy_analysis,
    compare_datasets,
    information_report,
    information_workflow,
)

__all__ = [
    # Subpackages
    "integration",
    "metrics",
    "workflow",
    # Modules
    "analysis",
    "continuous",
    "estimation",
    "semantic",
    "syntactic",
    "integration_module",
    "networks",
    "workflow_module",
    # Syntactic
    "shannon_entropy",
    "shannon_entropy_from_counts",
    "joint_entropy",
    "conditional_entropy",
    "mutual_information",
    "conditional_mutual_information",
    "kl_divergence",
    "cross_entropy",
    "jensen_shannon_divergence",
    "total_correlation",
    "transfer_entropy",
    "renyi_entropy",
    "tsallis_entropy",
    "normalized_mutual_information",
    "information_coefficient",
    # Semantic
    "information_content",
    "information_content_from_annotations",
    "semantic_entropy",
    "semantic_similarity",
    "semantic_similarity_matrix",
    "semantic_distance",
    "term_specificity",
    "ontology_complexity",
    "term_redundancy",
    "annotation_specificity",
    # Analysis
    "information_profile",
    "information_signature",
    "analyze_sequence_information",
    "compare_sequences_information",
    # Continuous
    "differential_entropy",
    "mutual_information_continuous",
    "kl_divergence_continuous",
    "entropy_estimation",
    "copula_entropy",
    "transfer_entropy_continuous",
    "conditional_entropy_continuous",
    "information_flow_network",
    # Estimation
    "entropy_estimator",
    "mutual_information_estimator",
    "kl_divergence_estimator",
    "bias_correction",
    "entropy_bootstrap_confidence",
    "panzeri_treves_bias_correction",
    "entropy_rate_estimator",
    # Advanced
    "fisher_information",
    "fisher_information_matrix",
    "relative_information_gain",
    "variation_of_information",
    "interaction_information",
    "binding_information",
    "lautum_information",
    # Networks
    "network_entropy",
    "information_flow",
    "information_community_detection",
    "network_information_centrality",
    "network_motif_information",
    "information_graph_distance",
    # Hypothesis testing
    "mi_permutation_test",
    "independence_test",
    "entropy_confidence_interval",
    "information_significance_filter",
    "entropy_rate_test",
    # Channel capacity
    "channel_capacity",
    "rate_distortion",
    "information_bottleneck",
    "channel_mutual_information",
    "noisy_channel_capacity",
    # Information geometry
    "fisher_rao_distance",
    "natural_gradient",
    "information_projection",
    "statistical_divergence",
    "exponential_family_entropy",
    "hellinger_distance",
    "entropy_power_inequality",
    "information_dimension",
    # Partial Information Decomposition
    "co_information",
    "partial_information_decomposition",
    "unique_information",
    "redundant_information",
    "synergistic_information",
    "dual_total_correlation",
    "o_information",
    # Workflows
    "batch_entropy_analysis",
    "information_workflow",
    "compare_datasets",
    "information_report",
    # Network information flow
    "network_info",
    "transfer_entropy_network",
    "granger_causality",
    "network_von_neumann_entropy",
    "build_information_flow_network",
    "mutual_information_network",
]
