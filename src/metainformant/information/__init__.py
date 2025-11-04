"""Information theory methods for biological data analysis.

This module provides comprehensive information theory tools for analyzing
biological data, including both syntactic (Shannon entropy, mutual information)
and semantic (information content, semantic similarity) information measures.

Modules:
- syntactic: Shannon entropy, mutual information, KL divergence, etc.
- semantic: Information content, semantic similarity, semantic entropy
- analysis: High-level analysis functions and workflows
- continuous: Continuous information theory methods
- estimation: Entropy and MI estimation with bias correction
- workflows: Batch processing and workflow functions
- integration: Integration with other modules (DNA, RNA, single-cell, etc.)
- visualization: Integration with visualization module
- networks: Integration with networks module for information-theoretic analysis
"""

from __future__ import annotations

from .analysis import (
    analyze_sequence_information,
    compare_sequences_information,
    information_profile,
    information_signature,
)
from .continuous import (
    differential_entropy,
    entropy_estimation,
    kl_divergence_continuous,
    mutual_information_continuous,
)
from .estimation import (
    bias_correction,
    entropy_estimator,
    kl_divergence_estimator,
    mutual_information_estimator,
)
from .integration import (
    dna_integration,
    ml_integration,
    multiomics_integration,
    rna_integration,
    singlecell_integration,
)
from .semantic import (
    information_content,
    semantic_entropy,
    semantic_similarity,
    semantic_similarity_matrix,
)
from .syntactic import (
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
    total_correlation,
    transfer_entropy,
    tsallis_entropy,
)
from .workflows import (
    batch_entropy_analysis,
    compare_datasets,
    information_report,
    information_workflow,
)

__all__ = [
    # Syntactic information theory
    "shannon_entropy",
    "joint_entropy",
    "conditional_entropy",
    "mutual_information",
    "conditional_mutual_information",
    "kl_divergence",
    "cross_entropy",
    "total_correlation",
    "transfer_entropy",
    "jensen_shannon_divergence",
    "renyi_entropy",
    "tsallis_entropy",
    "normalized_mutual_information",
    "information_coefficient",
    # Semantic information theory
    "information_content",
    "semantic_entropy",
    "semantic_similarity",
    "semantic_similarity_matrix",
    # Analysis functions
    "information_profile",
    "information_signature",
    "analyze_sequence_information",
    "compare_sequences_information",
    # Continuous information theory
    "differential_entropy",
    "mutual_information_continuous",
    "kl_divergence_continuous",
    "entropy_estimation",
    # Estimation methods
    "entropy_estimator",
    "mutual_information_estimator",
    "kl_divergence_estimator",
    "bias_correction",
    # Workflow functions
    "batch_entropy_analysis",
    "information_workflow",
    "compare_datasets",
    "information_report",
    # Integration functions
    "dna_integration",
    "rna_integration",
    "singlecell_integration",
    "multiomics_integration",
    "ml_integration",
]

