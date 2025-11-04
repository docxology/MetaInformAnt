"""Information theory methods for biological data analysis.

This module provides comprehensive information theory tools for analyzing
biological data, including both syntactic (Shannon entropy, mutual information)
and semantic (information content, semantic similarity) information measures.

Modules:
- syntactic: Shannon entropy, mutual information, KL divergence, etc.
- semantic: Information content, semantic similarity, semantic entropy
- analysis: High-level analysis functions and workflows
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
    joint_entropy,
    kl_divergence,
    mutual_information,
    shannon_entropy,
    total_correlation,
    transfer_entropy,
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
]

