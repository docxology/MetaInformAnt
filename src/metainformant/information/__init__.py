"""Information theory analysis module for METAINFORMANT.

This module provides comprehensive information-theoretic methods for biological data,
including syntactic information (Shannon entropy, mutual information),
semantic information (information content, semantic similarity),
continuous information theory, and analysis workflows.
"""

from __future__ import annotations

# Import all information theory submodules
from . import (
    analysis,
    continuous,
    estimation,
    integration,
    networks,
    semantic,
    syntactic,
    visualization,
    workflows,
)

# Direct imports of commonly used functions
from .syntactic import (
    shannon_entropy,
    shannon_entropy_from_counts,
    joint_entropy,
    conditional_entropy,
    mutual_information,
    conditional_mutual_information,
    kl_divergence,
    cross_entropy,
    jensen_shannon_divergence,
    total_correlation,
    transfer_entropy,
    renyi_entropy,
    tsallis_entropy,
    normalized_mutual_information,
    information_coefficient,
)

from .semantic import (
    information_content,
    information_content_from_annotations,
    semantic_entropy,
    semantic_similarity,
    semantic_similarity_matrix,
)

from .analysis import (
    information_profile,
    information_signature,
    analyze_sequence_information,
    compare_sequences_information,
)

from .continuous import (
    differential_entropy,
    mutual_information_continuous,
    kl_divergence_continuous,
    entropy_estimation,
)

from .estimation import (
    entropy_estimator,
    mutual_information_estimator,
    kl_divergence_estimator,
    bias_correction,
)

from .workflows import (
    batch_entropy_analysis,
    information_workflow,
    compare_datasets,
    information_report,
)

from .networks import (
    network_entropy,
    information_flow,
)

# Optional imports with graceful fallbacks
try:
    from . import integration
except ImportError:
    integration = None

try:
    from . import networks
except ImportError:
    networks = None

try:
    from . import visualization
except ImportError:
    visualization = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core information measures
    "syntactic",
    "semantic",
    "continuous",

    # Analysis and estimation
    "analysis",
    "estimation",
    "workflows",

    # Specialized applications
    "integration",
    "networks",
    "visualization",

    # Syntactic information functions
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

    # Semantic information functions
    "information_content",
    "information_content_from_annotations",
    "semantic_entropy",
    "semantic_similarity",
    "semantic_similarity_matrix",

    # Analysis functions
    "information_profile",
    "information_signature",
    "analyze_sequence_information",
    "compare_sequences_information",

    # Continuous information functions
    "differential_entropy",
    "mutual_information_continuous",
    "kl_divergence_continuous",
    "entropy_estimation",

    # Estimation functions
    "entropy_estimator",
    "mutual_information_estimator",
    "kl_divergence_estimator",
    "bias_correction",

    # Workflow functions
    "batch_entropy_analysis",
    "information_workflow",
    "compare_datasets",
    "information_report",

    # Network functions
    "network_entropy",
    "information_flow",
]



