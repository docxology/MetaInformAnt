"""RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

This module provides comprehensive tools for RNA-seq analysis, including
amalgkit workflow orchestration, expression quantification, quality control,
differential expression analysis, and multi-species comparative studies.
"""

from __future__ import annotations

# Import all RNA analysis submodules
from . import (
    amalgkit,
    cleanup,
    configs,
    deps,
    discovery,
    environment,
    genome_prep,
    metadata_filter,
    monitoring,
    orchestration,
    pipeline,
    progress_tracker,
    protein_integration,
    validation,
    workflow,
)

# Direct imports of commonly used classes and functions
from .configs import RNAPipelineConfig, AmalgkitRunLayout
from .progress_tracker import ProgressTracker
from .pipeline import summarize_curate_tables
from .workflow import AmalgkitWorkflowConfig

# All modules imported above are available

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core workflow management
    "workflow",
    "orchestration",
    "pipeline",
    "monitoring",
    "progress_tracker",

    # Amalgkit integration
    "amalgkit",

    # Configuration and setup
    "configs",

    # Environment and discovery
    "discovery",
    "environment",
    "genome_prep",
    "protein_integration",

    # Utilities
    "cleanup",
    "deps",
    "metadata_filter",
    "validation",

    # Direct exports
    "RNAPipelineConfig",
    "AmalgkitRunLayout",
    "ProgressTracker",
    "summarize_curate_tables",
    "AmalgkitWorkflowConfig",
]

