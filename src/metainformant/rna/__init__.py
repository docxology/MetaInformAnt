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
    metadata_filter,
    monitoring,
    orchestration,
    pipeline,
    progress_tracker,
    workflow,
)

# Import steps submodule
from . import steps

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
    "steps",

    # Configuration and setup
    "configs",

    # Utilities
    "cleanup",
    "deps",
    "metadata_filter",

    # Direct exports
    "RNAPipelineConfig",
    "AmalgkitRunLayout",
    "ProgressTracker",
    "summarize_curate_tables",
    "AmalgkitWorkflowConfig",
]

