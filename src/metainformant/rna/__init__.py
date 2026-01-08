"""RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

This module provides comprehensive tools for RNA-seq analysis, including
amalgkit workflow orchestration, expression quantification, quality control,
differential expression analysis, and multi-species comparative studies.
"""

from __future__ import annotations

# Import subpackages
from . import core        # Base configs and utils
from . import analysis    # Analysis modules
from . import amalgkit    # Amalgkit integration
from . import engine      # Orchestration engine (depends on others)

# Import modules from subpackages for backward compatibility
from .core import (
    cleanup,
    configs,
    deps,
    environment,
)
from .engine import (
    discovery,
    monitoring,
    orchestration,
    pipeline,
    progress_tracker,
    workflow,
)
from .amalgkit import (
    amalgkit as amalgkit_module,
    genome_prep,
    metadata_filter,
)
from .analysis import (
    protein_integration,
    validation,
)

# Direct imports of commonly used classes and functions
from .core.configs import RNAPipelineConfig, AmalgkitRunLayout
from .engine.progress_tracker import ProgressTracker
from .engine.pipeline import summarize_curate_tables
from .engine.workflow import AmalgkitWorkflowConfig

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "amalgkit",
    "analysis",
    "core",
    "engine",

    # Core workflow management
    "workflow",
    "orchestration",
    "pipeline",
    "monitoring",
    "progress_tracker",

    # Amalgkit integration
    "amalgkit_module",

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

