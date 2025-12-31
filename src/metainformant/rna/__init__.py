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
    monitoring,
    orchestration,
    pipeline,
    progress_tracker,
    protein_integration,
    workflow,
)

# Import steps submodule
from . import steps

# Optional imports with graceful fallbacks
try:
    from . import protein_integration
except ImportError:
    protein_integration = None

try:
    from . import discovery
except ImportError:
    discovery = None

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
    "environment",
    "genome_prep",

    # Utilities
    "cleanup",
    "deps",

    # Optional/advanced modules
    "protein_integration",
    "discovery",
]
