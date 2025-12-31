"""Biological simulation and synthetic data generation module for METAINFORMANT.

This module provides tools for generating synthetic biological data and running
agent-based simulations, including sequence evolution, population dynamics,
RNA-seq data simulation, and workflow testing data.
"""

from __future__ import annotations

# Import all simulation submodules
from . import (
    agents,
    popgen,
    rna,
    sequences,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import agents
except ImportError:
    agents = None

try:
    from . import workflow
except ImportError:
    workflow = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Sequence and molecular simulations
    "sequences",
    "rna",

    # Population and evolutionary simulations
    "popgen",
    "agents",

    # Workflow and testing
    "workflow",
]


