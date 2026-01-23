"""Evolutionary and Population Genetics Simulation module for METAINFORMANT.

This module provides tools for simulating evolutionary processes, including
agent-based modeling, population genetics simulation (forward and backward time),
molecular evolution, and sequence simulation.
"""

from __future__ import annotations

# Import subpackages
from . import models
from . import visualization
from . import workflow

# Import modules from subpackages for backward compatibility
from .models import (
    agents,
    popgen,
    rna,
    sequences,
)
from .models.agents import Agent, GridWorld
from .models.sequences import (
    generate_random_dna,
    generate_random_protein,
    mutate_sequence,
)
from .models.rna import simulate_counts_negative_binomial
from .workflow import workflow as workflow_module
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .models import popgen
except ImportError:
    popgen = None

try:
    from .models import agents
except ImportError:
    agents = None

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
    # Visualization
    "visualization",
    # Direct exports
    "Agent",
    "GridWorld",
    "generate_random_dna",
    "generate_random_protein",
    "mutate_sequence",
    "simulate_counts_negative_binomial",
]
