"""Evolutionary and Population Genetics Simulation module for METAINFORMANT.

This module provides comprehensive tools for simulating evolutionary processes,
including agent-based modeling, population genetics simulation, molecular evolution,
and sequence generation.

Capabilities:
    - **Sequence Simulation**: Generate random DNA/protein sequences and apply mutations
    - **RNA-seq Simulation**: Simulate count data using negative binomial distributions
    - **Agent-based Models**: Create agents in grid-based environments for ecological modeling
    - **Population Genetics**: Forward and backward-time coalescent simulations
    - **Molecular Evolution**: Simulate substitution models and phylogenetic trees

Key Classes:
    - Agent: Individual agent for agent-based modeling
    - GridWorld: 2D grid environment for agent simulations

Key Functions:
    - generate_random_dna: Generate random DNA sequences
    - generate_random_protein: Generate random protein sequences
    - mutate_sequence: Apply mutations to sequences
    - simulate_counts_negative_binomial: Simulate RNA-seq count data

Submodules:
    - models.sequences: DNA/protein sequence simulation
    - models.rna: RNA-seq count simulation
    - models.popgen: Population genetics simulations
    - models.agents: Agent-based modeling
    - workflow: Simulation workflow management
    - visualization: Plotting simulation results

Example:
    >>> from metainformant.simulation import generate_random_dna, mutate_sequence
    >>> seq = generate_random_dna(length=100)
    >>> mutated = mutate_sequence(seq, mutation_rate=0.01)

    >>> from metainformant.simulation import simulate_counts_negative_binomial
    >>> counts = simulate_counts_negative_binomial(n_genes=1000, n_samples=10)

See Also:
    - docs/simulation/ for detailed simulation documentation
    - metainformant.math for population genetics theory
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
