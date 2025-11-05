"""Simulation and synthetic data generation utilities.

This subpackage provides modular generators across domains and toy
agent-based simulations to support benchmarking and method development.
"""

from .agents import Agent, GridWorld
from .popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)
from .rna import simulate_counts_negative_binomial
from .sequences import generate_random_dna, generate_random_protein, mutate_sequence
from .workflow import (
    run_agent_simulation_workflow,
    run_popgen_simulation_workflow,
    run_sequence_simulation_workflow,
)

__all__ = [
    "generate_random_dna",
    "mutate_sequence",
    "generate_random_protein",
    "simulate_counts_negative_binomial",
    "GridWorld",
    "Agent",
    # Population genetics
    "generate_population_sequences",
    "generate_two_populations",
    "generate_genotype_matrix",
    "simulate_bottleneck_population",
    "simulate_population_expansion",
    "generate_site_frequency_spectrum",
    "generate_linkage_disequilibrium_data",
    # Workflows
    "run_sequence_simulation_workflow",
    "run_agent_simulation_workflow",
    "run_popgen_simulation_workflow",
]
