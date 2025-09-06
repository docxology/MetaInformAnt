"""Simulation and synthetic data generation utilities.

This subpackage provides modular generators across domains and toy
agent-based simulations to support benchmarking and method development.
"""

from .agents import Agent, GridWorld
from .rna import simulate_counts_negative_binomial
from .sequences import generate_random_dna, generate_random_protein, mutate_sequence

__all__ = [
    "generate_random_dna",
    "mutate_sequence",
    "generate_random_protein",
    "simulate_counts_negative_binomial",
    "GridWorld",
    "Agent",
]
