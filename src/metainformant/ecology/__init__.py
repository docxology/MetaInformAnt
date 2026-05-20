"""Ecology and community analysis module for METAINFORMANT.

This module provides tools for ecological analysis, including community
structure, diversity metrics, ordination, indicator species, functional
ecology, macroecology, and ecological visualization.
"""

from __future__ import annotations

from . import analysis, phylogenetic, traits, visualization
from .analysis.community import beta_diversity, calculate_diversity, species_richness

__all__ = [
    "analysis",
    "phylogenetic",
    "traits",
    "visualization",
    "calculate_diversity",
    "species_richness",
    "beta_diversity",
]
