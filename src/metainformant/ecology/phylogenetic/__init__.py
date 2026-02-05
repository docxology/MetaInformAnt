"""Phylogenetic ecology subpackage.

Provides phylogenetic diversity metrics (Faith's PD, UniFrac),
phylogenetic community structure (NRI/NTI), phylogenetic signal
tests, and simple tree construction.
"""

from __future__ import annotations

from .diversity import (
    build_simple_tree,
    compute_unifrac,
    faiths_pd,
    nri_nti,
    phylogenetic_beta_diversity,
    phylogenetic_signal,
)

__all__ = [
    "faiths_pd",
    "phylogenetic_beta_diversity",
    "compute_unifrac",
    "nri_nti",
    "phylogenetic_signal",
    "build_simple_tree",
]
