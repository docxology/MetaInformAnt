"""Information integration subpackage.

Provides cross-platform integration functions for combining information-theoretic
analysis with DNA, RNA, single-cell, multi-omics, ML, and network data.
"""

from __future__ import annotations

from .integration import (
    dna_integration,
    rna_integration,
    singlecell_integration,
    multiomics_integration,
    ml_integration,
)
from .networks import (
    network_entropy,
    information_flow,
    network_information_centrality,
    network_motif_information,
    information_graph_distance,
    information_community_detection,
)

__all__ = [
    "dna_integration",
    "rna_integration",
    "singlecell_integration",
    "multiomics_integration",
    "ml_integration",
    "network_entropy",
    "information_flow",
    "network_information_centrality",
    "network_motif_information",
    "information_graph_distance",
    "information_community_detection",
]
