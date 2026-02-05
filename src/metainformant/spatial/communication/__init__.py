"""Cell-cell communication analysis for spatial transcriptomics.

Provides ligand-receptor interaction scoring, spatial communication networks,
and communication pattern discovery using spatial proximity information.
"""

from __future__ import annotations

from . import cell_communication

from .cell_communication import (
    build_communication_network,
    communication_pattern_analysis,
    compute_ligand_receptor_interactions,
    default_lr_database,
    spatial_interaction_score,
)

__all__ = [
    "cell_communication",
    "build_communication_network",
    "communication_pattern_analysis",
    "compute_ligand_receptor_interactions",
    "default_lr_database",
    "spatial_interaction_score",
]
