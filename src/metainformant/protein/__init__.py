"""Protein sequence and structure analysis module for METAINFORMANT."""
from __future__ import annotations

from . import database, domains, function, sequence, structure, visualization, workflow

# alphafold is genuinely optional (external API dependency)
try:
    from .structure import alphafold
except ImportError:
    alphafold = None  # type: ignore[assignment]

__all__ = [
    "database",
    "domains",
    "function",
    "sequence",
    "structure",
    "visualization",
    "workflow",
    "alphafold",
]
