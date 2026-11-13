"""DNA sequence analysis and genomics module for METAINFORMANT."""

from __future__ import annotations

from . import (
    alignment,
    annotation,
    expression,
    external,
    integration,
    io,
    phylogeny,
    population,
    sequence,
    variation,
)
from .alignment import msa

__all__ = [
    "alignment",
    "annotation",
    "expression",
    "external",
    "integration",
    "io",
    "phylogeny",
    "population",
    "sequence",
    "variation",
    "msa",
]
