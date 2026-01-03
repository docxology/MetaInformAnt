"""Protein sequence and structure analysis module for METAINFORMANT.

This module provides comprehensive tools for protein analysis, including
sequence processing, structure prediction, functional annotation,
protein-protein interactions, and proteome analysis.
"""

from __future__ import annotations

# Import all protein analysis submodules
from . import (
    alignment,
    alphafold,
    contacts,
    interpro,
    pdb,
    proteomes,
    secondary,
    sequences,
    structure,
    structure_analysis,
    structure_io,
    uniprot,
    visualization,
)

# Optional imports with graceful fallbacks
try:
    from . import alphafold
except ImportError:
    alphafold = None

try:
    from . import interpro
except ImportError:
    interpro = None

try:
    from . import uniprot
except ImportError:
    uniprot = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core protein analysis
    "sequences",
    "alignment",
    "structure",
    "structure_analysis",
    "structure_io",

    # Structure prediction and modeling
    "alphafold",
    "contacts",
    "secondary",

    # Functional annotation and databases
    "interpro",
    "uniprot",
    "map_ids_uniprot",
    "pdb",

    # Proteome analysis
    "proteomes",

    # Visualization
    "visualization",
]



