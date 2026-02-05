"""Protein sequence and structure analysis module for METAINFORMANT.

This module provides comprehensive tools for protein analysis, including
sequence processing, structure prediction, functional annotation,
protein-protein interactions, proteome analysis, and orchestrated workflows.
"""

from __future__ import annotations

# Import subpackages
from . import database
from . import sequence
from . import structure
from . import visualization

# Import modules from subpackages for backward compatibility
from .sequence import (
    alignment,
    proteomes,
    sequences,
)
from .structure import (
    analysis as structure_analysis,
    io as structure_io,
    general as structure_general,
)

# Optional imports with graceful fallbacks
try:
    from .structure import alphafold
except ImportError:
    alphafold = None  # type: ignore[assignment]

try:
    from .structure import contacts
except ImportError:
    contacts = None  # type: ignore[assignment]

try:
    from .structure import secondary
except ImportError:
    secondary = None  # type: ignore[assignment]

try:
    from .structure import pdb
except ImportError:
    pdb = None  # type: ignore[assignment]

try:
    from .database import interpro
except ImportError:
    interpro = None  # type: ignore[assignment]

try:
    from .database import uniprot
except ImportError:
    uniprot = None  # type: ignore[assignment]

# Visualization
from .visualization import general as visualization_general

# Orchestration
try:
    from . import orchestration
except ImportError:
    orchestration = None  # type: ignore[assignment]

__all__ = [
    # Subpackages
    "database",
    "sequence",
    "structure",
    "visualization",
    # Core protein analysis
    "sequences",
    "alignment",
    "structure_general",
    "structure_analysis",
    "structure_io",
    # Structure prediction and modeling
    "alphafold",
    "contacts",
    "secondary",
    # Functional annotation and databases
    "interpro",
    "uniprot",
    "pdb",
    # Proteome analysis
    "proteomes",
    # Visualization
    "visualization_general",
    # Orchestration
    "orchestration",
]
