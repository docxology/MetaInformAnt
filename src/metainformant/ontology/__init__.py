"""Gene ontology and functional annotation module for METAINFORMANT.

This module provides tools for working with gene ontologies, functional annotation,
semantic similarity analysis, and biological knowledge representation.
"""

from __future__ import annotations

# Import all ontology submodules
from . import (
    go,
    obo,
    query,
    serialize,
    types,
)

# Optional imports with graceful fallbacks
try:
    from . import serialize
except ImportError:
    serialize = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core ontology functionality
    "go",
    "obo",
    "query",

    # Data structures and utilities
    "types",
    "serialize",
]
