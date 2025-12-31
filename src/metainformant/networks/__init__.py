"""Biological network analysis module for METAINFORMANT.

This module provides comprehensive tools for biological network analysis,
including graph construction, community detection, centrality measures,
pathway analysis, protein-protein interactions, and regulatory networks.
"""

from __future__ import annotations

# Import all network analysis submodules
from . import (
    community,
    graph,
    pathway,
    ppi,
    regulatory,
)

# Optional imports with graceful fallbacks
try:
    from . import regulatory
except ImportError:
    regulatory = None

try:
    from . import pathway
except ImportError:
    pathway = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core network functionality
    "graph",
    "community",

    # Biological network types
    "ppi",
    "pathway",
    "regulatory",
]
