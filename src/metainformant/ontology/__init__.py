"""Gene ontology and functional annotation module for METAINFORMANT.

This module provides tools for working with gene ontologies, functional annotation,
semantic similarity analysis, and biological knowledge representation.
"""

from __future__ import annotations

# Import all ontology submodules
from . import (
    go,
    obo,
    pathway_enrichment,
    query,
    serialize,
    types,
    visualization,
)

# Pathway enrichment imports
from .pathway_enrichment.enrichment import (
    compare_enrichments,
    compute_enrichment_score,
    gsea,
    over_representation_analysis,
    pathway_network,
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
    "visualization",
    # Pathway enrichment
    "pathway_enrichment",
    "over_representation_analysis",
    "gsea",
    "compute_enrichment_score",
    "pathway_network",
    "compare_enrichments",
]
