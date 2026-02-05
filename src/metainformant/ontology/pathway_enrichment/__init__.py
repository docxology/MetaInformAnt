"""Pathway enrichment analysis subpackage.

Provides over-representation analysis (ORA), Gene Set Enrichment
Analysis (GSEA), pathway similarity networks, and cross-condition
enrichment comparison.
"""

from __future__ import annotations

from .enrichment import (
    compare_enrichments,
    compute_enrichment_score,
    gsea,
    over_representation_analysis,
    pathway_network,
)

__all__ = [
    "over_representation_analysis",
    "gsea",
    "compute_enrichment_score",
    "pathway_network",
    "compare_enrichments",
]
