"""Community diversity metrics for metagenomic data.

Provides alpha diversity (Shannon, Simpson, Chao1, etc.), beta diversity
(Bray-Curtis, Jaccard, Aitchison), rarefaction analysis, PERMANOVA
testing, and ordination methods (PCoA, NMDS).
"""

from __future__ import annotations

from . import metrics

__all__ = [
    "metrics",
]
