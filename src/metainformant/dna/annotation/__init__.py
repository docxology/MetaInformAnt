"""DNA annotation submodule for gene prediction and functional annotation.

This submodule provides tools for gene prediction (ORF finding, regulatory element
detection, CpG island annotation, codon usage analysis, splice site prediction)
and functional annotation (variant classification, impact prediction, conservation
scoring, domain identification).
"""

from __future__ import annotations

from . import functional, gene_prediction

__all__ = [
    "gene_prediction",
    "functional",
]
