"""Gene annotation and prediction utilities.

This module re-exports all public symbols from the split submodules
``gene_finding`` and ``gene_annotation`` so that existing import paths
continue to work unchanged.
"""

from __future__ import annotations

# Re-export everything from the split modules
from .gene_finding import (
    GENETIC_CODES,
    IUPAC_MAP,
    REGULATORY_MOTIFS,
    START_CODONS,
    STOP_CODONS,
    predict_orfs,
)
from .gene_annotation import (
    annotate_coding_regions,
    annotate_cpg_islands,
    compute_codon_usage,
    find_regulatory_elements,
    find_splice_sites,
    mask_repeats,
)

__all__ = [
    "GENETIC_CODES",
    "IUPAC_MAP",
    "REGULATORY_MOTIFS",
    "START_CODONS",
    "STOP_CODONS",
    "predict_orfs",
    "annotate_coding_regions",
    "annotate_cpg_islands",
    "compute_codon_usage",
    "find_regulatory_elements",
    "find_splice_sites",
    "mask_repeats",
]
