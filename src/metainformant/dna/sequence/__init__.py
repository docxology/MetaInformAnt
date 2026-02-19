"""DNA sequence analysis sub-package (composition, k-mers, motifs, restriction)."""
from __future__ import annotations

from . import composition, consensus, core, kmer, motifs, restriction

__all__ = [
    "composition",
    "consensus",
    "core",
    "kmer",
    "motifs",
    "restriction",
]
