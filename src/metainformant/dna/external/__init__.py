"""External database interfaces (NCBI, Entrez, genome downloads)."""

from __future__ import annotations

from . import entrez, genomes, ncbi, ncbi_datasets

__all__ = [
    "entrez",
    "genomes",
    "ncbi",
    "ncbi_datasets",
]
