"""Functional metagenomics analysis.

Provides gene annotation, ORF prediction, gene family classification,
and metabolic pathway reconstruction from metagenomic sequences.
"""

from __future__ import annotations

from . import annotation, pathways

__all__ = [
    "annotation",
    "pathways",
]
