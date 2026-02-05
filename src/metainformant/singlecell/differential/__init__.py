"""Differential expression analysis for single-cell data.

Provides statistical testing for differential gene expression between
cell groups, pseudobulk aggregation, volcano plot preparation, and
gene set scoring.
"""

from __future__ import annotations

from . import expression

__all__ = [
    "expression",
]
