"""GWAS analysis modules."""
from __future__ import annotations

from . import annotation, association, calling, correction, heritability, ld_pruning, mixed_model, quality, structure, summary_stats, utils
from .utils import compute_r_squared

__all__ = [
    'annotation', 'association', 'calling', 'correction', 'heritability',
    'ld_pruning', 'mixed_model', 'quality', 'structure', 'summary_stats', 'utils',
    'compute_r_squared',
]
