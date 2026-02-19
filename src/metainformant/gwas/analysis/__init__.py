"""GWAS analysis modules."""
from __future__ import annotations

from . import annotation, association, benchmarking, calling, correction, eqtl, heritability, ld_pruning, mixed_model, quality, structure, summary_stats, utils
from .utils import compute_r_squared

__all__ = [
    'annotation', 'association', 'benchmarking', 'calling', 'correction', 'eqtl',
    'heritability', 'ld_pruning', 'mixed_model', 'quality', 'structure',
    'summary_stats', 'utils', 'compute_r_squared',
]

