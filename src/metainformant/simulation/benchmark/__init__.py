"""Benchmark dataset generation subpackage for simulation.

Provides synthetic benchmark dataset generators for classification,
regression, clustering, differential expression, GWAS, and method
comparison/evaluation."""
from __future__ import annotations

from . import generators

__all__ = ['generators']
