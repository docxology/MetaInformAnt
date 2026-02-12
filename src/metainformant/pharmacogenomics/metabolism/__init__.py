"""Metabolizer status prediction subpackage for pharmacogenomics.

Provides CPIC-style metabolizer phenotype prediction from genotype data,
activity score computation, and dose adjustment recommendations."""
from __future__ import annotations

from . import metabolizer_status

__all__ = ['metabolizer_status']
