"""Spatial-scRNA-seq integration module for METAINFORMANT.

Provides methods for mapping single-cell RNA-seq data to spatial transcriptomics,
including label transfer, correlation-based mapping, and gene imputation."""
from __future__ import annotations

from . import scrna_mapping

__all__ = ['scrna_mapping']
