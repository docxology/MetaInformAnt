"""Spatial-scRNA-seq integration module for METAINFORMANT.

Provides methods for mapping single-cell RNA-seq data to spatial transcriptomics,
including label transfer, correlation-based mapping, and gene imputation.
"""

from __future__ import annotations

from . import scrna_mapping
from .scrna_mapping import (
    anchor_based_transfer,
    correlation_mapping,
    impute_spatial_genes,
    map_scrna_to_spatial,
)

__all__ = [
    "scrna_mapping",
    "anchor_based_transfer",
    "correlation_mapping",
    "impute_spatial_genes",
    "map_scrna_to_spatial",
]
