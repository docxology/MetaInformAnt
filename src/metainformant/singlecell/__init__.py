"""Single-cell analysis module for METAINFORMANT.

This module provides tools for single-cell data analysis, including
dimensionality reduction, clustering, trajectory inference, multi-sample
integration, cell type annotation, differential expression, and RNA
velocity estimation.
"""

from __future__ import annotations

from . import analysis, celltyping, data, differential, doublet, io, velocity, visualization
from .io import (
    fetch_atlas_datasets,
    filter_regulons_post_scenic,
    run_differential_expression,
    run_go_enrichment,
    run_pseudotime_trajectory,
    run_salmon_alevin,
    run_scenic_aucell,
    run_scenic_cistarget,
    run_scenic_grnboost2,
    run_seurat_integration,
    run_shinycell_export,
)

__all__ = [
    "analysis",
    "celltyping",
    "data",
    "differential",
    "doublet",
    "io",
    "velocity",
    "visualization",
    "fetch_atlas_datasets",
    "run_salmon_alevin",
    "run_seurat_integration",
    "run_differential_expression",
    "run_scenic_grnboost2",
    "run_scenic_cistarget",
    "run_scenic_aucell",
    "filter_regulons_post_scenic",
    "run_pseudotime_trajectory",
    "run_go_enrichment",
    "run_shinycell_export",
]
