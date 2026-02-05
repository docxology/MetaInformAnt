"""Spatial deconvolution submodule for advanced cell type deconvolution.

Extends the core spatial analysis deconvolution with reference profile building,
spatial cell type mapping, deconvolution validation, and tissue niche identification.
"""

from __future__ import annotations

from . import spatial_deconvolution
from .spatial_deconvolution import (
    build_reference_profiles,
    deconvolve_spots,
    niche_identification,
    spatial_cell_type_mapping,
    validate_deconvolution,
)

__all__ = [
    "spatial_deconvolution",
    "build_reference_profiles",
    "deconvolve_spots",
    "niche_identification",
    "spatial_cell_type_mapping",
    "validate_deconvolution",
]
