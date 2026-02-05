"""Spatial data I/O module for METAINFORMANT.

Provides loaders for spatial transcriptomics platforms including
10x Visium, MERFISH, and 10x Xenium.
"""

from __future__ import annotations

from . import merfish, visium, xenium

# Visium
from .visium import (
    create_spatial_dataset,
    filter_tissue_spots,
    load_visium,
    read_spatial_image,
    read_tissue_positions,
)

# MERFISH
from .merfish import (
    aggregate_to_cells,
    load_merfish,
    load_transcript_spots,
    parse_cell_metadata,
)

# Xenium
from .xenium import (
    load_cell_boundaries,
    load_xenium,
    read_cell_features,
    read_transcripts,
)

__all__ = [
    # Submodules
    "merfish",
    "visium",
    "xenium",
    # Visium
    "create_spatial_dataset",
    "filter_tissue_spots",
    "load_visium",
    "read_spatial_image",
    "read_tissue_positions",
    # MERFISH
    "aggregate_to_cells",
    "load_merfish",
    "load_transcript_spots",
    "parse_cell_metadata",
    # Xenium
    "load_cell_boundaries",
    "load_xenium",
    "read_cell_features",
    "read_transcripts",
]
