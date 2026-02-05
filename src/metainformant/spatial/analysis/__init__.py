"""Spatial analysis module for METAINFORMANT.

Provides spatial statistics, clustering, deconvolution, and neighborhood
analysis algorithms for spatial transcriptomics data.
"""

from __future__ import annotations

from . import autocorrelation, clustering, deconvolution, neighborhood

# Clustering
from .clustering import (
    build_spatial_graph,
    leiden_clustering,
    louvain_clustering,
    spatial_cluster,
    spatial_domains,
)

# Deconvolution
from .deconvolution import (
    create_reference_profiles,
    deconvolve_spots,
    enrichment_score,
    estimate_cell_fractions,
    nnls_deconvolution,
)

# Neighborhood
from .neighborhood import (
    compute_interaction_matrix,
    ligand_receptor_spatial,
    neighborhood_enrichment,
    niche_detection,
    ripley_k,
)

# Autocorrelation
from .autocorrelation import (
    gearys_c,
    getis_ord_g,
    local_morans_i,
    morans_i,
    spatial_variogram,
    spatial_weights_matrix,
)

__all__ = [
    # Submodules
    "autocorrelation",
    "clustering",
    "deconvolution",
    "neighborhood",
    # Clustering
    "build_spatial_graph",
    "leiden_clustering",
    "louvain_clustering",
    "spatial_cluster",
    "spatial_domains",
    # Deconvolution
    "create_reference_profiles",
    "deconvolve_spots",
    "enrichment_score",
    "estimate_cell_fractions",
    "nnls_deconvolution",
    # Neighborhood
    "compute_interaction_matrix",
    "ligand_receptor_spatial",
    "neighborhood_enrichment",
    "niche_detection",
    "ripley_k",
    # Autocorrelation
    "gearys_c",
    "getis_ord_g",
    "local_morans_i",
    "morans_i",
    "spatial_variogram",
    "spatial_weights_matrix",
]
