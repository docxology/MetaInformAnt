"""Spatial transcriptomics module for METAINFORMANT.

Provides comprehensive tools for spatial transcriptomics data analysis,
supporting major platforms (10x Visium, MERFISH, 10x Xenium) with spatial
statistics, clustering, deconvolution, integration, and visualization.

Key Features:
    - Multi-platform I/O: Load data from Visium, MERFISH, and Xenium.
    - Spatial clustering: Leiden/Louvain with spatial graph construction.
    - Cell type deconvolution: NNLS and NMF-based spot deconvolution.
    - Spatial statistics: Moran's I, Geary's C, Ripley's K, variograms.
    - Neighborhood analysis: Co-localization, niches, ligand-receptor signaling.
    - scRNA-seq integration: Label transfer and gene imputation.
    - Visualization: Spatial scatter, tissue overlay, LISA maps.

Config prefix: SPATIAL_

Example usage::

    from metainformant.spatial import io, analysis, integration, visualization

    # Load Visium data
    dataset = io.load_visium("path/to/spaceranger_output/")

    # Spatial clustering
    result = analysis.spatial_cluster(
        dataset.expression, dataset.coordinates,
        method="leiden", spatial_weight=0.5,
    )

    # Moran's I for a gene
    weights = analysis.spatial_weights_matrix(dataset.coordinates, method="knn", k=6)
    gene_idx = dataset.gene_names.index("CD3E")
    morans = analysis.morans_i(dataset.expression[:, gene_idx].toarray().flatten(), weights)

    # Visualize
    visualization.plot_spatial_scatter(
        dataset.coordinates, result.labels,
        "output/spatial/clusters.png",
    )
"""

from __future__ import annotations

# Import subpackages
from . import analysis, communication, deconvolution, integration, io, visualization

# IO - Data loading
from .io import (
    # Visium
    create_spatial_dataset,
    filter_tissue_spots,
    load_visium,
    read_spatial_image,
    read_tissue_positions,
    # MERFISH
    aggregate_to_cells,
    load_merfish,
    load_transcript_spots,
    parse_cell_metadata,
    # Xenium
    load_cell_boundaries,
    load_xenium,
    read_cell_features,
    read_transcripts,
)

# Analysis - Spatial statistics and clustering
from .analysis import (
    # Clustering
    build_spatial_graph,
    leiden_clustering,
    louvain_clustering,
    spatial_cluster,
    spatial_domains,
    # Deconvolution
    create_reference_profiles,
    deconvolve_spots,
    enrichment_score,
    estimate_cell_fractions,
    nnls_deconvolution,
    # Neighborhood
    compute_interaction_matrix,
    ligand_receptor_spatial,
    neighborhood_enrichment,
    niche_detection,
    ripley_k,
    # Autocorrelation
    gearys_c,
    getis_ord_g,
    local_morans_i,
    morans_i,
    spatial_variogram,
    spatial_weights_matrix,
)

# Integration
from .integration import (
    anchor_based_transfer,
    correlation_mapping,
    impute_spatial_genes,
    map_scrna_to_spatial,
)

# Visualization
from .visualization import (
    plot_cell_type_map,
    plot_deconvolution_pie,
    plot_gene_expression_map,
    plot_interaction_heatmap,
    plot_neighborhood_graph,
    plot_spatial_autocorrelation,
    plot_spatial_scatter,
    plot_tissue_overlay,
)

# Data classes re-exported for convenience
from .io.visium import SpatialDataset, TissuePosition
from .io.merfish import CellMetadata, MERFISHDataset, TranscriptSpot
from .io.xenium import CellBoundary, XeniumDataset, XeniumTranscript
from .analysis.clustering import SpatialClusterResult
from .analysis.deconvolution import DeconvolutionResult
from .analysis.neighborhood import (
    InteractionResult,
    NeighborhoodEnrichmentResult,
    NicheResult,
    RipleyKResult,
)
from .analysis.autocorrelation import (
    GearyCResult,
    GetisOrdResult,
    LocalMoransResult,
    MoransIResult,
    VariogramResult,
)
from .integration.scrna_mapping import ImputationResult, MappingResult

# Deconvolution submodule
from .deconvolution import (
    build_reference_profiles as deconv_build_reference_profiles,
    deconvolve_spots as deconv_deconvolve_spots,
    niche_identification,
    spatial_cell_type_mapping,
    validate_deconvolution,
)

# Communication submodule
from .communication import (
    build_communication_network,
    communication_pattern_analysis,
    compute_ligand_receptor_interactions,
    default_lr_database,
    spatial_interaction_score,
)

__all__ = [
    # Subpackages
    "analysis",
    "communication",
    "deconvolution",
    "integration",
    "io",
    "visualization",
    # IO functions
    "create_spatial_dataset",
    "filter_tissue_spots",
    "load_visium",
    "read_spatial_image",
    "read_tissue_positions",
    "aggregate_to_cells",
    "load_merfish",
    "load_transcript_spots",
    "parse_cell_metadata",
    "load_cell_boundaries",
    "load_xenium",
    "read_cell_features",
    "read_transcripts",
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
    # Integration
    "anchor_based_transfer",
    "correlation_mapping",
    "impute_spatial_genes",
    "map_scrna_to_spatial",
    # Visualization
    "plot_cell_type_map",
    "plot_deconvolution_pie",
    "plot_gene_expression_map",
    "plot_interaction_heatmap",
    "plot_neighborhood_graph",
    "plot_spatial_autocorrelation",
    "plot_spatial_scatter",
    "plot_tissue_overlay",
    # Deconvolution submodule
    "deconv_build_reference_profiles",
    "deconv_deconvolve_spots",
    "niche_identification",
    "spatial_cell_type_mapping",
    "validate_deconvolution",
    # Communication submodule
    "build_communication_network",
    "communication_pattern_analysis",
    "compute_ligand_receptor_interactions",
    "default_lr_database",
    "spatial_interaction_score",
    # Data classes
    "CellBoundary",
    "CellMetadata",
    "DeconvolutionResult",
    "GearyCResult",
    "GetisOrdResult",
    "ImputationResult",
    "InteractionResult",
    "LocalMoransResult",
    "MERFISHDataset",
    "MappingResult",
    "MoransIResult",
    "NeighborhoodEnrichmentResult",
    "NicheResult",
    "RipleyKResult",
    "SpatialClusterResult",
    "SpatialDataset",
    "TissuePosition",
    "TranscriptSpot",
    "VariogramResult",
    "XeniumDataset",
    "XeniumTranscript",
]
