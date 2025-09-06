"""Multi-omics data integration and analysis.

This module provides tools for integrating and analyzing multiple omics datasets:
- Data integration across genomics, transcriptomics, proteomics, metabolomics
- Multi-view learning and joint dimensionality reduction
- Cross-omics correlation and association analysis
- Pathway-level integration and enrichment
- Biomarker discovery across omics layers
"""

from .association import cross_omics_correlation, multi_omics_gwas, omics_qtl_analysis, trans_omics_network
from .biomarkers import cross_omics_signatures, integrated_feature_selection, multi_omics_biomarker_discovery
from .clustering import consensus_clustering_omics, multi_view_clustering, omics_subtype_discovery
from .integration import MultiOmicsData, canonical_correlation, integrate_omics_data, joint_nmf, joint_pca
from .pathway import multi_omics_pathway_analysis, omics_pathway_crosstalk, pathway_activity_inference
from .visualization import circos_plot_omics, plot_multi_omics_network, plot_omics_correlation, plot_pathway_heatmap

__all__ = [
    # Data integration
    "MultiOmicsData",
    "integrate_omics_data",
    "joint_pca",
    "joint_nmf",
    "canonical_correlation",
    # Association analysis
    "cross_omics_correlation",
    "trans_omics_network",
    "omics_qtl_analysis",
    "multi_omics_gwas",
    # Pathway analysis
    "pathway_activity_inference",
    "multi_omics_pathway_analysis",
    "omics_pathway_crosstalk",
    # Biomarker discovery
    "multi_omics_biomarker_discovery",
    "cross_omics_signatures",
    "integrated_feature_selection",
    # Clustering
    "multi_view_clustering",
    "consensus_clustering_omics",
    "omics_subtype_discovery",
    # Visualization
    "plot_omics_correlation",
    "plot_pathway_heatmap",
    "plot_multi_omics_network",
    "circos_plot_omics",
]
