"""Multi-omics data integration and analysis.

This module provides tools for integrating and analyzing multiple omics datasets:
- Data integration across genomics, transcriptomics, proteomics, metabolomics
- Multi-view learning and joint dimensionality reduction
- Cross-omics correlation and association analysis
- Pathway-level integration and enrichment
- Biomarker discovery across omics layers
"""

from .integration import (
    MultiOmicsData,
    canonical_correlation,
    from_dna_variants,
    from_protein_abundance,
    from_rna_expression,
    integrate_omics_data,
    joint_nmf,
    joint_pca,
)

__all__ = [
    # Data integration
    "MultiOmicsData",
    "integrate_omics_data",
    "joint_pca",
    "joint_nmf",
    "canonical_correlation",
    # Integration helpers
    "from_dna_variants",
    "from_rna_expression",
    "from_protein_abundance",
]
