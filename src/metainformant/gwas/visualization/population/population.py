"""Population structure visualization for GWAS.

This module provides plots for visualizing population stratification and structure.

This module re-exports all public symbols from :mod:`population_pca` and
:mod:`population_admixture` for backward compatibility.
"""

from __future__ import annotations

from metainformant.gwas.visualization.population.population_pca import (
    _load_pca_from_file,
    pca_3d,
    pca_multi_panel,
    pca_plot,
    pca_scree_plot,
)
from metainformant.gwas.visualization.population.population_admixture import (
    _load_admixture_from_file,
    _load_kinship_from_file,
    admixture_plot,
    kinship_clustermap,
    kinship_dendrogram,
    kinship_heatmap,
)

__all__ = [
    # PCA functions
    "_load_pca_from_file",
    "pca_plot",
    "pca_scree_plot",
    "pca_multi_panel",
    "pca_3d",
    # Admixture/kinship functions
    "_load_kinship_from_file",
    "kinship_heatmap",
    "_load_admixture_from_file",
    "admixture_plot",
    "kinship_dendrogram",
    "kinship_clustermap",
]
