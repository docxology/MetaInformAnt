"""Multi-omics integration module for METAINFORMANT.

This module provides comprehensive tools for integrating multiple omics data layers,
enabling cross-platform analysis of genomics, transcriptomics, epigenomics, and
proteomics data.

Capabilities:
    - **Data Integration**: Combine multi-omic datasets with matched samples
    - **Feature Mapping**: Map features across omic layers (genes, proteins, etc.)
    - **Dimensionality Reduction**: Multi-omic PCA, CCA, and factor analysis
    - **Pathway Integration**: Integrate pathway-level signals across omics
    - **Network Integration**: Build multi-layer biological networks
    - **Integration Methods**: Joint NMF, MOFA, tensor decomposition, SNF, CCA
    - **Clustering**: Multi-omic clustering, consensus clustering, spectral methods
    - **Pathway Analysis**: Multi-omic enrichment, active modules, topology analysis
    - **Survival Analysis**: Cox PH, Kaplan-Meier, log-rank, Lasso-Cox models

Submodules:
    - analysis.integration: Core integration algorithms and methods
    - visualization: Multi-omic visualization tools
    - methods: Matrix factorization, clustering, and network fusion
    - pathways: Multi-omic pathway enrichment and concordance
    - survival: Cox regression, Kaplan-Meier, risk stratification

Integration Methods:
    - Early integration: Concatenate features from all omics
    - Late integration: Combine predictions from omic-specific models
    - Joint integration: Simultaneous modeling of multiple omics

Typical Workflow:
    1. Load data from multiple omic layers
    2. Normalize and preprocess each layer
    3. Match samples across layers
    4. Apply integration method (MOFA, CCA, etc.)
    5. Interpret integrated factors/components
    6. Visualize integrated results

Example:
    >>> from metainformant.multiomics import analysis
    >>> # Load transcriptomics and proteomics data
    >>> rna_data = load_rna_data("expression.csv")
    >>> protein_data = load_protein_data("proteome.csv")
    >>> # Integrate using available methods
    >>> integrated = analysis.integration.integrate_omics(
    ...     {"rna": rna_data, "protein": protein_data},
    ...     method="early"
    ... )

Dependencies:
    - numpy, pandas for data handling
    - scipy for statistical methods
    - Optional: scikit-learn for ML integration methods

See Also:
    - metainformant.rna for RNA-seq analysis
    - metainformant.protein for proteomics analysis
    - metainformant.networks for network-based integration
"""

from __future__ import annotations

# Import subpackages
from . import analysis, visualization

# Import modules from subpackages for backward compatibility
from .analysis import (
    integration,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .analysis import integration
except ImportError:
    integration = None

try:
    from . import methods
except ImportError:
    methods = None  # type: ignore[assignment]

try:
    from . import pathways
except ImportError:
    pathways = None  # type: ignore[assignment]

try:
    from . import survival
except ImportError:
    survival = None  # type: ignore[assignment]

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "visualization",
    "methods",
    "pathways",
    "survival",
]
