"""Information theory visualization integration.

This module provides unified access to information theory visualization functions
from the information module, integrating them into the central visualization package.
"""

from __future__ import annotations

# Re-export information theory visualization functions
try:
    from ..information.visualization import (
        plot_entropy_distribution,
        plot_entropy_landscape,
        plot_information_network,
        plot_information_profile,
        plot_mi_network,
        plot_mutual_information_matrix,
        plot_renyi_spectrum,
        plot_semantic_similarity_network,
    )

    INFORMATION_VISUALIZATION_AVAILABLE = True
except ImportError:
    INFORMATION_VISUALIZATION_AVAILABLE = False
    # Define placeholder functions
    def _not_available(*args, **kwargs):
        raise ImportError("Information theory visualization functions not available. Install information module dependencies.")
    
    plot_entropy_distribution = _not_available
    plot_entropy_landscape = _not_available
    plot_information_network = _not_available
    plot_information_profile = _not_available
    plot_mi_network = _not_available
    plot_mutual_information_matrix = _not_available
    plot_renyi_spectrum = _not_available
    plot_semantic_similarity_network = _not_available

__all__ = [
    "plot_entropy_distribution",
    "plot_entropy_landscape",
    "plot_information_network",
    "plot_information_profile",
    "plot_mi_network",
    "plot_mutual_information_matrix",
    "plot_renyi_spectrum",
    "plot_semantic_similarity_network",
    "INFORMATION_VISUALIZATION_AVAILABLE",
]

