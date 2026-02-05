"""Core multi-omic integration algorithms.

Provides matrix factorization, clustering, and network fusion methods
for integrating data across multiple omic layers.

Submodules:
    - factorization: Joint NMF, MOFA, tensor decomposition, SNF, CCA
    - clustering: Multi-omic clustering, consensus clustering, spectral methods
"""

from __future__ import annotations

from . import factorization
from . import clustering

__all__ = [
    "factorization",
    "clustering",
]
