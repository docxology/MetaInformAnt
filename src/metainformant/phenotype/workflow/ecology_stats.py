"""Phenotype-facing adapter for general ecological statistics.

PERMANOVA and PCoA are implemented in the ecology domain and are also useful
for multivariate phenotype matrices.  This adapter is the deliberate
integration seam between those domains.
"""

from __future__ import annotations

from metainformant.ecology.analysis.indicators import permanova
from metainformant.ecology.analysis.ordination import pcoa

__all__ = ["permanova", "pcoa"]
