"""Multi-omic pathway analysis.

Provides methods for combining pathway-level signals across multiple omic
layers, including enrichment analysis, active module detection, and
cross-omic concordance assessment.

Submodules:
    - enrichment: Multi-omic enrichment, active modules, topology analysis
"""

from __future__ import annotations

from . import enrichment

__all__ = [
    "enrichment",
]
