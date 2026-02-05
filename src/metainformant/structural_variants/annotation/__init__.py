"""Structural variant annotation subpackage.

Provides tools for annotating structural variants with gene overlaps,
regulatory element overlaps, and functional impact predictions including
dosage sensitivity assessment and TAD boundary disruption analysis.
"""

from __future__ import annotations

from . import overlap
from . import functional_impact

__all__ = [
    "overlap",
    "functional_impact",
]
