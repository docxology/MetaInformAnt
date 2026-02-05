"""Multi-omic survival analysis.

Provides Cox regression, Kaplan-Meier estimation, log-rank tests,
and multi-omic survival models with regularized feature selection.

Submodules:
    - analysis: Cox PH, Kaplan-Meier, log-rank, multi-omic survival models
"""

from __future__ import annotations

from . import analysis

__all__ = [
    "analysis",
]
