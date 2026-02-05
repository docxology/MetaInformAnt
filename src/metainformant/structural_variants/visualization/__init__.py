"""Structural variant visualization subpackage.

Provides publication-quality plots for structural variant analysis including
Circos-style genome-wide views, coverage tracks with SV overlays,
size distribution histograms, SV type summaries, breakpoint detail views,
and genome-wide CNV profiles.
"""

from __future__ import annotations

from . import plots

__all__ = [
    "plots",
]
