"""Long-read visualization module for read length, quality, and alignment plots.

Generates publication-quality plots for read length distributions, quality
versus length scatter plots, sequence dotplots, alignment views, methylation
tracks, and phasing block diagrams.
"""

from __future__ import annotations

from . import plots

__all__ = [
    "plots",
]
