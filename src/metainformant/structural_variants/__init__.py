"""Structural Variant (SV) analysis module for METAINFORMANT.

This module provides comprehensive tools for detecting, annotating, filtering,
and visualizing structural variants including copy number variations (CNVs),
deletions, duplications, inversions, translocations, and insertions.

Config prefix: SV_
"""

from __future__ import annotations

from . import annotation, detection, filtering, population, visualization

__all__ = [
    "annotation",
    "detection",
    "filtering",
    "population",
    "visualization",
]
