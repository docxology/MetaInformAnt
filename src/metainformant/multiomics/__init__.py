"""Multi-omics integration module for METAINFORMANT.

This module provides comprehensive tools for integrating multiple omics data layers,
enabling cross-platform analysis of genomics, transcriptomics, epigenomics, and
proteomics data.
"""

from __future__ import annotations

from . import analysis, methods, pathways, survival, visualization

__all__ = ["analysis", "methods", "pathways", "survival", "visualization"]
