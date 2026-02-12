"""Spatial data I/O module for METAINFORMANT.

Provides loaders for spatial transcriptomics platforms including
10x Visium, MERFISH, and 10x Xenium."""
from __future__ import annotations

from . import merfish, visium, xenium

__all__ = ['merfish', 'visium', 'xenium']
