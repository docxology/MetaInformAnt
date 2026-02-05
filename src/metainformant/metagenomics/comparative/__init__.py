"""Comparative metagenomics analysis.

Provides differential abundance testing (ALDEx2-like, ANCOM-like),
indicator species analysis, effect size ranking (LEfSe-style), and
ML-based biomarker discovery for microbiome studies.
"""

from __future__ import annotations

from . import differential_abundance

__all__ = [
    "differential_abundance",
]
