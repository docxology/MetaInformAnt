"""Metagenomics analysis module for METAINFORMANT.

This module provides comprehensive tools for metagenomic analysis, including
amplicon-based community profiling (16S/ITS), shotgun metagenomics (assembly,
binning, profiling), functional annotation (gene prediction, pathway reconstruction),
specialized visualization for microbial ecology data, community diversity metrics,
and comparative/differential abundance analysis.
"""

from __future__ import annotations

from . import amplicon, comparative, diversity, functional, shotgun, visualization

__all__ = [
    "amplicon",
    "comparative",
    "diversity",
    "functional",
    "shotgun",
    "visualization",
]
