"""Shotgun metagenomics analysis.

Provides metagenome assembly, metagenomic binning, and community
profiling from whole-genome shotgun sequencing data.
"""

from __future__ import annotations

from . import assembly, binning, profiling

__all__ = [
    "assembly",
    "binning",
    "profiling",
]
