"""Haplotype phasing submodule for long-read sequencing.

Provides read-based haplotype phasing, phase block construction,
switch error computation, read haplotagging, and allele-specific analysis.
"""

from __future__ import annotations

from . import haplotyping
from .haplotyping import (
    allele_specific_analysis,
    build_phase_blocks,
    compute_switch_errors,
    haplotag_reads,
    phase_reads,
)

__all__ = [
    "haplotyping",
    "allele_specific_analysis",
    "build_phase_blocks",
    "compute_switch_errors",
    "haplotag_reads",
    "phase_reads",
]
