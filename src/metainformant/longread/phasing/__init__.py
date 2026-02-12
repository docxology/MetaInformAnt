"""Haplotype phasing submodule for long-read sequencing.

Provides read-based haplotype phasing, phase block construction,
switch error computation, read haplotagging, and allele-specific analysis."""
from __future__ import annotations

from . import haplotyping

__all__ = ['haplotyping']
