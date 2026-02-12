"""Phylogenetic ecology subpackage.

Provides phylogenetic diversity metrics (Faith's PD, UniFrac),
phylogenetic community structure (NRI/NTI), phylogenetic signal
tests, and simple tree construction."""
from __future__ import annotations

from . import diversity

__all__ = ['diversity']
