"""Chromatin state learning submodule for epigenome analysis.

This submodule provides tools for learning and analyzing chromatin states
from histone modification data, including ChromHMM-style state discovery
using Gaussian mixture models, state assignment, biological interpretation,
enrichment analysis, and cross-condition comparison."""
from __future__ import annotations

from . import state_learning

__all__ = ['state_learning']
