"""Information theory metrics subpackage.

Provides syntactic, semantic, continuous, estimation, analysis, advanced,
hypothesis testing, channel capacity, and information geometry modules.

Subpackages:
    core/       - syntactic, continuous, estimation
    advanced/   - semantic, geometry, channel, hypothesis, decomposition
    analysis/   - analysis, advanced_analysis"""
from __future__ import annotations

from . import advanced, analysis, core

__all__ = ['advanced', 'analysis', 'core']
