"""Long-read methylation calling and analysis submodule.

Provides signal-level methylation calling, per-site aggregation,
differentially methylated region detection, single-read pattern analysis,
and summary statistics for nanopore-derived methylation data."""
from __future__ import annotations

from . import calling

__all__ = ['calling']
