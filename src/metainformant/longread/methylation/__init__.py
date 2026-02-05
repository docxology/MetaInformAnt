"""Long-read methylation calling and analysis submodule.

Provides signal-level methylation calling, per-site aggregation,
differentially methylated region detection, single-read pattern analysis,
and summary statistics for nanopore-derived methylation data.
"""

from __future__ import annotations

from . import calling
from .calling import (
    aggregate_methylation,
    call_methylation_from_signal,
    compute_methylation_stats,
    detect_dmrs,
    methylation_pattern_analysis,
)

__all__ = [
    "calling",
    "aggregate_methylation",
    "call_methylation_from_signal",
    "compute_methylation_stats",
    "detect_dmrs",
    "methylation_pattern_analysis",
]
