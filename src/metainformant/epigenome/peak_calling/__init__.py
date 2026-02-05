"""Peak calling submodule for epigenome analysis.

This submodule provides tools for identifying enriched regions (peaks)
in ChIP-seq and ATAC-seq signal data, including narrow peak calling,
broad domain detection, peak merging, filtering, and differential
peak analysis between conditions.
"""

from __future__ import annotations

from .peak_detection import (
    call_peaks_broad,
    call_peaks_simple,
    compute_frip,
    compute_local_lambda,
    differential_peaks,
    filter_peaks,
    merge_peaks,
    peak_summit_refinement,
)

__all__ = [
    "call_peaks_simple",
    "call_peaks_broad",
    "compute_local_lambda",
    "peak_summit_refinement",
    "merge_peaks",
    "filter_peaks",
    "compute_frip",
    "differential_peaks",
]
