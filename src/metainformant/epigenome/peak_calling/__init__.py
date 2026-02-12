"""Peak calling submodule for epigenome analysis.

This submodule provides tools for identifying enriched regions (peaks)
in ChIP-seq and ATAC-seq signal data, including narrow peak calling,
broad domain detection, peak merging, filtering, and differential
peak analysis between conditions."""
from __future__ import annotations

from . import peak_detection

__all__ = ['peak_detection']
