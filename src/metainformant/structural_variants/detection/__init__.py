"""Structural variant detection subpackage.

Provides algorithms for detecting structural variants from sequencing data,
including copy number variation detection via read depth analysis,
structural variant calling from split/discordant reads, and
breakpoint refinement to base-pair resolution.
"""

from __future__ import annotations

from . import cnv
from . import sv_calling
from . import breakpoints

__all__ = [
    "cnv",
    "sv_calling",
    "breakpoints",
]
