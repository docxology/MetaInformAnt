"""Structural variant filtering subpackage.

Provides quality-based filtering, size filtering, population frequency filtering,
blacklist region filtering, multi-caller consensus merging, and deduplication
of structural variant calls.
"""

from __future__ import annotations

from . import quality_filter
from . import merge

__all__ = [
    "quality_filter",
    "merge",
]
