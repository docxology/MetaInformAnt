"""Structural variant filtering subpackage.

Provides quality-based filtering, size filtering, population frequency filtering,
blacklist region filtering, multi-caller consensus merging, and deduplication
of structural variant calls.
"""

from __future__ import annotations

from . import merge, quality_filter

__all__ = [
    "quality_filter",
    "merge",
]
