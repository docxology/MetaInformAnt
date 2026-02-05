"""Long-read assembly module for overlap computation, consensus, and hybrid assembly.

Provides minimizer-based overlap detection, partial-order alignment consensus
generation, iterative polishing, and hybrid assembly combining long and short reads.
"""

from __future__ import annotations

from . import overlap
from . import consensus
from . import hybrid

__all__ = [
    "overlap",
    "consensus",
    "hybrid",
]
