"""Ecological and environmental analysis module for METAINFORMANT.

This module provides tools for ecological data analysis, community ecology,
environmental correlations, species interactions, and biodiversity assessment.
"""

from __future__ import annotations

# Import all ecology submodules
from . import (
    community,
    environmental,
    interactions,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import interactions
except ImportError:
    interactions = None

try:
    from . import workflow
except ImportError:
    workflow = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core ecological analysis
    "community",
    "environmental",

    # Advanced analysis
    "interactions",
    "workflow",
]



