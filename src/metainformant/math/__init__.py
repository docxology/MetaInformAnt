"""Mathematical biology and theoretical modeling module for METAINFORMANT.

This module provides mathematical and computational tools for biological modeling,
including population genetics theory, evolutionary dynamics, epidemiology,
coalescent theory, and quantitative genetics.
"""

from __future__ import annotations

# Import all mathematical biology submodules
from . import (
    coalescent,
    ddm,
    demography,
    dynamics,
    effective_size,
    egt,
    epidemiology,
    fst,
    ld,
    popgen,
    popgen_stats,
    price,
    quantgen,
    selection,
)

# Import selection_experiments submodule
from . import selection_experiments

# Optional imports with graceful fallbacks
try:
    from . import selection_experiments
except ImportError:
    selection_experiments = None

try:
    from . import ddm
except ImportError:
    ddm = None

try:
    from . import egt
except ImportError:
    egt = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Population genetics
    "popgen",
    "popgen_stats",
    "fst",
    "ld",
    "effective_size",

    # Evolutionary theory
    "coalescent",
    "selection",
    "price",

    # Population dynamics and ecology
    "demography",
    "dynamics",
    "epidemiology",

    # Decision theory and behavior
    "ddm",
    "egt",

    # Quantitative genetics
    "quantgen",

    # Specialized experiments
    "selection_experiments",
]
