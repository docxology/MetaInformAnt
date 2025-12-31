"""Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

This module provides comprehensive GWAS analysis capabilities, including
quality control, population structure analysis, association testing,
multiple testing correction, and publication-quality visualization.
"""

from __future__ import annotations

# Import all GWAS analysis submodules
from . import (
    association,
    calling,
    config,
    correction,
    download,
    quality,
    sra_download,
    structure,
    visualization,
    visualization_comparison,
    visualization_effects,
    visualization_genome,
    visualization_population,
    visualization_regional,
    visualization_statistical,
    visualization_suite,
    visualization_variants,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import visualization_comparison
except ImportError:
    visualization_comparison = None

try:
    from . import visualization_effects
except ImportError:
    visualization_effects = None

try:
    from . import visualization_genome
except ImportError:
    visualization_genome = None

try:
    from . import visualization_population
except ImportError:
    visualization_population = None

try:
    from . import visualization_regional
except ImportError:
    visualization_regional = None

try:
    from . import visualization_statistical
except ImportError:
    visualization_statistical = None

try:
    from . import visualization_suite
except ImportError:
    visualization_suite = None

try:
    from . import visualization_variants
except ImportError:
    visualization_variants = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core GWAS analysis
    "quality",
    "structure",
    "association",
    "correction",
    "workflow",

    # Data acquisition
    "download",
    "sra_download",
    "calling",

    # Visualization suite
    "visualization",
    "visualization_comparison",
    "visualization_effects",
    "visualization_genome",
    "visualization_population",
    "visualization_regional",
    "visualization_statistical",
    "visualization_suite",
    "visualization_variants",

    # Configuration
    "config",
]


