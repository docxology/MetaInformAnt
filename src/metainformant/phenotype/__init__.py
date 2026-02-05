"""
Phenotype module for MetaInformAnt.

This module encapsulates phenotypic trait analysis, including:
- Morphological measurements (from AntWiki or direct measurement)
- Behavioral sequences and ethograms
- Chemical profiles (GC-MS)
- Acoustic signals (Stridulation)
- Electronic tracking data (RFID, video)

Submodules:
    antwiki: Legacy AntWiki data integration (loading/scraping).
    analysis: Life course and temporal analysis.
    behavior: Behavioral phenotype analysis.
    chemical: Chemical phenotype analysis.
    electronic: Electronic/Sensor phenotype analysis.
    morphological: Morphometric phenotype analysis.
    sonic: Acoustic phenotype analysis.
    visualization: Phenotype visualization tools.
    data: Data loading utilities.
"""

# Export new submodules
from . import behavior
from . import chemical
from . import electronic
from . import morphological
from . import sonic

# Export legacy/existing components
from . import analysis
from . import data
from . import visualization

__all__ = [
    "behavior",
    "chemical",
    "electronic",
    "morphological",
    "sonic",
    "analysis",
    "data",
    "visualization",
]
