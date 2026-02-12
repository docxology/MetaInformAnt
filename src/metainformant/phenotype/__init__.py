"""Phenotype module for MetaInformAnt.

Comprehensive phenotypic trait analysis including:
- Morphological measurements and allometric analysis
- Behavioral sequences, ethograms, and diversity metrics
- Chemical profiles (GC-MS, CHC) with distance and marker detection
- Acoustic signals with spectral and temporal analysis
- Electronic tracking (RFID, video, GPS) with movement ecology
- Life course trajectory analysis
- Cross-omic integration (phenotype-genotype, trait-expression)
- Configurable analysis pipelines
"""

from __future__ import annotations

from . import (
    analysis,
    behavior,
    chemical,
    data,
    electronic,
    gwas_integration,
    integration,
    morphological,
    sonic,
    visualization,
    workflow,
)

__all__ = [
    "analysis",
    "behavior",
    "chemical",
    "data",
    "electronic",
    "gwas_integration",
    "integration",
    "morphological",
    "sonic",
    "visualization",
    "workflow",
]
