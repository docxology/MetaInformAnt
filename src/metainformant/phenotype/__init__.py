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

import importlib


def _optional_import(name: str):
    try:
        return importlib.import_module(f"{__name__}.{name}")
    except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional scientific plotting stack
        if exc.name and not exc.name.startswith(f"{__name__}.{name}"):
            return None
        raise


analysis = _optional_import("analysis")
behavior = _optional_import("behavior")
chemical = _optional_import("chemical")
data = _optional_import("data")
electronic = _optional_import("electronic")
gwas_integration = _optional_import("gwas_integration")
integration = _optional_import("integration")
morphological = _optional_import("morphological")
sonic = _optional_import("sonic")
visualization = _optional_import("visualization")
workflow = _optional_import("workflow")

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
