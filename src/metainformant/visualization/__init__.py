"""Visualization and plotting utilities module for METAINFORMANT."""

from __future__ import annotations

from . import analysis, config, dashboards, genomics, interactive_dashboards, plots
from .config.conventions import OKABE_ITO  # noqa: F401  (re-export; single canonical source)

# Wong colorblind-safe palette (commonly used default): derived from the
# canonical Okabe-Ito definition in visualization.config.conventions
# (single-source contract; palettes.WONG orders black first).
WONG = list(config.palettes.WONG)

__all__ = ["analysis", "config", "dashboards", "genomics", "interactive_dashboards", "plots", "WONG", "OKABE_ITO"]
