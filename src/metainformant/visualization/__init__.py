"""Visualization and plotting utilities module for METAINFORMANT."""

from __future__ import annotations

from . import analysis, config, dashboards, genomics, interactive_dashboards, plots
from .config.conventions import OKABE_ITO  # noqa: F401  (re-export; single canonical source)

# Wong colorblind-safe palette (commonly used default)
WONG = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# Canonical Okabe-Ito palette: re-exported from the single definition in
# visualization.config.conventions (shared chart-convention contract).

__all__ = ["analysis", "config", "dashboards", "genomics", "interactive_dashboards", "plots", "WONG", "OKABE_ITO"]
