"""Visualization and plotting utilities module for METAINFORMANT."""

from __future__ import annotations

from . import analysis, config, dashboards, genomics, interactive_dashboards, plots

# Wong colorblind-safe palette (commonly used default)
WONG = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# Canonical Okabe-Ito palette (shared chart-convention contract, see
# visualization.config.conventions). Kept in sync by
# tests/visualization/test_conventions.py.
OKABE_ITO = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
)

__all__ = ["analysis", "config", "dashboards", "genomics", "interactive_dashboards", "plots", "WONG", "OKABE_ITO"]
