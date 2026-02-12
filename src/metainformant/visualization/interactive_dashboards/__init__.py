"""Interactive dashboard subpackage for visualization.

Provides interactive scatter plots, heatmaps, genome browser tracks,
volcano plots, HTML export, and dashboard composition -- all as
JSON-serializable data structures suitable for web rendering."""
from __future__ import annotations

from . import dashboards

__all__ = ['dashboards']
