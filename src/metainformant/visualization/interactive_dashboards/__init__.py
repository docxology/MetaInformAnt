"""Interactive dashboard subpackage for visualization.

Provides interactive scatter plots, heatmaps, genome browser tracks,
volcano plots, HTML export, and dashboard composition -- all as
JSON-serializable data structures suitable for web rendering.
"""

from __future__ import annotations

from .dashboards import (
    create_dashboard,
    create_genome_browser_track,
    create_interactive_heatmap,
    create_interactive_scatter,
    create_interactive_volcano,
    export_to_html,
)

__all__ = [
    "create_interactive_scatter",
    "create_interactive_heatmap",
    "create_genome_browser_track",
    "create_interactive_volcano",
    "export_to_html",
    "create_dashboard",
]
