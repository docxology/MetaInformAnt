"""Visualization adapter for RNA campaign observation bundles.

Figure rendering belongs to the visualization domain. RNA engine modules import
it through this integration module, the package-boundary surface for cross-domain
dependencies, rather than reaching into ``metainformant.visualization`` directly.
"""

from __future__ import annotations

from metainformant.visualization.plots.rna_campaign import plot_acquisition_observation

__all__ = ["plot_acquisition_observation"]
