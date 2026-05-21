"""Compatibility facade for general GWAS visualization utilities.

The implementation lives in :mod:`metainformant.gwas.visualization._general_impl`
so this public path can remain stable while the visualization package is split
into smaller modules.
"""

from __future__ import annotations

from metainformant.gwas.visualization import _general_impl as _impl
from metainformant.gwas.visualization._general_impl import *  # noqa: F403

__all__ = [name for name in vars(_impl) if not name.startswith("_")]
