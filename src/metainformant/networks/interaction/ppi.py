"""Compatibility facade for protein-protein interaction analysis.

The implementation lives in :mod:`metainformant.networks.interaction._ppi_impl`
so existing imports from ``metainformant.networks.interaction.ppi`` remain
stable while interaction analysis is split into clearer implementation modules.
"""

from __future__ import annotations

from metainformant.networks.interaction import _ppi_impl as _impl
from metainformant.networks.interaction._ppi_impl import *  # noqa: F403

__all__ = [name for name in vars(_impl) if not name.startswith("_")]
