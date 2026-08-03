"""Public facade for Amalgkit CLI integration.

The implementation lives in :mod:`metainformant.rna.amalgkit._amalgkit_impl`
while this module exposes the stable public CLI-wrapper namespace.
"""

from __future__ import annotations

from metainformant.rna.amalgkit import _amalgkit_impl as _impl
from metainformant.rna.amalgkit._amalgkit_impl import *  # noqa: F401,F403

__all__ = [name for name in vars(_impl) if not name.startswith("_")]
