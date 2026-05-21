"""Compatibility facade for amalgkit CLI integration.

The implementation lives in :mod:`metainformant.rna.amalgkit._amalgkit_impl`
so existing imports from ``metainformant.rna.amalgkit.amalgkit`` remain valid
while RNA workflow internals continue to be decomposed.
"""

from __future__ import annotations

from metainformant.rna.amalgkit import _amalgkit_impl as _impl
from metainformant.rna.amalgkit._amalgkit_impl import *  # noqa: F403

__all__ = [name for name in vars(_impl) if not name.startswith("_")]
