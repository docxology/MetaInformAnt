"""Configuration subpackage for visualization themes and palettes.

Re-exports all public APIs from :mod:`themes` and :mod:`palettes` for
convenience."""
from __future__ import annotations

import sys

from . import palettes, themes

# Self-reference so `from metainformant.visualization.config import config` works
config = sys.modules[__name__]

__all__ = ['palettes', 'themes', 'config']
