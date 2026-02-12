"""Configuration subpackage for visualization themes and palettes.

Re-exports all public APIs from :mod:`themes` and :mod:`palettes` for
convenience."""
from __future__ import annotations

from . import palettes, themes

__all__ = ['palettes', 'themes']
