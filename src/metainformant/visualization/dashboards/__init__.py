"""Dashboards subpackage for composite figures and interactive plots.

Re-exports all public APIs from :mod:`composite` and :mod:`interactive`
for convenience."""
from __future__ import annotations

from . import composite, interactive

__all__ = ['composite', 'interactive']
