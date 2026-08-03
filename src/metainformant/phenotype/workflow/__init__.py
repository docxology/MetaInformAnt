"""Phenotype workflow adapters for cross-domain statistical routines."""

from __future__ import annotations

from .ecology_stats import pcoa, permanova

__all__ = ["permanova", "pcoa"]
