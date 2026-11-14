"""Ontology querying, traversal, and serialization.

The functions live in the :mod:`.query` (traversal, search, statistics) and :mod:`.serialize` (OBO, JSON and
NetworkX I/O) modules, for example ``from metainformant.ontology.query.query import ancestors``."""

from __future__ import annotations

from . import query, serialize

__all__ = ["query", "serialize"]
