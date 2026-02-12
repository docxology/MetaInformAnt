"""Ontology querying, traversal, and serialization.

Re-exports all public symbols from query and serialize modules for backward compatibility."""
from __future__ import annotations

from . import query, serialize

__all__ = ['query', 'serialize']
