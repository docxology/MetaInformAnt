"""Ontology and functional annotation utilities."""

from .go import load_go_obo, write_go_summary
from .query import ancestors, descendants, subgraph

__all__ = [
    "load_go_obo",
    "write_go_summary",
    "ancestors",
    "descendants",
    "subgraph",
]
