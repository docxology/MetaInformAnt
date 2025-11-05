"""Ontology and functional annotation utilities."""

from .go import enrich_genes, load_go_obo, semantic_similarity, validate_go_ontology, write_go_summary
from .query import (
    ancestors,
    clear_cache,
    common_ancestors,
    descendants,
    distance,
    filter_by_namespace,
    find_term_by_name,
    get_leaves,
    get_roots,
    path_to_root,
    set_cache_enabled,
    set_cache_ttl,
    subgraph,
)
from .serialize import load_ontology, save_ontology
from .types import Ontology, Term

__all__ = [
    "load_go_obo",
    "write_go_summary",
    "validate_go_ontology",
    "enrich_genes",
    "semantic_similarity",
    "ancestors",
    "descendants",
    "subgraph",
    "common_ancestors",
    "path_to_root",
    "distance",
    "find_term_by_name",
    "filter_by_namespace",
    "get_roots",
    "get_leaves",
    "save_ontology",
    "load_ontology",
    "clear_cache",
    "set_cache_enabled",
    "set_cache_ttl",
    "Ontology",
    "Term",
]
