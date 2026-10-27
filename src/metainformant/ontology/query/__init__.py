"""Ontology querying, traversal, and serialization.

Re-exports every public function of the :mod:`.query` and :mod:`.serialize` modules, so
``from metainformant.ontology.query import ancestors`` works as it did when ``query`` was a single module. The
submodules themselves stay importable (``from metainformant.ontology.query import query, serialize``)."""

from __future__ import annotations

from . import query, serialize
from .query import (
    ancestors,
    calculate_ic_map,
    clear_cache,
    common_ancestors,
    descendants,
    distance,
    filter_by_namespace,
    find_term_by_name,
    find_terms_by_name,
    find_terms_by_namespace,
    get_leaves,
    get_roots,
    get_subontology,
    get_subontology_stats,
    information_content,
    most_informative_common_ancestor,
    path_to_root,
    set_cache_enabled,
    shortest_path,
    subgraph,
    validate_ontology_integrity,
)
from .serialize import (
    export_ontology_stats,
    graph_to_ontology,
    load_ontology,
    merge_ontologies,
    ontology_to_graph,
    save_ontology,
)

__all__ = [
    "query",
    "serialize",
    # query
    "ancestors",
    "calculate_ic_map",
    "clear_cache",
    "common_ancestors",
    "descendants",
    "distance",
    "filter_by_namespace",
    "find_term_by_name",
    "find_terms_by_name",
    "find_terms_by_namespace",
    "get_leaves",
    "get_roots",
    "get_subontology",
    "get_subontology_stats",
    "information_content",
    "most_informative_common_ancestor",
    "path_to_root",
    "set_cache_enabled",
    "shortest_path",
    "subgraph",
    "validate_ontology_integrity",
    # serialize
    "export_ontology_stats",
    "graph_to_ontology",
    "load_ontology",
    "merge_ontologies",
    "ontology_to_graph",
    "save_ontology",
]
