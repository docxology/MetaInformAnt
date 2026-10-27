"""Regression tests for ``filter_by_namespace``, ``subgraph`` and the ancestors/descendants query cache.

* ``filter_by_namespace`` and ``subgraph`` built the new ``Ontology`` from ``terms`` and ``relationships`` only.
  ``Ontology`` answers ``get_roots``, ``get_parents``, ``get_children`` and ``get_ancestors`` for ``is_a`` from its
  ``parents_of`` / ``children_of`` maps, so on the result every term was a root and every ancestor set was empty. The
  maps are now carried over for the terms that are kept.
* The query cache behind ``ancestors`` / ``descendants`` was keyed by ``id(onto)`` and never invalidated: results went
  stale as soon as terms or relationships were added to a cached ontology, and a new ontology that reused the address
  of a garbage-collected one was served the old one's answers. Entries are now tied to the ontology object (weak
  reference) and to its structure, and they disappear with it.

The oracle for the sub-ontologies is an independent ``networkx`` digraph (edges child -> parent): roots are nodes
without outgoing edges, ancestors are ``nx.descendants`` in that digraph. No test doubles.
"""

from __future__ import annotations

import gc
import random

import networkx as nx
import pytest

from metainformant.ontology.core.types import Ontology, Relationship, Term
from metainformant.ontology.query import query as query_module
from metainformant.ontology.query.query import (
    ancestors,
    clear_cache,
    descendants,
    filter_by_namespace,
    set_cache_enabled,
    subgraph,
)


@pytest.fixture(autouse=True)
def fresh_cache():
    clear_cache()
    set_cache_enabled(True)
    yield
    set_cache_enabled(True)
    clear_cache()


def _ontology(spec: dict[str, tuple[str, list[str]]]) -> Ontology:
    """Build an ontology from ``{term_id: (namespace, [is_a parents])}`` through the real ``add_term``."""
    onto = Ontology()
    for term_id, (namespace, parents) in spec.items():
        onto.add_term(Term(term_id, f"name of {term_id}", namespace=namespace, is_a_parents=list(parents)))
    return onto


GO_LIKE = {
    "BP:root": ("biological_process", []),
    "BP:a": ("biological_process", ["BP:root"]),
    "BP:b": ("biological_process", ["BP:a"]),
    "MF:root": ("molecular_function", []),
    "MF:a": ("molecular_function", ["MF:root"]),
    "BP:cross": ("biological_process", ["BP:b", "MF:a"]),  # one parent in another namespace
    "BP:orphan": ("biological_process", ["MF:a"]),  # its only parent is in another namespace
}


def _induced_digraph(spec: dict[str, tuple[str, list[str]]], keep: set[str]) -> nx.DiGraph:
    """Independent oracle: the child -> parent digraph restricted to ``keep``."""
    graph = nx.DiGraph()
    graph.add_nodes_from(spec)
    graph.add_edges_from((term, parent) for term, (_, parents) in spec.items() for parent in parents)
    return graph.subgraph(keep).copy()


def _assert_matches_oracle(sub: Ontology, oracle: nx.DiGraph) -> None:
    assert set(sub.terms) == set(oracle.nodes)
    assert sub.get_roots() == {node for node in oracle if oracle.out_degree(node) == 0}
    assert sub.get_leaves() == {node for node in oracle if oracle.in_degree(node) == 0}
    for node in oracle:
        assert set(sub.get_parents(node)) == set(oracle.successors(node)), node
        assert set(sub.get_children(node)) == set(oracle.predecessors(node)), node
        assert ancestors(sub, node) == nx.descendants(oracle, node), node
        assert descendants(sub, node) == nx.ancestors(oracle, node), node


# ---------------------------------------------------------------------------
# filter_by_namespace / subgraph keep the hierarchy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("namespace", ["biological_process", "molecular_function"])
def test_filter_by_namespace_keeps_parents_children_roots_and_ancestors(namespace: str) -> None:
    onto = _ontology(GO_LIKE)
    keep = {term for term, (ns, _) in GO_LIKE.items() if ns == namespace}

    _assert_matches_oracle(filter_by_namespace(onto, namespace), _induced_digraph(GO_LIKE, keep))


def test_filter_by_namespace_makes_terms_with_only_foreign_parents_roots() -> None:
    sub = filter_by_namespace(_ontology(GO_LIKE), "biological_process")

    assert sub.get_roots() == {"BP:root", "BP:orphan"}
    assert ancestors(sub, "BP:cross") == {"BP:b", "BP:a", "BP:root"}  # the molecular_function parent is gone
    assert ancestors(sub, "BP:orphan") == set()


def test_filter_by_namespace_does_not_change_the_source_ontology() -> None:
    onto = _ontology(GO_LIKE)
    roots_before = onto.get_roots()
    parents_before = {term: set(parents) for term, parents in onto.parents_of.items()}

    filter_by_namespace(onto, "biological_process")

    assert onto.get_roots() == roots_before == {"BP:root", "MF:root"}
    assert onto.parents_of == parents_before


def test_subgraph_keeps_the_links_between_selected_terms() -> None:
    onto = _ontology(GO_LIKE)
    keep = {"BP:root", "BP:a", "BP:b", "BP:cross"}

    _assert_matches_oracle(subgraph(onto, sorted(keep)), _induced_digraph(GO_LIKE, keep))


def test_subgraph_does_not_bridge_over_a_dropped_intermediate_term() -> None:
    sub = subgraph(_ontology(GO_LIKE), ["BP:b", "BP:root"])  # BP:a, the term between them, is not selected

    assert sub.get_roots() == {"BP:b", "BP:root"}
    assert ancestors(sub, "BP:b") == set()


def test_subgraph_with_another_relation_type_carries_only_that_relation() -> None:
    onto = _ontology({"A": ("ns", []), "B": ("ns", ["A"]), "C": ("ns", [])})
    onto.add_relationship(Relationship("C", "A", "part_of"))

    sub = subgraph(onto, ["A", "B", "C"], relation_type="part_of")

    assert [(r.source, r.target, r.relation_type) for r in sub.relationships] == [("C", "A", "part_of")]
    assert sub.get_parents("C", "part_of") == {"A"}
    assert sub.parents_of == {} and sub.children_of == {}  # is_a was not requested


def test_random_ontologies_match_the_networkx_oracle() -> None:
    rng = random.Random(20261027)
    namespaces = ["biological_process", "molecular_function", "cellular_component"]
    for trial in range(25):
        spec: dict[str, tuple[str, list[str]]] = {}
        for index in range(18):
            earlier = list(spec)
            parents = rng.sample(earlier, rng.randint(0, min(3, len(earlier)))) if earlier else []
            spec[f"T{index:02d}"] = (rng.choice(namespaces), parents)
        onto = _ontology(spec)
        chosen = rng.sample(list(spec), 9)

        clear_cache()
        _assert_matches_oracle(subgraph(onto, chosen), _induced_digraph(spec, set(chosen)))
        namespace = rng.choice(namespaces)
        keep = {term for term, (ns, _) in spec.items() if ns == namespace}
        if keep:
            _assert_matches_oracle(filter_by_namespace(onto, namespace), _induced_digraph(spec, keep))


# ---------------------------------------------------------------------------
# query cache
# ---------------------------------------------------------------------------


def _chain() -> Ontology:
    return _ontology({"A": ("ns", []), "B": ("ns", ["A"]), "C": ("ns", ["B"])})


def test_cached_results_follow_terms_added_later() -> None:
    onto = _chain()
    assert descendants(onto, "A") == {"B", "C"}
    assert ancestors(onto, "C") == {"A", "B"}

    onto.add_term(Term("D", "grandgrandchild", is_a_parents=["C"]))

    assert descendants(onto, "A") == {"B", "C", "D"}
    assert ancestors(onto, "D") == {"A", "B", "C"}


def test_cached_results_follow_relationships_added_later() -> None:
    onto = _ontology({"A": ("ns", []), "B": ("ns", ["A"]), "Z": ("ns", [])})
    assert ancestors(onto, "B") == {"A"}
    assert descendants(onto, "Z") == set()

    onto.add_relationship(Relationship("A", "Z", "is_a"))  # A is_a Z

    assert ancestors(onto, "B") == {"A", "Z"}
    assert descendants(onto, "Z") == {"A", "B"}


def test_cached_results_follow_other_relation_types_too() -> None:
    onto = _ontology({"A": ("ns", []), "B": ("ns", []), "C": ("ns", [])})
    onto.add_relationship(Relationship("B", "A", "part_of"))
    assert ancestors(onto, "B", "part_of") == {"A"}

    onto.add_relationship(Relationship("A", "C", "part_of"))

    assert ancestors(onto, "B", "part_of") == {"A", "C"}


def test_an_ontology_that_reuses_a_collected_ontology_address_gets_its_own_answers() -> None:
    reused = False
    for _ in range(100):
        first = _chain()
        assert ancestors(first, "C") == {"A", "B"}
        assert descendants(first, "A") == {"B", "C"}
        first_address = id(first)
        del first  # CPython frees it at once; the next object of the same size usually takes over its address

        # same term ids, different structure: stale answers from the collected ontology would be wrong here
        second = _ontology({"A": ("ns", []), "B": ("ns", []), "C": ("ns", ["A"])})
        assert ancestors(second, "C") == {"A"}
        assert descendants(second, "A") == {"C"}
        assert descendants(second, "B") == set()
        if id(second) == first_address:
            reused = True
            break
        del second
    assert reused, "the address of the collected ontology was never reused; the scenario did not occur"


def test_cache_entries_do_not_outlive_their_ontology() -> None:
    baseline = len(query_module._query_cache)

    for _ in range(25):
        onto = _chain()
        ancestors(onto, "C")
        descendants(onto, "A")
        del onto

    gc.collect()
    assert len(query_module._query_cache) == baseline


def test_cache_is_used_for_repeated_queries_and_can_be_disabled_and_cleared() -> None:
    onto = _chain()

    first = ancestors(onto, "C")
    first.add("junk")  # a caller changing its copy must not corrupt the cache
    assert ancestors(onto, "C") == {"A", "B"}
    assert len(query_module._query_cache) == 1

    clear_cache()
    assert len(query_module._query_cache) == 0

    set_cache_enabled(False)
    assert ancestors(onto, "C") == {"A", "B"}
    assert len(query_module._query_cache) == 0
