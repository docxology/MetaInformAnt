"""Regression tests: ``get_subontology`` keeps only the hierarchy that lies inside the sub-ontology.

``get_subontology`` built its result with ``create_ontology``, i.e. through ``Ontology.add_term``, which takes every
parent from ``Term.is_a_parents`` whether or not that parent is part of the sub-ontology. For a seed that is not a
root of the source ontology (every seed but the root of a namespace) the result therefore contained links to terms it
does not hold: its ``parents_of`` marked the seed as having a parent, so ``get_roots`` returned nothing,
``get_ancestors`` wandered out of the sub-ontology, and ``relationships`` held dangling and duplicated ``is_a``
entries. It now builds the maps with the same helper as ``filter_by_namespace`` and ``subgraph`` (links with both ends
inside, nothing else). A ``root_terms`` generator is also read once, from the materialised list, instead of being
exhausted by the validation.

The oracle is an independent ``networkx`` digraph (edges child -> parent): the sub-ontology is the seeds plus
everything that reaches a seed, roots are nodes without outgoing edges, ancestors are ``nx.descendants`` in that
digraph. No test doubles.
"""

from __future__ import annotations

import random
from collections import Counter

import networkx as nx
import pytest

from metainformant.ontology.core.types import Ontology, Relationship, Term
from metainformant.ontology.query.query import ancestors, clear_cache, descendants, get_subontology

SPEC: dict[str, list[str]] = {
    "BP:root": [],
    "BP:a": ["BP:root"],
    "BP:b": ["BP:a"],
    "MF:root": [],
    "MF:a": ["MF:root"],
    "BP:cross": ["BP:b", "MF:a"],  # a parent outside any sub-ontology rooted in the BP branch
    "BP:leaf": ["BP:cross"],
}


@pytest.fixture(autouse=True)
def fresh_cache():
    clear_cache()
    yield
    clear_cache()


def _ontology(spec: dict[str, list[str]]) -> Ontology:
    onto = Ontology()
    for term_id, parents in spec.items():
        onto.add_term(Term(term_id, f"name of {term_id}", is_a_parents=list(parents)))
    return onto


def _oracle(spec: dict[str, list[str]], seeds: list[str]) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(spec)
    graph.add_edges_from((term, parent) for term, parents in spec.items() for parent in parents)
    keep = set(seeds)
    for seed in seeds:
        keep |= nx.ancestors(graph, seed)  # everything that is_a-reaches the seed, i.e. its descendants in the ontology
    return graph.subgraph(keep).copy()


def _assert_matches(sub: Ontology, oracle: nx.DiGraph) -> None:
    assert set(sub.terms) == set(oracle.nodes)
    assert sub.get_roots() == {node for node in oracle if oracle.out_degree(node) == 0}
    assert sub.get_leaves() == {node for node in oracle if oracle.in_degree(node) == 0}
    for node in oracle:
        assert set(sub.get_parents(node)) == set(oracle.successors(node)), node
        assert set(sub.get_children(node)) == set(oracle.predecessors(node)), node
        assert ancestors(sub, node) == nx.descendants(oracle, node), node
        assert descendants(sub, node) == nx.ancestors(oracle, node), node
    # every relationship joins two terms of the sub-ontology, and each is_a link is recorded exactly once
    assert all(rel.source in sub.terms and rel.target in sub.terms for rel in sub.relationships)
    assert Counter((r.source, r.target, r.relation_type) for r in sub.relationships) == Counter(
        (source, target, "is_a") for source, target in oracle.edges
    )


def test_a_non_root_seed_is_a_root_of_its_sub_ontology() -> None:
    sub = get_subontology(_ontology(SPEC), ["BP:a"])

    assert set(sub.terms) == {"BP:a", "BP:b", "BP:cross", "BP:leaf"}
    assert sub.get_roots() == {
        "BP:a"
    }  # it used to be empty: BP:a still pointed at BP:root, which is not in the sub-ontology
    assert ancestors(sub, "BP:leaf") == {"BP:cross", "BP:b", "BP:a"}  # not BP:root, not MF:a, not MF:root


@pytest.mark.parametrize(
    "seeds", [["BP:root"], ["BP:a"], ["BP:b"], ["BP:cross"], ["BP:leaf"], ["MF:a", "BP:b"], ["BP:leaf", "MF:root"]]
)
def test_sub_ontology_matches_the_networkx_oracle(seeds: list[str]) -> None:
    _assert_matches(get_subontology(_ontology(SPEC), seeds), _oracle(SPEC, seeds))


def test_a_generator_of_seeds_is_not_exhausted_by_the_validation() -> None:
    sub = get_subontology(_ontology(SPEC), (seed for seed in ["BP:b", "MF:a"]))

    assert set(sub.terms) == set(_oracle(SPEC, ["BP:b", "MF:a"]).nodes)


def test_random_ontologies_match_the_networkx_oracle() -> None:
    rng = random.Random(20261027)
    for _ in range(25):
        spec: dict[str, list[str]] = {}
        for index in range(18):
            earlier = list(spec)
            spec[f"T{index:02d}"] = rng.sample(earlier, rng.randint(0, min(3, len(earlier)))) if earlier else []
        seeds = rng.sample(list(spec), rng.randint(1, 3))

        clear_cache()
        _assert_matches(get_subontology(_ontology(spec), seeds), _oracle(spec, seeds))


def test_other_relation_types_carry_only_that_relation() -> None:
    onto = _ontology({"A": [], "B": ["A"], "C": [], "D": []})
    onto.add_relationship(Relationship("C", "A", "part_of"))
    onto.add_relationship(Relationship("D", "C", "part_of"))

    sub = get_subontology(onto, ["A", "C"], relation_type="part_of")

    assert set(sub.terms) == {"A", "C", "D"}
    assert Counter((r.source, r.target, r.relation_type) for r in sub.relationships) == Counter(
        {("C", "A", "part_of"): 1, ("D", "C", "part_of"): 1}
    )
    assert sub.parents_of == {} and sub.children_of == {}  # is_a was not requested
    assert sub.get_ancestors("D", "part_of") == {"C", "A"}


def test_the_source_ontology_and_its_terms_are_unchanged() -> None:
    onto = _ontology(SPEC)
    relationships_before = Counter((r.source, r.target, r.relation_type) for r in onto.relationships)
    parents_before = {term: set(parents) for term, parents in onto.parents_of.items()}

    get_subontology(onto, ["BP:b"])

    assert Counter((r.source, r.target, r.relation_type) for r in onto.relationships) == relationships_before
    assert onto.parents_of == parents_before
    assert onto.get_roots() == {"BP:root", "MF:root"}


def test_metadata_is_copied() -> None:
    onto = _ontology({"A": [], "B": ["A"]})
    onto.metadata.update({"id": "demo", "terms": "a key that is also a keyword argument of create_ontology"})

    sub = get_subontology(onto, ["A"])

    assert sub.metadata == onto.metadata and sub.metadata is not onto.metadata
