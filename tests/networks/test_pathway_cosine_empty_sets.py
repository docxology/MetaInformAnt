"""Regression tests: pathway similarity with an empty gene set.

``pathway_similarity_analysis(method="cosine")`` divided the overlap by ``sqrt(|A|) * sqrt(|B|)``, so one pathway
without genes (an empty list, a dict with an empty ``genes`` entry, or an entry that is neither a list nor a dict with
``genes``) raised ``ZeroDivisionError`` for the whole analysis, and for ``create_pathway_network(method="cosine")``
too. Jaccard and overlap already returned 0.0 for that case. The cosine of an empty gene set is undefined; it is now
defined as 0.0, like the other two measures, and the docstring says so.

The oracle is the cosine of the 0/1 indicator vectors computed with NumPy over the union of all genes. No test
doubles.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from metainformant.networks.analysis.pathway import create_pathway_network, pathway_similarity_analysis

PATHWAYS = {
    "A": ["g1", "g2", "g3"],
    "B": ["g2", "g3", "g4", "g5"],
    "C": ["g5", "g6"],
    "EMPTY_LIST": [],
    "EMPTY_DICT": {"genes": []},
    "NO_GENES": {"description": "an entry without a gene list"},
}
EMPTY = {"EMPTY_LIST", "EMPTY_DICT", "NO_GENES"}


def _genes(entry) -> set[str]:
    if isinstance(entry, dict):
        return set(entry.get("genes", []))
    return set(entry)


def _cosine_oracle(first: str, second: str) -> float:
    vocabulary = sorted(set().union(*(_genes(entry) for entry in PATHWAYS.values())))
    a = np.array([g in _genes(PATHWAYS[first]) for g in vocabulary], dtype=float)
    b = np.array([g in _genes(PATHWAYS[second]) for g in vocabulary], dtype=float)
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / norm) if norm > 0 else 0.0


def test_cosine_with_an_empty_pathway_is_zero_instead_of_a_zero_division() -> None:
    similarities = pathway_similarity_analysis(PATHWAYS, method="cosine")["similarities"]

    for first, second in itertools.combinations(PATHWAYS, 2):
        value = similarities[f"{first}_{second}"]
        assert value == pytest.approx(_cosine_oracle(first, second)), (first, second)
        if first in EMPTY or second in EMPTY:
            assert value == 0.0


def test_cosine_of_nonempty_pathways_is_unchanged() -> None:
    similarities = pathway_similarity_analysis(PATHWAYS, method="cosine")["similarities"]

    assert similarities["A_B"] == pytest.approx(2 / np.sqrt(3 * 4))
    assert similarities["B_C"] == pytest.approx(1 / np.sqrt(4 * 2))
    assert similarities["A_C"] == 0.0


@pytest.mark.parametrize("method", ["jaccard", "overlap", "cosine"])
def test_every_method_gives_zero_for_pairs_with_an_empty_pathway(method: str) -> None:
    similarities = pathway_similarity_analysis(PATHWAYS, method=method)["similarities"]

    for first, second in itertools.combinations(PATHWAYS, 2):
        if first in EMPTY or second in EMPTY:
            assert similarities[f"{first}_{second}"] == 0.0, (first, second)


def test_two_empty_pathways_have_zero_similarity() -> None:
    assert pathway_similarity_analysis({"X": [], "Y": []}, method="cosine")["similarities"] == {"X_Y": 0.0}


def test_pathway_network_with_cosine_similarity_skips_the_empty_pathways() -> None:
    network = create_pathway_network(PATHWAYS, similarity_threshold=0.1, method="cosine")

    assert set(network.nodes()) == set(PATHWAYS)
    assert {frozenset(edge) for edge in network.edges()} == {frozenset(("A", "B")), frozenset(("B", "C"))}
    assert network["A"]["B"]["similarity"] == pytest.approx(_cosine_oracle("A", "B"))


def test_documentation_states_the_empty_pathway_convention() -> None:
    doc = " ".join((pathway_similarity_analysis.__doc__ or "").split())

    assert "no genes has similarity 0.0 with every other pathway" in doc
