"""``add_edges_from_dataframe(directed=True)`` cannot silently store undirected edges.

The ``directed`` argument ("Whether edges are directed") was never read, so ``directed=True`` with a plain
``nx.Graph`` stored every edge undirected and a caller who wanted the direction of a regulatory edge list lost it
without a message. The edges go into the graph the caller passes, so the only honest reading of the flag is as an
assertion about that graph: ``directed=True`` needs a ``DiGraph``. No test doubles.
"""

from __future__ import annotations

import networkx as nx
import pandas as pd
import pytest

from metainformant.networks.analysis.graph_core import add_edges_from_dataframe

EDGES = pd.DataFrame({"src": ["A", "B", "C"], "tgt": ["B", "C", "A"], "w": [1.0, 2.0, 3.0]})


def test_directed_edges_into_an_undirected_graph_are_an_error() -> None:
    graph = nx.Graph()

    with pytest.raises(ValueError, match="directed"):
        add_edges_from_dataframe(graph, EDGES, "src", "tgt", directed=True)

    assert graph.number_of_edges() == 0  # nothing was half-added


def test_directed_edges_keep_their_direction_in_a_digraph() -> None:
    graph = nx.DiGraph()

    add_edges_from_dataframe(graph, EDGES, "src", "tgt", ["w"], directed=True)

    assert set(graph.edges) == {("A", "B"), ("B", "C"), ("C", "A")}
    assert not graph.has_edge("B", "A")
    assert graph["B"]["C"]["w"] == 2.0


def test_the_default_still_fills_whatever_graph_is_given() -> None:
    undirected, directed = nx.Graph(), nx.DiGraph()

    add_edges_from_dataframe(undirected, EDGES, "src", "tgt")
    add_edges_from_dataframe(directed, EDGES, "src", "tgt")

    assert undirected.number_of_edges() == 3 and undirected.has_edge("B", "A")
    assert set(directed.edges) == {("A", "B"), ("B", "C"), ("C", "A")}
