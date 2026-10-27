"""Regression tests for ``ppi_network_analysis`` centrality and the PPI JSON round trip.

* ``ppi_network_analysis`` labelled ``dict(list(centrality.items())[:10])`` as the "top 10" centrality. That is the
  first ten nodes in insertion order, whatever their centrality, so hubs were routinely missing from the report. It
  now returns the ten most central nodes, most central first.
* ``save_ppi_network(format="json")`` writes a node-link file, but ``load_ppi_network`` and
  ``ProteinNetwork.from_file`` rejected ``format="json"`` (the ``from_file`` docstring listed it). The loader now
  reads what the saver writes (including the ``links`` key older networkx releases use) and the docstring lists
  exactly the formats that load.

The centrality oracle is ``networkx.degree_centrality`` / ``networkx.betweenness_centrality`` computed independently
of the code under test; the round trip compares against the graph that was saved. No test doubles.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import networkx as nx
import pytest

from metainformant.networks.interaction._ppi_impl import (
    ProteinNetwork,
    load_ppi_network,
    ppi_network_analysis,
    save_ppi_network,
)


def _hub_graph() -> nx.Graph:
    """30 proteins; the three hubs are inserted last, so they are not among the first ten nodes."""
    graph = nx.Graph()
    graph.add_nodes_from(f"n{i:02d}" for i in range(30))
    graph.add_edges_from((f"n{i:02d}", f"n{i + 1:02d}") for i in range(19))  # a path through n00..n19
    for hub, start in (("n27", 0), ("n28", 7), ("n29", 14)):
        graph.add_edges_from((hub, f"n{(start + j) % 27:02d}") for j in range(12))
    return graph


@pytest.mark.parametrize(
    ("key", "oracle"),
    [("degree_centrality", nx.degree_centrality), ("betweenness_centrality", nx.betweenness_centrality)],
)
def test_reported_centrality_is_the_ten_most_central_nodes(key: str, oracle) -> None:
    graph = _hub_graph()
    expected = oracle(graph)
    first_ten_inserted = {f"n{i:02d}" for i in range(10)}

    reported = ppi_network_analysis(graph)["centrality"][key]

    assert len(reported) == 10
    # values are the true centralities, listed from most to least central
    assert all(reported[node] == pytest.approx(expected[node]) for node in reported)
    values = list(reported.values())
    assert values == sorted(values, reverse=True)
    # and they are the ten largest centralities in the whole graph (ties at the cut-off may pick either node)
    assert values == pytest.approx(sorted(expected.values(), reverse=True)[:10])
    assert {"n27", "n28", "n29"} <= set(reported), "the hubs must be reported"
    assert set(reported) != first_ten_inserted


def test_reported_centrality_is_stable_when_centralities_tie() -> None:
    # every node has the same degree and the nodes are inserted from p14 down to p00
    graph = nx.cycle_graph([f"p{i:02d}" for i in reversed(range(15))])

    first = ppi_network_analysis(graph)["centrality"]["degree_centrality"]
    second = ppi_network_analysis(graph.copy())["centrality"]["degree_centrality"]

    # ties are broken by node name, not by insertion order
    assert list(first) == list(second) == [f"p{i:02d}" for i in range(10)]


def _protein_graph() -> nx.Graph:
    graph = nx.Graph()
    graph.add_edge("P_1", "P_2", weight=0.9, evidence="two-hybrid")
    graph.add_edge("P_2", "P_3", weight=0.4)
    graph.add_node("LONELY", organism="human")
    return graph


def test_saved_json_network_loads_back_unchanged(tmp_path: Path) -> None:
    graph = _protein_graph()
    path = tmp_path / "ppi.json"
    save_ppi_network(graph, path, format="json")

    loaded = load_ppi_network(path, format="json")

    assert isinstance(loaded, nx.Graph) and not loaded.is_directed()
    assert nx.utils.graphs_equal(loaded, graph)
    assert loaded["P_1"]["P_2"] == {"weight": 0.9, "evidence": "two-hybrid"}
    assert loaded.nodes["LONELY"] == {"organism": "human"}


def test_protein_network_from_file_reads_json(tmp_path: Path) -> None:
    path = tmp_path / "ppi.json"
    save_ppi_network(_protein_graph(), path, format="json")

    network = ProteinNetwork.from_file(path, format="json")

    assert set(network.get_proteins()) == {"P_1", "P_2", "P_3", "LONELY"}
    assert {frozenset(pair) for pair in network.get_interactions()} == {
        frozenset(("P_1", "P_2")),
        frozenset(("P_2", "P_3")),
    }
    assert network.metadata["format"] == "json"


def test_json_with_the_legacy_links_key_loads(tmp_path: Path) -> None:
    path = tmp_path / "legacy.json"
    path.write_text(
        json.dumps(
            {
                "directed": False,
                "multigraph": False,
                "graph": {},
                "nodes": [{"id": "A"}, {"id": "B"}],
                "links": [{"source": "A", "target": "B", "weight": 0.5}],
            }
        )
    )

    loaded = load_ppi_network(path, format="json")

    assert loaded["A"]["B"]["weight"] == 0.5


def test_json_edges_without_a_score_get_the_default_weight(tmp_path: Path) -> None:
    path = tmp_path / "unscored.json"
    path.write_text(json.dumps({"nodes": [{"id": "A"}, {"id": "B"}], "edges": [{"source": "A", "target": "B"}]}))

    assert load_ppi_network(path, format="json")["A"]["B"]["weight"] == 1.0


@pytest.mark.parametrize(
    "payload",
    [
        [{"source": "A", "target": "B"}],  # not a node-link object
        {"nodes": [{"name": "A"}], "edges": []},  # node without an id
        {"nodes": [{"id": "A"}], "edges": [{"source": "A"}]},  # edge without a target
    ],
)
def test_malformed_json_is_rejected_with_a_value_error(tmp_path: Path, payload: object) -> None:
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="PPI JSON"):
        load_ppi_network(path, format="json")


def test_from_file_docstring_only_lists_formats_the_loader_accepts(tmp_path: Path) -> None:
    match = re.search(r"format: File format \(([^)]*)\)", ProteinNetwork.from_file.__doc__ or "")
    assert match, "ProteinNetwork.from_file no longer documents its formats in the expected form"
    advertised = re.findall(r"'(\w+)'", match.group(1))
    assert "json" in advertised

    empty = tmp_path / "empty.txt"
    empty.write_text("")
    for fmt in advertised:
        try:
            load_ppi_network(empty, format=fmt)
        except ValueError as error:
            assert "Unsupported PPI format" not in str(error), f"{fmt!r} is documented but not loadable"
        except Exception:  # noqa: BLE001 - an empty file may be invalid for the format; only the format name matters
            pass
