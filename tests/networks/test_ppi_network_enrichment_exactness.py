"""``ppi_network_enrichment`` keeps tiny p-values and uses the interaction density of its background.

* The p-value was ``1 - poisson.cdf(observed - 1, expected)``. The subtraction cancels to exactly 0.0 for every upper
  tail below about 1e-16, so a 12-protein clique in a sparse network got p = 0.0 instead of its real 1e-186. The
  survival function is used now.
* The expected number of interactions used ``len(ppi_graph.edges())``, every edge of the graph, as the interactions
  of the background, even when an explicit background held only part of the network; and proteins outside the
  background were still counted as part of the tested draw.

The oracle is ``scipy.stats.poisson.sf`` with the expected count recomputed from the subgraphs in the test. No test
doubles.
"""

from __future__ import annotations

import math

import networkx as nx
import pytest

stats = pytest.importorskip("scipy.stats")

from metainformant.networks.interaction.ppi import ppi_network_enrichment  # noqa: E402


def _graph_with_clique() -> nx.Graph:
    graph = nx.gnp_random_graph(200, 0.02, seed=3)
    graph = nx.relabel_nodes(graph, {n: f"P{n}" for n in graph})
    clique = [f"C{i}" for i in range(12)]
    graph.add_edges_from((a, b) for i, a in enumerate(clique) for b in clique[i + 1 :])
    return graph


def test_a_dense_clique_keeps_its_tiny_p_value() -> None:
    graph = _graph_with_clique()
    clique = [f"C{i}" for i in range(12)]

    result = ppi_network_enrichment(clique, graph)

    observed = graph.subgraph(clique).number_of_edges()
    n_nodes = graph.number_of_nodes()
    expected = graph.number_of_edges() * (12 * 11 / 2) / (n_nodes * (n_nodes - 1) / 2)
    reference = float(stats.poisson.sf(observed - 1, expected))
    assert observed == 66
    assert 0.0 < reference < 1e-30
    assert result["p_value"] == pytest.approx(reference, rel=1e-9, abs=1e-300)  # the old value was 0.0
    assert result["observed_interactions"] == 66
    assert result["expected_interactions"] == pytest.approx(expected)


def test_an_explicit_background_uses_its_own_interaction_density() -> None:
    graph = _graph_with_clique()
    background = [n for n in graph if str(n).startswith("P")][:100] + [f"C{i}" for i in range(12)]
    test = [f"C{i}" for i in range(5)] + ["P0", "P1"]

    result = ppi_network_enrichment(test, graph, background_proteins=background)

    background_edges = graph.subgraph(background).number_of_edges()
    pairs = len(background) * (len(background) - 1) / 2
    drawn = [p for p in test if p in set(background)]
    expected = background_edges * (len(drawn) * (len(drawn) - 1) / 2) / pairs
    assert result["background_proteins"] == len(background)
    assert result["expected_interactions"] == pytest.approx(expected)
    # the old numerator was every edge of the graph, which is more than the background's
    assert graph.number_of_edges() > background_edges


def test_test_proteins_outside_the_background_are_not_drawn() -> None:
    graph = nx.path_graph(10)

    result = ppi_network_enrichment([0, 1, 2, 9], graph, background_proteins=[0, 1, 2, 3, 4])

    assert result["test_proteins"] == 3  # 9 is in the network but not in the background
    assert result["observed_interactions"] == 2


def test_duplicates_count_once() -> None:
    graph = nx.complete_graph(6)

    once = ppi_network_enrichment([0, 1, 2], graph)
    repeated = ppi_network_enrichment([0, 1, 2, 2, 1, 0], graph, background_proteins=list(range(6)) * 2)

    assert once == repeated


def test_no_enrichment_gives_a_p_value_near_one_and_an_empty_test_set_is_still_reported() -> None:
    graph = nx.cycle_graph(30)

    result = ppi_network_enrichment([0, 10, 20], graph)

    assert result["observed_interactions"] == 0
    assert result["p_value"] == pytest.approx(1.0)
    assert math.isfinite(result["enrichment_ratio"])
    assert ppi_network_enrichment(["x"], graph) == {"error": "No test proteins found in PPI network"}
