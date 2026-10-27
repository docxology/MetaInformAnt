"""Regression tests for feed-forward-loop detection and the advertised GRN inference methods.

* A feed-forward loop is three edges, ``TF1 -> TF2``, ``TF2 -> Target`` and the direct ``TF1 -> Target``.
  ``regulatory_motifs`` (and ``GeneRegulatoryNetwork.regulatory_motifs``, which reports the same motif) only checked
  the first two, so every plain cascade ``TF1 -> TF2 -> Target`` was reported as a feed-forward loop.
* ``infer_grn`` advertised ``'granger'`` and ``infer_regulatory_network_from_expression`` advertised ``'grnboost2'``;
  both raised when called. Granger causality needs time-ordered observations, which a samples x genes matrix does not
  guarantee, so it is not offered (time series go through ``information_flow_network``); GRNBoost2 needs the
  ``arboreto`` package. The documented method lists now match what the functions run.

The motif oracle is a brute-force count over node triples of an independent ``networkx`` graph. The documentation
check parses the method list out of each docstring and runs every advertised method. No test doubles.
"""

from __future__ import annotations

import itertools
import random
import re

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from metainformant.networks.interaction.regulatory_analysis import (
    infer_grn,
    infer_regulatory_network_from_expression,
    regulatory_motifs,
)
from metainformant.networks.interaction.regulatory_core import GeneRegulatoryNetwork


def _grn(edges: list[tuple[str, str]]) -> GeneRegulatoryNetwork:
    grn = GeneRegulatoryNetwork("ffl")
    for tf, target in edges:
        grn.add_regulation(tf, target)
    return grn


def _brute_force_ffl_triples(edges: list[tuple[str, str]]) -> set[tuple[str, str, str]]:
    """All (a, b, c), pairwise distinct, with a -> b, b -> c and a -> c, straight from the definition."""
    graph = nx.DiGraph(edges)
    return {
        (a, b, c)
        for a, b, c in itertools.permutations(graph.nodes, 3)
        if graph.has_edge(a, b) and graph.has_edge(b, c) and graph.has_edge(a, c)
    }


def test_a_plain_cascade_is_not_a_feed_forward_loop() -> None:
    grn = _grn([("TF1", "TF2"), ("TF2", "T")])

    assert [m for m in regulatory_motifs(grn) if m["motif_type"] == "feed_forward_loop"] == []
    assert grn.regulatory_motifs()["feed_forward"] == []


def test_adding_the_direct_edge_completes_the_loop() -> None:
    grn = _grn([("TF1", "TF2"), ("TF2", "T"), ("TF1", "T")])

    found = [m["genes"] for m in regulatory_motifs(grn) if m["motif_type"] == "feed_forward_loop"]

    assert found == [["TF1", "TF2", "T"]]
    assert grn.regulatory_motifs()["feed_forward"] == [["TF1", "TF2", "T"]]


def test_motifs_match_a_brute_force_triple_count_on_random_networks() -> None:
    rng = random.Random(20261027)
    nodes = [f"N{i}" for i in range(9)]
    checked_with_loops = 0
    for _ in range(40):
        edges = [(a, b) for a, b in itertools.permutations(nodes, 2) if rng.random() < 0.25]
        expected = _brute_force_ffl_triples(edges)
        grn = _grn(edges)

        standalone = {tuple(m["genes"]) for m in regulatory_motifs(grn) if m["motif_type"] == "feed_forward_loop"}
        from_class = {tuple(triple) for triple in grn.regulatory_motifs()["feed_forward"]}

        assert standalone == expected
        assert from_class == expected
        checked_with_loops += bool(expected)
    assert checked_with_loops >= 10, "the random networks must actually contain feed-forward loops"


def _advertised_methods(function) -> list[str]:
    """Method names listed in the ``method:`` entry of a docstring, e.g. ``('correlation' or 'mutual_info')``."""
    match = re.search(r"method: Inference method \(([^)]*)\)", function.__doc__ or "")
    assert match, f"{function.__name__} no longer documents its methods in the expected form"
    return re.findall(r"'(\w+)'", match.group(1))


def test_infer_grn_runs_every_method_it_advertises() -> None:
    rng = np.random.default_rng(3)
    expression = rng.normal(size=(30, 4))
    genes = ["G0", "G1", "G2", "G3"]

    advertised = _advertised_methods(infer_grn)

    assert {"correlation", "mutual_info"} <= set(advertised)
    for method in advertised:
        network = infer_grn(expression, genes, method=method, threshold=0.2, tf_genes=["G0"])
        assert isinstance(network, GeneRegulatoryNetwork)


def test_infer_grn_granger_is_rejected_with_the_supported_methods_listed() -> None:
    expression = np.random.default_rng(3).normal(size=(30, 4))

    with pytest.raises(ValueError, match="Unknown GRN inference method: granger.*correlation.*mutual_info"):
        infer_grn(expression, ["G0", "G1", "G2", "G3"], method="granger")


def test_expression_inference_runs_every_method_it_advertises() -> None:
    rng = np.random.default_rng(3)
    frame = pd.DataFrame(rng.normal(size=(4, 30)), index=["TF1", "T1", "T2", "T3"])

    advertised = _advertised_methods(infer_regulatory_network_from_expression)

    assert {"correlation", "mutual_info"} <= set(advertised)
    for method in advertised:
        graph = infer_regulatory_network_from_expression(frame, regulators=["TF1"], method=method, threshold=0.0)
        assert isinstance(graph, nx.DiGraph)


def test_grnboost2_is_not_advertised_and_still_fails_loudly() -> None:
    frame = pd.DataFrame(np.random.default_rng(3).normal(size=(3, 20)), index=["TF1", "T1", "T2"])

    assert "grnboost2" not in _advertised_methods(infer_regulatory_network_from_expression)
    with pytest.raises(NotImplementedError, match="arboreto"):
        infer_regulatory_network_from_expression(frame, regulators=["TF1"], method="grnboost2")
