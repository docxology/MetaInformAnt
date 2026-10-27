"""Regression tests: ontology network layouts take a seed, and the ontology workflow can seed its GSEA permutations.

* ``plot_go_dag``, ``plot_ontology_network`` and ``create_interactive_go_network`` laid their graphs out with
  ``nx.spring_layout`` (and, when that failed, ``nx.random_layout``) without a seed, so the same graph was drawn
  differently on every call and the call consumed NumPy's global random stream. Each now takes ``layout_seed`` and
  hands it to networkx, which builds a local generator from an int.
* ``run_ontology_analysis`` called ``gsea`` with no way to seed its permutations. The ``ontology`` section of the
  config may now carry ``random_state``, which is passed on to ``gsea``.

The layouts are compared through the node positions matplotlib actually draws. The workflow cannot run offline (it
queries QuickGO), so its seed plumbing is checked on the parsed source plus the config parser. The
``create_interactive_go_network`` layout needs plotly, which is not installed everywhere, so it is covered by the
source check only. No test doubles.
"""

from __future__ import annotations

import ast
import inspect
import textwrap

import matplotlib
import networkx as nx
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from metainformant.ontology.pathway_enrichment.enrichment import gsea  # noqa: E402
from metainformant.ontology.visualization import visualization as visualization_module  # noqa: E402
from metainformant.ontology.visualization.visualization import plot_go_dag, plot_ontology_network  # noqa: E402
from metainformant.ontology.workflow import run_ontology  # noqa: E402


def _dag(weight: object = None) -> nx.DiGraph:
    graph = nx.DiGraph()
    edges = [("root", "a"), ("root", "b"), ("a", "c"), ("a", "d"), ("b", "d"), ("b", "e"), ("c", "f"), ("e", "f")]
    for source, target in edges:
        if weight is None:
            graph.add_edge(source, target)
        else:
            graph.add_edge(source, target, weight=weight)
    return graph


def _node_positions(draw, graph: nx.Graph, **kwargs) -> np.ndarray:
    ax = draw(graph, **kwargs)
    try:
        return np.array(ax.collections[0].get_offsets())  # the node markers are the first collection drawn
    finally:
        plt.close(ax.figure)


@pytest.mark.parametrize("draw", [plot_go_dag, plot_ontology_network], ids=["plot_go_dag", "plot_ontology_network"])
def test_layout_is_reproducible_for_a_seed_and_differs_between_seeds(draw) -> None:
    first = _node_positions(draw, _dag(), layout_seed=7)

    assert first.shape == (7, 2)
    assert np.array_equal(first, _node_positions(draw, _dag(), layout_seed=7))
    assert not np.array_equal(first, _node_positions(draw, _dag(), layout_seed=8))


@pytest.mark.parametrize("draw", [plot_go_dag, plot_ontology_network], ids=["plot_go_dag", "plot_ontology_network"])
def test_a_seeded_layout_leaves_the_global_numpy_state_alone(draw) -> None:
    np.random.seed(20261027)
    before = np.random.get_state()

    _node_positions(draw, _dag(), layout_seed=3)

    after = np.random.get_state()
    assert after[0] == before[0] and after[2:] == before[2:] and np.array_equal(after[1], before[1])


def test_the_random_layout_fallback_of_plot_go_dag_is_seeded_too() -> None:
    # non-numeric edge weights make spring_layout raise ValueError, which sends plot_go_dag to nx.random_layout
    graph = _dag(weight="not-a-number")
    with pytest.raises(ValueError):
        nx.spring_layout(graph, k=1, iterations=50)

    first = _node_positions(plot_go_dag, graph, layout_seed=5)

    assert np.array_equal(first, _node_positions(plot_go_dag, graph, layout_seed=5))
    assert not np.array_equal(first, _node_positions(plot_go_dag, graph, layout_seed=6))


def _layout_calls(source: str) -> list[ast.Call]:
    tree = ast.parse(textwrap.dedent(source))
    return [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in {"spring_layout", "random_layout"}
    ]


def test_every_layout_call_in_the_visualization_module_is_seeded() -> None:
    calls = _layout_calls(inspect.getsource(visualization_module))

    assert len(calls) == 4, "plot_go_dag (two), plot_ontology_network and create_interactive_go_network"
    unseeded = [ast.unparse(call) for call in calls if "seed" not in {keyword.arg for keyword in call.keywords}]
    assert unseeded == []


def test_the_workflow_reads_the_seed_from_the_ontology_config() -> None:
    assert run_ontology._gsea_random_state({}) is None
    assert run_ontology._gsea_random_state({"random_state": None}) is None
    assert run_ontology._gsea_random_state({"random_state": 7}) == 7
    assert run_ontology._gsea_random_state({"random_state": "11"}) == 11  # YAML or env values may arrive as text
    assert run_ontology._gsea_random_state({"random_state": 0}) == 0


def test_the_workflow_hands_the_seed_to_gsea() -> None:
    tree = ast.parse(textwrap.dedent(inspect.getsource(run_ontology.run_ontology_analysis)))
    gsea_calls = [
        node for node in ast.walk(tree) if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "gsea"
    ]

    assert len(gsea_calls) == 1
    assert "random_state" in {keyword.arg for keyword in gsea_calls[0].keywords}


def test_a_configured_seed_makes_the_gsea_permutations_reproducible() -> None:
    names = [f"G{i:03d}" for i in range(80)]
    ranked = [(name, 4.0 - 0.1 * i) for i, name in enumerate(names)]
    gene_sets = {"S1": names[:12], "S2": names[30:45], "S3": names[60:75]}
    seed = run_ontology._gsea_random_state({"random_state": "21"})

    first = gsea(ranked, gene_sets, n_permutations=50, min_size=5, random_state=seed)
    second = gsea(ranked, gene_sets, n_permutations=50, min_size=5, random_state=seed)

    assert first == second
