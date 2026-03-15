"""Tests for ontology visualization: GO DAG, similarity matrices, enrichment plots.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import pytest

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

import numpy as np

from metainformant.ontology.visualization.visualization import (
    create_interactive_go_network,
    plot_functional_annotation_heatmap,
    plot_go_dag,
    plot_go_enrichment_barplot,
    plot_go_enrichment_dotplot,
    plot_go_term_hierarchy,
    plot_information_content_profile,
    plot_ontology_network,
    plot_semantic_similarity_clustermap,
    plot_semantic_similarity_matrix,
)


def _make_go_graph():
    if not HAS_NETWORKX:
        pytest.skip("networkx required")
    G = nx.DiGraph()
    G.add_node("GO:0008150", name="biological_process")
    G.add_node("GO:0009987", name="cellular_process")
    G.add_node("GO:0006950", name="response_to_stress")
    G.add_node("GO:0007154", name="cell_communication")
    G.add_edge("GO:0009987", "GO:0008150", relation="is_a")
    G.add_edge("GO:0006950", "GO:0008150", relation="is_a")
    G.add_edge("GO:0007154", "GO:0009987", relation="is_a")
    return G


def _make_enrichment_results():
    return [
        {"term_id": "GO:0006950", "name": "response_to_stress", "p_value": 0.001, "adjusted_p": 0.005, "count": 15, "fold_enrichment": 2.5},
        {"term_id": "GO:0009987", "name": "cellular_process", "p_value": 0.01, "adjusted_p": 0.03, "count": 30, "fold_enrichment": 1.8},
        {"term_id": "GO:0007154", "name": "cell_communication", "p_value": 0.05, "adjusted_p": 0.08, "count": 10, "fold_enrichment": 1.5},
    ]


@pytest.fixture(autouse=True)
def close_plots():
    yield
    if HAS_MATPLOTLIB:
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_go_dag
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NETWORKX, reason="matplotlib/networkx required")
class TestPlotGODAG:
    def test_basic_dag(self):
        G = _make_go_graph()
        fig, ax = plt.subplots()
        plot_go_dag(G, ax=ax)
        assert ax is not None

    def test_with_specific_terms(self):
        G = _make_go_graph()
        fig, ax = plt.subplots()
        plot_go_dag(G, terms=["GO:0009987", "GO:0008150"], ax=ax)
        assert ax is not None

    def test_output_to_file(self, tmp_path):
        G = _make_go_graph()
        fig, ax = plt.subplots()
        plot_go_dag(G, ax=ax)
        output_path = tmp_path / "go_dag.png"
        fig.savefig(output_path)
        assert output_path.exists()


# ---------------------------------------------------------------------------
# plot_semantic_similarity_matrix
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotSemanticSimilarityMatrix:
    def test_basic_matrix(self):
        matrix = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.6], [0.3, 0.6, 1.0]])
        fig, ax = plt.subplots()
        plot_semantic_similarity_matrix(matrix, ax=ax)
        assert ax is not None

    def test_with_labels(self):
        matrix = np.array([[1.0, 0.8], [0.8, 1.0]])
        fig, ax = plt.subplots()
        plot_semantic_similarity_matrix(matrix, term_labels=["term_a", "term_b"], ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_go_enrichment_barplot
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotGOEnrichmentBarplot:
    def test_basic_barplot(self):
        results = _make_enrichment_results()
        fig, ax = plt.subplots()
        plot_go_enrichment_barplot(results, ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_go_enrichment_dotplot
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotGOEnrichmentDotplot:
    def test_basic_dotplot(self):
        results = _make_enrichment_results()
        fig, ax = plt.subplots()
        plot_go_enrichment_dotplot(results, ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_ontology_network
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NETWORKX, reason="matplotlib/networkx required")
class TestPlotOntologyNetwork:
    def test_basic_network(self):
        G = _make_go_graph()
        fig, ax = plt.subplots()
        plot_ontology_network(G, ax=ax)
        assert ax is not None

    def test_with_node_colors(self):
        G = _make_go_graph()
        colors = {node: i * 0.25 for i, node in enumerate(G.nodes())}
        fig, ax = plt.subplots()
        plot_ontology_network(G, node_colors=colors, ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_information_content_profile
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotInformationContentProfile:
    def test_basic_ic_profile(self):
        term_ic = {"GO:0008150": 0.5, "GO:0009987": 2.0, "GO:0006950": 3.5}
        fig, ax = plt.subplots()
        plot_information_content_profile(term_ic, ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_go_term_hierarchy
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NETWORKX, reason="matplotlib/networkx required")
class TestPlotGOTermHierarchy:
    def test_basic_hierarchy(self):
        G = _make_go_graph()
        fig, ax = plt.subplots()
        plot_go_term_hierarchy(G, root_term="GO:0008150", ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_functional_annotation_heatmap
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotFunctionalAnnotationHeatmap:
    def test_basic_heatmap(self):
        matrix = np.array([[1, 0, 1], [0, 1, 0], [1, 1, 0]])
        fig, ax = plt.subplots()
        plot_functional_annotation_heatmap(matrix, term_labels=["t1", "t2", "t3"], ax=ax)
        assert ax is not None


# ---------------------------------------------------------------------------
# plot_semantic_similarity_clustermap
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotSemanticSimilarityClustermap:
    def test_basic_clustermap(self):
        try:
            import seaborn  # noqa: F401
        except ImportError:
            pytest.skip("seaborn required for clustermap")
        matrix = np.array([[1.0, 0.7, 0.2], [0.7, 1.0, 0.5], [0.2, 0.5, 1.0]])
        result = plot_semantic_similarity_clustermap(matrix, term_labels=["a", "b", "c"])
        assert result is not None


# ---------------------------------------------------------------------------
# create_interactive_go_network
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_NETWORKX, reason="networkx required")
class TestCreateInteractiveGONetwork:
    def test_basic_interactive(self):
        try:
            import plotly  # noqa: F401
        except ImportError:
            pytest.skip("plotly required for interactive GO network")
        G = _make_go_graph()
        result = create_interactive_go_network(G)
        assert result is not None
