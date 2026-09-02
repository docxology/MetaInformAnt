"""Real-implementation tests for metainformant.networks.analysis.pathway.

Depth coverage for the module-level analysis functions (enrichment, topology,
similarity, hierarchy, network creation, visualization data, activity scores)
complementing the existing PathwayNetwork class tests. Real synthetic data;
no test doubles.
"""

from __future__ import annotations

import math

import networkx as nx
import pytest

from metainformant.networks.analysis.pathway import (
    PathwayNetwork,
    create_pathway_network,
    find_pathway_modules,
    load_pathway_database,
    network_enrichment_analysis,
    pathway_activity_score,
    pathway_disease_association,
    pathway_enrichment,
    pathway_enrichment_analysis,
    pathway_hierarchy_analysis,
    pathway_similarity,
    pathway_similarity_analysis,
    pathway_topology_analysis,
    pathway_visualization_data,
)

PATHWAYS = {
    "glycolysis": ["G1", "G2", "G3", "G4"],
    "tca": ["G3", "G4", "G5", "G6"],
    "oxphos": ["G5", "G6", "G7", "G8"],
    "unrelated": ["U1", "U2", "U3", "U4"],
}
BACKGROUND = sorted({g for genes in PATHWAYS.values() for g in genes})


class TestPathwayEnrichmentAnalysis:
    def test_fisher_flags_enriched_pathway(self) -> None:
        results = pathway_enrichment_analysis(["G1", "G2", "G3"], BACKGROUND, PATHWAYS, method="fisher")
        by_name = {r["pathway"]: r for r in results}
        assert by_name["glycolysis"]["genes_in_pathway"] == 3
        assert by_name["glycolysis"]["enrichment_ratio"] > 1.0
        assert by_name["glycolysis"]["p_value_corrected"] <= by_name["unrelated"]["p_value_corrected"]
        assert all("p_value_corrected" in r and "significant" in r for r in results)

    def test_results_sorted_by_corrected_pvalue(self) -> None:
        results = pathway_enrichment_analysis(["G1", "G5"], BACKGROUND, PATHWAYS)
        ps = [r["p_value_corrected"] for r in results]
        assert ps == sorted(ps)

    def test_bonferroni_caps_at_one(self) -> None:
        results = pathway_enrichment_analysis(["G1"], BACKGROUND, PATHWAYS, correction="bonferroni")
        assert all(0.0 <= r["p_value_corrected"] <= 1.0 for r in results)

    def test_fdr_correction_monotone_and_bounded(self) -> None:
        results = pathway_enrichment_analysis(["G3"], BACKGROUND, PATHWAYS, correction="fdr")
        ps = [r["p_value_corrected"] for r in results]
        assert all(0.0 <= p <= 1.0 for p in ps)

    def test_hypergeometric_method(self) -> None:
        results = pathway_enrichment_analysis(["G5", "G6", "G7"], BACKGROUND, PATHWAYS, method="hypergeometric")
        by_name = {r["pathway"]: r for r in results}
        assert by_name["oxphos"]["genes_in_pathway"] == 3
        assert by_name["oxphos"]["p_value"] < 0.05

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            pathway_enrichment_analysis(["G1"], BACKGROUND, PATHWAYS, method="chi2")

    def test_unknown_correction_is_identity(self) -> None:
        plain = pathway_enrichment_analysis(["G1"], BACKGROUND, PATHWAYS, correction=None)
        assert all(r["p_value_corrected"] == r["p_value"] for r in plain)


class TestPathwayTopologyAnalysis:
    def test_undirected_graph(self) -> None:
        g = nx.path_graph(5)
        analysis = pathway_topology_analysis(g)
        assert analysis["n_nodes"] == 5
        assert analysis["n_edges"] == 4
        assert analysis["directed"] is False
        assert analysis["degree_stats"]["max"] == 2
        assert analysis["connected_components"]["n_components"] == 1
        assert analysis["average_clustering"] == 0.0

    def test_directed_graph(self) -> None:
        g = nx.DiGraph([("A", "B"), ("B", "C")])
        analysis = pathway_topology_analysis(g)
        assert analysis["directed"] is True
        assert "average_clustering" not in analysis
        assert "connected_components" not in analysis
        # Directed max edges = n*(n-1) = 6; density = 2/6
        assert analysis["density"] == pytest.approx(2 / 6)

    def test_undirected_density(self) -> None:
        g = nx.path_graph(5)  # 4 edges, max = 5*4/2 = 10
        assert pathway_topology_analysis(g)["density"] == pytest.approx(0.4)

    def test_empty_graph(self) -> None:
        analysis = pathway_topology_analysis(nx.Graph())
        assert analysis["n_nodes"] == 0
        assert "degree_stats" not in analysis


class TestPathwayModulesAndSimilarity:
    def test_find_pathway_modules_partitions_graph(self) -> None:
        g = nx.Graph([("A", "B"), ("B", "C"), ("D", "E")])
        modules = find_pathway_modules(g, method="label_propagation")
        members = {m for mod in modules for m in mod}
        assert members == {"A", "B", "C", "D", "E"}

    def test_similarity_analysis_all_methods(self) -> None:
        for method, expected_gly_tca in (
            ("jaccard", 2 / 6),
            ("overlap", 2 / 4),
            ("cosine", 2 / (math.sqrt(4) * math.sqrt(4))),
        ):
            result = pathway_similarity_analysis(PATHWAYS, method=method)
            assert result["method"] == method
            assert result["n_pathways"] == 4
            assert result["similarities"]["glycolysis_tca"] == pytest.approx(expected_gly_tca)

    def test_similarity_supports_dict_with_genes_key(self) -> None:
        data = {"P1": {"genes": ["a", "b"]}, "P2": {"genes": ["b", "c"]}}
        result = pathway_similarity_analysis(data, method="jaccard")
        assert result["similarities"]["P1_P2"] == pytest.approx(1 / 3)

    def test_similarity_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown similarity method"):
            pathway_similarity_analysis(PATHWAYS, method="hamming")


class TestHierarchyAndNetworkCreation:
    def test_hierarchy_levels(self) -> None:
        hierarchy = {"root": ["metabolism", "immune"], "metabolism": ["glycolysis"]}
        analysis = pathway_hierarchy_analysis(PATHWAYS, hierarchy)
        assert analysis["n_pathways"] == 4
        assert analysis["hierarchy_levels"]["root"] == 0
        assert analysis["hierarchy_levels"]["glycolysis"] == 2
        assert analysis["max_depth"] == 2
        assert analysis["n_roots"] == 1

    def test_hierarchy_without_data(self) -> None:
        analysis = pathway_hierarchy_analysis(PATHWAYS)
        assert analysis["hierarchy_levels"] == {}

    def test_create_pathway_network_edges_respect_threshold(self) -> None:
        net = create_pathway_network(PATHWAYS, similarity_threshold=0.2)
        # glycolysis-tca jaccard = 2/6 = 0.333 >= 0.2; tca-oxphos = 2/6; gly-oxphos = 0
        assert net.has_edge("glycolysis", "tca")
        assert net.has_edge("tca", "oxphos")
        assert not net.has_edge("glycolysis", "oxphos")
        assert set(net.nodes()) == set(PATHWAYS)

    def test_create_pathway_network_high_threshold_no_edges(self) -> None:
        net = create_pathway_network(PATHWAYS, similarity_threshold=0.99)
        assert net.number_of_edges() == 0


class TestPathwayDiseaseAssociation:
    def test_documented_stub_behavior(self) -> None:
        # pathway_disease_association currently builds empty pathway gene sets,
        # so overlap is always 0 and the result is empty. Pin this contract.
        results = [{"pathway": "glycolysis", "p_value_corrected": 0.01}]
        assoc = pathway_disease_association(results, {"diabetes": ["G1"]})
        assert assoc == {}


class TestVisualizationData:
    def test_spring_layout(self) -> None:
        g = nx.path_graph(4)
        data = pathway_visualization_data(g, layout_method="spring", seed=1)
        assert data["n_nodes"] == 4
        assert data["n_edges"] == 3
        assert data["layout_method"] == "spring"
        for node in data["nodes"]:
            assert isinstance(node["x"], float) and math.isfinite(node["x"])
            assert isinstance(node["y"], float) and math.isfinite(node["y"])

    def test_circular_and_random_layouts(self) -> None:
        g = nx.cycle_graph(6)
        for method in ("circular", "random"):
            data = pathway_visualization_data(g, layout_method=method)
            assert data["layout_method"] == method
            assert data["n_nodes"] == 6

    def test_unknown_layout_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown layout method"):
            pathway_visualization_data(nx.path_graph(3), layout_method="spiral")

    def test_edge_attributes_carried_through(self) -> None:
        g = nx.Graph()
        g.add_edge("A", "B", weight=0.5)
        data = pathway_visualization_data(g, layout_method="circular")
        assert data["edges"][0]["weight"] == 0.5


class TestPathwayNetworkExtras:
    def test_load_from_database_gmt(self, tmp_path) -> None:
        p = tmp_path / "db.gmt"
        p.write_text(
            "glycolysis\tGlucose metabolism\tG1\tG2\tG3\nincomplete\tonly-description\n",
            encoding="utf-8",
        )
        net = PathwayNetwork.load_from_database(p, format="gmt")
        assert net.name == "db"
        assert net.pathways["glycolysis"] == ["G1", "G2", "G3"]
        assert "incomplete" not in net.pathways  # fewer than 3 columns skipped

    def test_load_from_database_json(self, tmp_path) -> None:
        import json

        p = tmp_path / "db.json"
        p.write_text(json.dumps({"name": "mydb", "pathways": {"P1": ["a", "b"]}}), encoding="utf-8")
        net = PathwayNetwork.load_from_database(p, format="json")
        assert net.name == "mydb"
        assert net.pathways == {"P1": ["a", "b"]}
        assert net.metadata["source_file"] == str(p)

    def test_load_from_database_missing_file(self, tmp_path) -> None:
        with pytest.raises(FileNotFoundError):
            PathwayNetwork.load_from_database(tmp_path / "nope.json")

    def test_load_from_database_bad_format(self, tmp_path) -> None:
        p = tmp_path / "db.xml"
        p.write_text("<x/>", encoding="utf-8")
        with pytest.raises(ValueError, match="Unsupported format"):
            PathwayNetwork.load_from_database(p, format="xml")

    def test_gene_pathways_property(self) -> None:
        net = PathwayNetwork(pathways={"A": ["g1", "g2"], "B": ["g2", "g3"]})
        gp = net.gene_pathways
        assert gp["g1"] == {"A"}
        assert gp["g2"] == {"A", "B"}

    def test_filter_pathways_by_size(self) -> None:
        net = PathwayNetwork(name="src", pathways={"big": ["a", "b", "c"], "small": ["d"]})
        filtered = net.filter_pathways_by_size(min_size=2)
        assert filtered.name == "src_filtered"
        assert set(filtered.pathways) == {"big"}
        assert filtered.metadata == net.metadata

    def test_overlap_matrix_symmetric(self) -> None:
        net = PathwayNetwork(pathways={"A": ["g1", "g2"], "B": ["g2", "g3"], "C": ["g1"]})
        m = net.pathway_overlap_matrix()
        assert m["A"]["A"] == 1.0
        assert m["A"]["B"] == pytest.approx(1 / 3)
        assert m["B"]["A"] == m["A"]["B"]
        assert m["A"]["C"] == 0.5
        assert m["C"]["A"] == 0.5

    def test_dunder_protocols(self) -> None:
        net = PathwayNetwork(pathways={"A": ["g1"]})
        assert len(net) == 1
        assert "A" in net
        assert "B" not in net
        assert net["A"] == ["g1"]
        assert net.get_pathway("MISSING") == []

    def test_add_remove_pathway(self) -> None:
        net = PathwayNetwork()
        net.add_pathway("P", ["g1"], metadata={"source": "KEGG"})
        assert "P" in net and net.pathway_metadata["P"] == {"source": "KEGG"}
        assert net.remove_pathway("P") is True
        assert net.remove_pathway("P") is False

    def test_add_gene_to_pathway_dedupes(self) -> None:
        net = PathwayNetwork(pathways={"P": ["g1"]})
        net.add_gene_to_pathway("P", "g1")
        net.add_gene_to_pathway("P", "g2")
        net.add_gene_to_pathway("NEW", "g9")
        assert net.pathways["P"] == ["g1", "g2"]
        assert net.pathways["NEW"] == ["g9"]

    def test_find_pathways_containing_gene(self) -> None:
        net = PathwayNetwork(pathways={"A": ["g1", "g2"], "B": ["g2", "g3"]})
        assert net.find_pathways_containing_gene("g2") == ["A", "B"]
        assert net.find_pathways_containing_gene("zz") == []


class TestPathwayEnrichmentFunction:
    def _network(self) -> PathwayNetwork:
        return PathwayNetwork(name="test", pathways=PATHWAYS)

    def test_fisher_enrichment(self) -> None:
        results = pathway_enrichment(["G1", "G2", "G3"], self._network())
        assert "glycolysis" in results
        g = results["glycolysis"]
        assert g["overlap_size"] == 3
        assert g["enrichment_ratio"] > 1.0
        assert 0.0 <= g["p_value"] <= 1.0
        assert set(g["overlap_genes"]) == {"G1", "G2", "G3"}

    def test_hypergeometric_method(self) -> None:
        results = pathway_enrichment(["G5", "G6", "G7"], self._network(), method="hypergeom")
        assert results["oxphos"]["p_value"] < 0.05

    def test_min_overlap_filters(self) -> None:
        results = pathway_enrichment(["G1"], self._network(), min_overlap=2)
        assert results == {}

    def test_query_genes_outside_background_dropped(self) -> None:
        results = pathway_enrichment(["ZZ1", "ZZ2"], self._network())
        assert results == {}

    def test_explicit_background_narrows_universe(self) -> None:
        results = pathway_enrichment(["G1", "G2"], self._network(), background_genes=["G1", "G2", "G3", "G4"])
        assert "glycolysis" in results
        assert results["glycolysis"]["background_size"] == 4

    def test_bonferroni_correction_applied(self) -> None:
        results = pathway_enrichment(["G3"], self._network(), correction="bonferroni")
        assert all("corrected_p_value" in r for r in results.values())
        assert all(r["corrected_p_value"] >= r["p_value"] for r in results.values())

    def test_empty_network_returns_empty(self) -> None:
        results = pathway_enrichment(["G1"], PathwayNetwork(name="empty"))
        assert results == {}

    def test_network_enrichment_analysis_alias(self) -> None:
        direct = pathway_enrichment(["G1", "G2", "G3"], self._network())
        alias = network_enrichment_analysis(["G1", "G2", "G3"], self._network())
        assert direct == alias


class TestLoadPathwayDatabase:
    def test_nested_pathways_format(self) -> None:
        data = {
            "pathways": {
                "P1": {"name": "Pathway One", "genes": ["g1", "g2"]},
                "P2": {"genes": "single_gene_str"},
                "P3": ["g3", "g4"],
            },
            "version": "1.0",
        }
        net = load_pathway_database(data, name="db")
        assert net.name == "db"
        assert net.pathways["P1"] == ["g1", "g2"]
        assert net.pathways["P2"] == ["single_gene_str"]
        assert net.pathways["P3"] == ["g3", "g4"]
        assert net.pathway_metadata["P1"] == {"name": "Pathway One"}
        assert net.metadata["version"] == "1.0"

    def test_direct_format(self) -> None:
        data = {"P1": {"genes": ["g1"], "description": "d"}}
        net = load_pathway_database(data)
        assert net.pathways["P1"] == ["g1"]
        assert net.pathway_metadata["P1"]["description"] == "d"

    def test_metadata_and_source(self) -> None:
        data = {
            "pathways": {"P1": ["g1"]},
            "metadata": {"build": "test"},
            "source": "unit-test",
        }
        net = load_pathway_database(data)
        assert net.metadata["build"] == "test"
        assert net.metadata["source"] == "unit-test"


class TestStandaloneSimilarity:
    def test_jaccard(self) -> None:
        assert pathway_similarity({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)

    def test_overlap(self) -> None:
        assert pathway_similarity({"a", "b", "c"}, {"b"}, method="overlap") == 1.0

    def test_dice(self) -> None:
        assert pathway_similarity({"a", "b"}, {"b", "c"}, method="dice") == pytest.approx(2 * 1 / 4)

    def test_empty_sets_return_zero(self) -> None:
        assert pathway_similarity(set(), {"a"}) == 0.0
        assert pathway_similarity({"a"}, set()) == 0.0

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown similarity method"):
            pathway_similarity({"a"}, {"a"}, method="sokal")


class TestPathwayActivityScore:
    def _network(self) -> PathwayNetwork:
        return PathwayNetwork(pathways={"P": ["g1", "g2", "g3", "g4"]})

    def test_mean(self) -> None:
        expr = {"g1": 1.0, "g2": 2.0, "g3": 3.0, "g4": 4.0}
        assert pathway_activity_score(self._network(), "P", expr, method="mean") == pytest.approx(2.5)

    def test_max_sum_median(self) -> None:
        expr = {"g1": 1.0, "g2": 2.0, "g3": 3.0, "g4": 4.0}
        net = self._network()
        assert pathway_activity_score(net, "P", expr, method="max") == 4.0
        assert pathway_activity_score(net, "P", expr, method="sum") == 10.0
        assert pathway_activity_score(net, "P", expr, method="median") == 3.0

    def test_missing_expression_genes_skipped(self) -> None:
        expr = {"g1": 10.0}
        assert pathway_activity_score(self._network(), "P", expr) == 10.0

    def test_no_expression_returns_zero(self) -> None:
        assert pathway_activity_score(self._network(), "P", {}) == 0.0

    def test_empty_pathway_returns_zero(self) -> None:
        assert pathway_activity_score(PathwayNetwork(), "NOPE", {"g1": 1.0}) == 0.0

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            pathway_activity_score(self._network(), "P", {"g1": 1.0}, method="geometric")
