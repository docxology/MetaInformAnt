"""Real-implementation tests for metainformant.networks.interaction._ppi_impl.

Covers file I/O (load/save across formats), module-level analysis functions
(hubs, clustering, disease associations, comparison, enrichment), the
ProteinNetwork class, and interaction prediction. All inputs are real files in
tmp_path or real synthetic graphs; no test doubles.
"""

from __future__ import annotations

import pandas as pd
import pytest

from metainformant.networks.interaction._ppi_impl import (
    ProteinNetwork,
    analyze_ppi_disease_associations,
    construct_ppi_network_from_interactions,
    detect_complexes,
    export_to_string_format,
    find_ppi_hubs,
    functional_enrichment_ppi,
    load_ppi_network,
    load_string_interactions,
    ppi_network_analysis,
    ppi_network_clustering,
    ppi_network_comparison,
    ppi_network_enrichment,
    predict_interactions,
    protein_similarity,
    save_ppi_network,
)


def _write(tmp_path, name: str, text: str):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _sample_graph():
    interactions = [
        ("A", "B", 0.9),
        ("A", "C", 0.8),
        ("A", "D", 0.7),
        ("B", "C", 0.6),
        ("E", "F", 0.5),
    ]
    return construct_ppi_network_from_interactions(interactions)


class TestLoadPpiNetwork:
    def test_tsv_with_scores(self, tmp_path) -> None:
        f = _write(tmp_path, "ppi.tsv", "#comment line\nA\tB\t0.9\nB\tC\t0.5\n")
        g = load_ppi_network(f, format="tsv")
        assert set(g.nodes()) == {"A", "B", "C"}
        assert g["A"]["B"]["weight"] == 0.9
        assert g["B"]["C"]["weight"] == 0.5

    def test_tsv_without_scores_defaults_one(self, tmp_path) -> None:
        f = _write(tmp_path, "ppi.tsv", "A\tB\n")
        g = load_ppi_network(f, format="tsv")
        assert g["A"]["B"]["weight"] == 1.0

    def test_tsv_duplicate_keeps_higher_score(self, tmp_path) -> None:
        f = _write(tmp_path, "ppi.tsv", "A\tB\t0.3\nA\tB\t0.8\n")
        g = load_ppi_network(f, format="tsv")
        assert g["A"]["B"]["weight"] == 0.8

    def test_csv_with_header(self, tmp_path) -> None:
        f = _write(tmp_path, "ppi.csv", "protein1,protein2,score\nA,B,0.7\nC,D,0.4\n")
        g = load_ppi_network(f, format="csv")
        assert set(g.nodes()) == {"A", "B", "C", "D"}
        assert g["A"]["B"]["weight"] == 0.7

    def test_bioplex_skips_header(self, tmp_path) -> None:
        f = _write(tmp_path, "bioplex.tsv", "GeneA\tGeneB\tother\nA\tB\tx\n")
        g = load_ppi_network(f, format="bioplex")
        assert set(g.nodes()) == {"A", "B"}
        assert g["A"]["B"]["weight"] == 1.0  # BioPlex carries no scores

    def test_intact_extracts_ids_from_prefixed_fields(self, tmp_path) -> None:
        f = _write(
            tmp_path,
            "intact.tsv",
            "#comment\nuniprotkb:P12345\tuniprotkb:P67890\tx\tx\tx\tx\n",
        )
        g = load_ppi_network(f, format="intact")
        assert set(g.nodes()) == {"P12345", "P67890"}

    def test_unsupported_format_raises(self, tmp_path) -> None:
        f = _write(tmp_path, "ppi.tsv", "A\tB\n")
        with pytest.raises(ValueError, match="Unsupported PPI format"):
            load_ppi_network(f, format="graphml")

    def test_roundtrip_tsv_save_load(self, tmp_path) -> None:
        g = _sample_graph()
        out = tmp_path / "out" / "saved.tsv"
        save_ppi_network(g, out, format="tsv")
        g2 = load_ppi_network(out, format="tsv")
        assert set(g2.edges()) == set(g.edges())
        assert g2["A"]["B"]["weight"] == 0.9


class TestSavePpiNetwork:
    def test_csv_save(self, tmp_path) -> None:
        g = _sample_graph()
        out = tmp_path / "saved.csv"
        save_ppi_network(g, out, format="csv")
        text = out.read_text()
        assert text.splitlines()[0] == "protein1,protein2,score"
        assert "A,B,0.9" in text

    def test_json_save(self, tmp_path) -> None:
        from metainformant.core import io as core_io

        g = _sample_graph()
        out = tmp_path / "saved.json"
        save_ppi_network(g, out, format="json")
        data = core_io.load_json(out)
        assert data["nodes"], "node-link JSON should contain nodes"
        # networkx renamed 'links' -> 'edges' across major versions
        edge_key = "links" if "links" in data else "edges"
        assert len(data[edge_key]) == len(g.edges())

    def test_unsupported_format_raises(self, tmp_path) -> None:
        g = _sample_graph()
        with pytest.raises(ValueError, match="Unsupported output format"):
            save_ppi_network(g, tmp_path / "x.parquet", format="parquet")


class TestModuleLevelAnalysis:
    def test_construct_from_interactions(self) -> None:
        g = _sample_graph()
        assert g.number_of_nodes() == 6
        assert g.number_of_edges() == 5
        assert g["A"]["B"]["weight"] == 0.9

    def test_construct_dedupes_with_higher_score(self) -> None:
        g = construct_ppi_network_from_interactions([("A", "B", 0.2), ("A", "B", 0.6)])
        assert g.number_of_edges() == 1
        assert g["A"]["B"]["weight"] == 0.6

    def test_ppi_network_analysis_basic(self) -> None:
        analysis = ppi_network_analysis(_sample_graph())
        assert analysis["basic_stats"]["n_proteins"] == 6
        assert analysis["basic_stats"]["n_interactions"] == 5
        assert analysis["degree_distribution"]["max_degree"] == 3  # A
        assert analysis["connectivity"]["n_components"] == 2
        assert analysis["connectivity"]["largest_component_size"] == 4
        assert "density" in analysis["topology"]
        assert "degree_centrality" in analysis["centrality"]

    def test_ppi_network_analysis_empty_graph(self) -> None:
        import networkx as nx

        analysis = ppi_network_analysis(nx.Graph())
        assert analysis["basic_stats"]["n_proteins"] == 0

    def test_ppi_network_analysis_skips_centrality_when_large(self) -> None:
        import networkx as nx

        g = nx.path_graph(50)
        analysis = ppi_network_analysis(g, max_nodes_for_centrality=10)
        assert analysis["centrality"]["skipped"] is True
        assert "too large" in analysis["centrality"]["reason"]

    def test_find_ppi_hubs_explicit_threshold(self) -> None:
        hubs = find_ppi_hubs(_sample_graph(), degree_threshold=3)
        assert hubs == [("A", 3)]

    def test_find_ppi_hubs_percentile(self) -> None:
        hubs = find_ppi_hubs(_sample_graph(), percentile=70.0)
        assert hubs[0][0] == "A"
        degrees = [d for _, d in hubs]
        assert degrees == sorted(degrees, reverse=True)

    def test_ppi_network_clustering(self) -> None:
        clusters = ppi_network_clustering(_sample_graph(), method="label_propagation")
        all_members = {p for cluster in clusters for p in cluster}
        assert all_members == {"A", "B", "C", "D", "E", "F"}

    def test_disease_associations(self) -> None:
        g = _sample_graph()
        assoc = analyze_ppi_disease_associations(g, {"diseaseX": ["A", "B", "Z"], "diseaseY": ["W"]})
        assert "diseaseX" in assoc
        x = assoc["diseaseX"]
        assert x["proteins_in_network"] == 2  # Z not in network
        assert x["total_proteins"] == 3
        assert x["subgraph_edges"] == 1  # A-B
        assert x["external_connections"] == 2  # C, D
        assert "diseaseY" not in assoc  # no proteins in network

    def test_network_comparison(self) -> None:
        g1 = construct_ppi_network_from_interactions([("A", "B", 1.0), ("C", "D", 1.0)])
        g2 = construct_ppi_network_from_interactions([("A", "B", 1.0), ("E", "F", 1.0)])
        comp = ppi_network_comparison(g1, g2)
        assert comp["protein_overlap"]["common_proteins"] == 2
        assert comp["protein_overlap"]["unique_to_1"] == 2
        assert comp["protein_overlap"]["jaccard_similarity"] == pytest.approx(2 / 6)
        assert comp["interaction_overlap"]["common_interactions"] == 1
        assert comp["interaction_overlap"]["jaccard_similarity"] == pytest.approx(1 / 3)

    def test_ppi_enrichment_clique_enriched(self) -> None:
        import networkx as nx

        # Background: sparse ring; test set: dense triangle
        g = nx.cycle_graph(20)
        g.add_edges_from([("T1", "T2"), ("T2", "T3"), ("T1", "T3")])
        result = ppi_network_enrichment(["T1", "T2", "T3"], g)
        assert result["test_proteins"] == 3
        assert result["observed_interactions"] == 3
        assert result["enrichment_ratio"] > 1.0
        assert 0.0 <= result["p_value"] <= 1.0

    def test_ppi_enrichment_no_test_proteins(self) -> None:
        g = _sample_graph()
        result = ppi_network_enrichment(["NOPE1", "NOPE2"], g)
        assert result == {"error": "No test proteins found in PPI network"}


class TestProteinNetworkClass:
    def test_init_empty_and_len_contains(self) -> None:
        net = ProteinNetwork(name="empty")
        assert len(net) == 0
        assert "A" not in net
        net.add_interaction("A", "B", 0.5)
        assert len(net) == 2
        assert "A" in net

    def test_init_from_string_sets_name(self) -> None:
        net = ProteinNetwork("named_graph")
        assert net.name == "named_graph"
        assert net.graph is not None

    def test_init_from_existing_graph_hydrates_interactions(self) -> None:
        import networkx as nx

        g = nx.Graph()
        g.add_edge("A", "B", weight=0.7)
        g.add_node("C", gene_name="geneC")
        net = ProteinNetwork(graph=g)
        assert len(net.interactions) == 1
        assert net.interactions[0][2]["confidence"] == 0.7
        assert net.protein_metadata["C"] == {"gene_name": "geneC"}

    def test_interactions_and_proteins_properties(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.9)
        net.add_interaction("B", "C", 0.6)
        assert ("A", "B") == net.interactions[0][:2]
        assert net.proteins == {"A", "B", "C"}

    def test_from_interactions(self) -> None:
        net = ProteinNetwork.from_interactions([("A", "B"), ("B", "C")], name="test")
        assert net.name == "test"
        assert net.get_proteins() == ["A", "B", "C"] or set(net.get_proteins()) == {"A", "B", "C"}
        assert set(net.get_interactions()) == {("A", "B"), ("B", "C")}
        assert net.metadata["n_interactions"] == 2

    def test_from_file(self, tmp_path) -> None:
        f = _write(tmp_path, "my_net.tsv", "A\tB\t0.5\n")
        net = ProteinNetwork.from_file(f)
        assert net.name == "my_net"
        assert net.metadata["source_file"] == str(f)
        assert len(net) == 2

    def test_degree_distribution_and_neighbors(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_interaction("A", "C", 1.0)
        assert net.degree_distribution()["A"] == 2
        assert set(net.find_neighbors("A")) == {"B", "C"}
        assert net.find_neighbors("MISSING") == []

    def test_shortest_path_found_and_missing(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_interaction("B", "C", 1.0)
        assert net.shortest_path("A", "C") == ["A", "B", "C"]
        assert net.shortest_path("A", "MISSING") is None

    def test_clustering_coefficient_single_and_all(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_interaction("A", "C", 1.0)
        net.add_interaction("B", "C", 1.0)  # triangle: cc = 1.0
        assert net.clustering_coefficient("A") == 1.0
        all_cc = net.clustering_coefficient()
        assert all(v == 1.0 for v in all_cc.values())

    def test_betweenness_centrality_bridge(self) -> None:
        net = ProteinNetwork()
        for a, b in [("A", "B"), ("B", "C"), ("C", "D")]:
            net.add_interaction(a, b, 1.0)
        bc = net.betweenness_centrality()
        assert bc["B"] > 0 and bc["C"] > 0
        assert bc["B"] == max(bc.values())

    def test_components_and_largest(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_interaction("C", "D", 1.0)
        net.add_interaction("E", "F", 1.0)
        comps = net.connected_components()
        assert sorted(len(c) for c in comps) == [2, 2, 2]
        assert len(net.largest_component()) == 2

    def test_network_summary(self) -> None:
        net = ProteinNetwork(name="sumnet")
        net.add_interaction("A", "B", 1.0)
        net.metadata["source"] = "unit-test"
        summary = net.network_summary()
        assert summary["n_proteins"] == 2
        assert summary["n_interactions"] == 1
        assert summary["source"] == "unit-test"
        assert 0.0 < summary["density"] <= 1.0

    def test_protein_metadata_add_get_and_setter(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_protein_metadata("A", gene_name="geneA", function="kinase")
        md = net.get_protein_metadata("A")
        assert md["gene_name"] == "geneA" and md["function"] == "kinase"
        assert net.get_protein_metadata("MISSING") == {}
        net.protein_metadata = {"B": {"gene_name": "geneB"}}
        assert net.get_protein_metadata("B")["gene_name"] == "geneB"

    def test_get_protein_partners_with_confidence(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.9)
        net.add_interaction("A", "C", 0.3)
        assert set(net.get_protein_partners("A")) == {"B", "C"}
        assert net.get_protein_partners("A", min_confidence=0.5) == ["B"]

    def test_filter_by_confidence(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.9)
        net.add_interaction("A", "C", 0.3)
        filtered = net.filter_by_confidence(0.5)
        assert set(filtered.get_interactions()) == {("A", "B")}

    def test_filter_by_evidence(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.9, evidence_types=["experimental"])
        net.add_interaction("A", "C", 0.8, evidence_types=["textmining"])
        filtered = net.filter_by_evidence("experimental")
        assert set(filtered.get_interactions()) == {("A", "B")}

    def test_get_network_statistics(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.4)
        net.add_interaction("B", "C", 0.8)
        stats = net.get_network_statistics()
        assert stats["num_proteins"] == 3
        assert stats["num_interactions"] == 2
        assert stats["avg_confidence"] == pytest.approx(0.6)
        assert 0.0 < stats["density"] <= 1.0

    def test_to_biological_network(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.7)
        net.add_protein_metadata("A", gene_name="geneA")
        bio = net.to_biological_network()
        assert bio.graph.has_edge("A", "B")
        assert bio.graph.nodes["A"]["gene_name"] == "geneA"


class TestLoadStringInteractions:
    def test_basic_load_with_threshold(self) -> None:
        df = pd.DataFrame(
            {
                "protein1": ["A", "B"],
                "protein2": ["B", "C"],
                "combined_score": [900, 300],
            }
        )
        net = load_string_interactions(df, confidence_threshold=400)
        assert set(net.get_interactions()) == {("A", "B")}  # 300 filtered out
        assert net.name == "STRING_PPI"
        edge_attrs = net.interactions[0][2]
        assert edge_attrs["confidence"] == pytest.approx(0.9)
        assert edge_attrs["combined_score"] == 900

    def test_evidence_columns_become_evidence_types(self) -> None:
        df = pd.DataFrame(
            {
                "protein1": ["A"],
                "protein2": ["B"],
                "combined_score": [800],
                "experimental": [500],
                "database": [0],
                "textmining": [100],
            }
        )
        net = load_string_interactions(df)
        assert net.interactions[0][2]["evidence_types"] == ["experimental", "textmining"]

    def test_protein_metadata_frame(self) -> None:
        interactions = pd.DataFrame({"protein1": ["A"], "protein2": ["B"], "combined_score": [900]})
        proteins = pd.DataFrame(
            {
                "protein_id": ["A"],
                "gene_name": ["geneA"],
                "protein_name": ["Protein A"],
            }
        )
        net = load_string_interactions(interactions, proteins)
        assert net.get_protein_metadata("A")["gene_name"] == "geneA"
        assert net.get_protein_metadata("A")["protein_name"] == "Protein A"


class TestFunctionalEnrichment:
    def _annotated_network(self) -> ProteinNetwork:
        net = ProteinNetwork(name="annotated")
        edges = [
            ("K1", "K2"),
            ("K1", "K3"),
            ("K2", "K3"),  # kinase clique
            ("P1", "P2"),  # phosphatase pair
            ("K1", "P1"),  # bridge
        ]
        for a, b in edges:
            net.add_interaction(a, b, 1.0)
        for p, fn in [
            ("K1", "kinase"),
            ("K2", "kinase"),
            ("K3", "kinase"),
            ("P1", "phosphatase"),
            ("P2", "phosphatase"),
        ]:
            net.add_protein_metadata(p, function=fn)
        return net

    def test_kinase_enrichment(self) -> None:
        net = self._annotated_network()
        results = functional_enrichment_ppi(["K1", "K2", "K3"], net, min_overlap=2)
        assert "kinase" in results
        k = results["kinase"]
        assert k["count"] == 3
        assert k["enrichment_ratio"] > 1.0
        assert 0.0 <= k["p_value"] <= 1.0

    def test_min_overlap_filters(self) -> None:
        net = self._annotated_network()
        results = functional_enrichment_ppi(["K1"], net, min_overlap=2)
        assert results == {}

    def test_unannotated_only_test_set_gives_no_results(self) -> None:
        # A test set with no annotations produces no enrichment entries
        net = self._annotated_network()
        results = functional_enrichment_ppi(["UNANNOTATED1", "UNANNOTATED2"], net, min_overlap=1)
        assert results == {}

    def test_single_function_string_annotation(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        net.add_interaction("C", "D", 1.0)
        # function stored as plain string (not list)
        net.add_protein_metadata("A", function="strand")
        net.add_protein_metadata("B", function="strand")
        net.add_protein_metadata("C", function="other")
        net.add_protein_metadata("D", function="other")
        results = functional_enrichment_ppi(["A", "B"], net, min_overlap=2)
        assert "strand" in results


class TestPredictInteractions:
    def _network(self) -> ProteinNetwork:
        net = ProteinNetwork()
        # A and B share neighbors C, D; A and X share only C
        for a, b in [("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("X", "C"), ("E", "F")]:
            net.add_interaction(a, b, 1.0)
        return net

    def test_similarity_predicts_shared_neighbors(self) -> None:
        preds = predict_interactions(["A"], known_network=self._network(), method="similarity", threshold=0.3)
        assert "A" in preds
        partners = {p["partner"] for p in preds["A"]}
        assert "B" in partners  # B shares C, D with A
        assert all(p["method"] == "jaccard_similarity" for p in preds["A"])
        confs = [p["confidence"] for p in preds["A"]]
        assert confs == sorted(confs, reverse=True)

    def test_similarity_missing_protein_gives_empty(self) -> None:
        preds = predict_interactions(["NOT_PRESENT"], known_network=self._network(), method="similarity")
        assert preds["NOT_PRESENT"] == []

    def test_similarity_requires_network(self) -> None:
        assert predict_interactions(["A"], known_network=None, method="similarity") == {}

    def test_correlation_method(self) -> None:
        features = {
            "A": [1.0, 2.0, 3.0, 4.0],
            "B": [1.0, 2.0, 3.0, 3.9],
            "C": [4.0, 3.0, 2.0, 1.0],
        }
        preds = predict_interactions(["A"], features=features, method="correlation", threshold=0.5)
        partners = {p["partner"] for p in preds["A"]}
        assert "B" in partners  # highly correlated
        assert all(p["method"] == "feature_correlation" for p in preds["A"])

    def test_correlation_requires_features(self) -> None:
        assert predict_interactions(["A"], features=None, method="correlation") == {}

    def test_guilt_by_association(self) -> None:
        net = ProteinNetwork()
        # A-C, A-D, B-C, B-D, C-E: A and B share C, D; second-degree from A: B, E
        for a, b in [("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "E")]:
            net.add_interaction(a, b, 1.0)
        preds = predict_interactions(["A"], known_network=net, method="guilt-by-association", threshold=0.1)
        partners = {p["partner"] for p in preds["A"]}
        assert "B" in partners  # B shares C, D with A
        assert all(p["method"] == "guilt_by_association" for p in preds["A"])

    def test_guilt_by_association_requires_network(self) -> None:
        assert predict_interactions(["A"], known_network=None, method="guilt-by-association") == {}

    def test_method_aliases(self) -> None:
        features = {"A": [1.0, 2.0], "B": [1.0, 2.0]}
        coexpr = predict_interactions(["A"], features=features, method="coexpression", threshold=0.5)
        corr = predict_interactions(["A"], features=features, method="correlation", threshold=0.5)
        assert coexpr == corr

    def test_ml_method_falls_back_to_correlation(self) -> None:
        features = {"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0], "C": [9.0, 1.0, 5.0]}
        net = self._network()
        ml = predict_interactions(["A"], known_network=net, features=features, method="ml", threshold=0.5)
        corr = predict_interactions(["A"], features=features, method="correlation", threshold=0.5)
        assert ml == corr

    def test_ml_requires_both_inputs(self) -> None:
        assert predict_interactions(["A"], known_network=None, features=None, method="ml") == {}

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown prediction method"):
            predict_interactions(["A"], method="telepathy")

    def test_confidence_threshold_alias_kwarg(self) -> None:
        features = {"A": [1.0, 2.0], "B": [1.0, 2.0]}
        preds = predict_interactions(["A"], features=features, method="correlation", confidence_threshold=0.99)
        # threshold 0.99 still passes (correlation = 1.0)
        assert preds["A"][0]["partner"] == "B"


class TestUtilityFunctions:
    def test_protein_similarity(self) -> None:
        net = ProteinNetwork()
        for a, b in [("A", "C"), ("A", "D"), ("B", "C"), ("B", "E")]:
            net.add_interaction(a, b, 1.0)
        sim = protein_similarity(net, "A", "B")
        assert sim["common_neighbors"] == 1  # C
        assert sim["jaccard_similarity"] == pytest.approx(1 / 3)

    def test_protein_similarity_missing_protein(self) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 1.0)
        assert protein_similarity(net, "A", "MISSING") == {
            "jaccard_similarity": 0.0,
            "common_neighbors": 0,
        }

    def test_detect_complexes(self) -> None:
        net = ProteinNetwork()
        # Dense triangle (complex) + sparse chain (not a complex)
        for a, b in [("T1", "T2"), ("T2", "T3"), ("T1", "T3"), ("C1", "C2")]:
            net.add_interaction(a, b, 1.0)
        complexes = detect_complexes(net, min_size=3)
        assert complexes == [["T1", "T2", "T3"]]

    def test_detect_complexes_min_size_filters(self) -> None:
        net = ProteinNetwork()
        for a, b in [("T1", "T2"), ("T2", "T3"), ("T1", "T3")]:
            net.add_interaction(a, b, 1.0)
        assert detect_complexes(net, min_size=4) == []

    def test_export_to_string_format(self, tmp_path) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.95)
        net.add_interaction("C", "D", 0.10)
        out = tmp_path / "sub" / "string_export.tsv"
        export_to_string_format(net, out, score_threshold=500)
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "protein1\tprotein2\tcombined_score"
        assert lines[1] == "A\tB\t950"  # 0.95 * 1000; 100-score edge filtered

    def test_export_includes_all_when_threshold_zero(self, tmp_path) -> None:
        net = ProteinNetwork()
        net.add_interaction("A", "B", 0.95)
        net.add_interaction("C", "D", 0.10)
        out = tmp_path / "all.tsv"
        export_to_string_format(net, out)
        assert len(out.read_text().strip().splitlines()) == 3  # header + 2 edges
