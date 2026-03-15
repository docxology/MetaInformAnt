"""Tests for multiomics pathways: enrichment, active modules, topology, concordance.

Tests multi-omic enrichment (Fisher, Stouffer, min-p), active module detection,
pathway topology analysis, and cross-omic concordance scoring.

All tests use real implementations with realistic pathway/gene data (NO mocking).
"""

from __future__ import annotations

import math
from typing import Any

import pytest

from metainformant.multiomics.pathways.enrichment import (
    active_module_detection,
    cross_omic_pathway_concordance,
    multi_omic_enrichment,
    pathway_topology_analysis,
)

# ---------------------------------------------------------------------------
# Fixtures: realistic pathway and gene score data
# ---------------------------------------------------------------------------


def _make_gene_sets() -> dict[str, list[str]]:
    """Build example gene sets (pathways) with overlapping genes."""
    return {
        "KEGG_CELL_CYCLE": ["CDK1", "CDK2", "CCNB1", "CCND1", "TP53", "RB1", "E2F1"],
        "KEGG_APOPTOSIS": ["TP53", "BCL2", "BAX", "CASP3", "CASP9", "CYCS", "APAF1"],
        "KEGG_PI3K_AKT": ["PIK3CA", "AKT1", "MTOR", "PTEN", "TSC1", "TSC2"],
        "KEGG_MAPK": ["BRAF", "RAF1", "MAP2K1", "MAPK1", "MAPK3", "RAS"],
        "KEGG_WNT": ["WNT1", "CTNNB1", "APC", "GSK3B", "AXIN1"],
    }


def _make_omic_pvalues_significant() -> dict[str, dict[str, float]]:
    """Generate per-gene p-values from two omic layers, some genes significant."""
    return {
        "transcriptomics": {
            "CDK1": 0.001,
            "CDK2": 0.01,
            "CCNB1": 0.05,
            "CCND1": 0.8,
            "TP53": 0.002,
            "RB1": 0.1,
            "E2F1": 0.6,
            "BCL2": 0.003,
            "BAX": 0.02,
            "CASP3": 0.5,
            "CASP9": 0.7,
            "CYCS": 0.04,
            "APAF1": 0.3,
            "PIK3CA": 0.001,
            "AKT1": 0.005,
            "MTOR": 0.03,
            "PTEN": 0.9,
            "TSC1": 0.4,
            "TSC2": 0.5,
            "BRAF": 0.01,
            "RAF1": 0.02,
            "MAP2K1": 0.1,
            "MAPK1": 0.3,
            "MAPK3": 0.5,
            "RAS": 0.6,
            "WNT1": 0.7,
            "CTNNB1": 0.8,
            "APC": 0.9,
            "GSK3B": 0.6,
            "AXIN1": 0.5,
        },
        "proteomics": {
            "CDK1": 0.005,
            "CDK2": 0.02,
            "CCNB1": 0.1,
            "CCND1": 0.7,
            "TP53": 0.003,
            "RB1": 0.15,
            "E2F1": 0.5,
            "BCL2": 0.01,
            "BAX": 0.05,
            "CASP3": 0.6,
            "CASP9": 0.8,
            "CYCS": 0.08,
            "APAF1": 0.4,
            "PIK3CA": 0.002,
            "AKT1": 0.01,
            "MTOR": 0.05,
            "PTEN": 0.85,
            "TSC1": 0.45,
            "TSC2": 0.55,
            "BRAF": 0.015,
            "RAF1": 0.03,
            "MAP2K1": 0.12,
            "MAPK1": 0.35,
            "MAPK3": 0.55,
            "RAS": 0.65,
            "WNT1": 0.75,
            "CTNNB1": 0.85,
            "APC": 0.95,
            "GSK3B": 0.65,
            "AXIN1": 0.55,
        },
    }


def _make_pathway_graph() -> dict[str, list[str]]:
    """Build a small pathway interaction network."""
    return {
        "PIK3CA": ["AKT1"],
        "AKT1": ["MTOR", "GSK3B"],
        "PTEN": ["PIK3CA"],
        "MTOR": ["TSC1", "TSC2"],
        "TSC1": ["MTOR"],
        "TSC2": ["MTOR"],
        "GSK3B": [],
    }


def _make_network_adjacency() -> dict[str, list[str]]:
    """Build a small PPI-like network with gene names."""
    return {
        "TP53": ["MDM2", "BAX", "BCL2", "CDKN1A"],
        "MDM2": ["TP53"],
        "BAX": ["TP53", "BCL2"],
        "BCL2": ["BAX", "TP53"],
        "CDKN1A": ["TP53", "CDK2"],
        "CDK2": ["CDKN1A", "CDK1", "CCNE1"],
        "CDK1": ["CDK2", "CCNB1"],
        "CCNB1": ["CDK1"],
        "CCNE1": ["CDK2"],
        "BRCA1": ["TP53", "RAD51"],
        "RAD51": ["BRCA1"],
    }


# ===================================================================
# Multi-Omic Enrichment Tests
# ===================================================================


class TestMultiOmicEnrichment:
    """Tests for multi_omic_enrichment combining p-values across omics."""

    def test_fisher_combined_basic(self) -> None:
        """Fisher's method should return sorted enrichment results."""
        gene_sets = _make_gene_sets()
        omic_results = _make_omic_pvalues_significant()
        results = multi_omic_enrichment(gene_sets, omic_results, method="fisher_combined")

        assert len(results) == len(gene_sets)
        for r in results:
            assert "pathway_id" in r
            assert "combined_p" in r
            assert "per_omic_p" in r
            assert "n_genes" in r
            assert "leading_edge" in r
            assert 0.0 <= r["combined_p"] <= 1.0

    def test_results_sorted_by_p(self) -> None:
        """Results should be sorted by combined_p ascending."""
        gene_sets = _make_gene_sets()
        omic_results = _make_omic_pvalues_significant()
        results = multi_omic_enrichment(gene_sets, omic_results)
        pvals = [r["combined_p"] for r in results]
        assert pvals == sorted(pvals)

    def test_stouffer_method(self) -> None:
        """Stouffer's method should produce valid results."""
        gene_sets = _make_gene_sets()
        omic_results = _make_omic_pvalues_significant()
        results = multi_omic_enrichment(gene_sets, omic_results, method="stouffer")
        assert len(results) == len(gene_sets)
        for r in results:
            assert 0.0 <= r["combined_p"] <= 1.0

    def test_min_p_method(self) -> None:
        """Min-p with Bonferroni should produce valid results."""
        gene_sets = _make_gene_sets()
        omic_results = _make_omic_pvalues_significant()
        results = multi_omic_enrichment(gene_sets, omic_results, method="min_p")
        assert len(results) == len(gene_sets)
        for r in results:
            assert 0.0 <= r["combined_p"] <= 1.0

    def test_per_omic_p_present(self) -> None:
        """Each result should have per-omic p-values."""
        gene_sets = {"test_pathway": ["CDK1", "TP53"]}
        omic_results = {
            "expr": {"CDK1": 0.01, "TP53": 0.05},
            "prot": {"CDK1": 0.02, "TP53": 0.1},
        }
        results = multi_omic_enrichment(gene_sets, omic_results)
        assert "expr" in results[0]["per_omic_p"]
        assert "prot" in results[0]["per_omic_p"]

    def test_leading_edge_genes(self) -> None:
        """Leading edge should contain actual gene names from the pathway."""
        gene_sets = _make_gene_sets()
        omic_results = _make_omic_pvalues_significant()
        results = multi_omic_enrichment(gene_sets, omic_results)
        for r in results:
            pathway_genes = set(gene_sets[r["pathway_id"]])
            for gene in r["leading_edge"]:
                assert gene in pathway_genes

    def test_genes_not_in_omics(self) -> None:
        """Genes not present in omic results should not break enrichment."""
        gene_sets = {"pathway1": ["UNKNOWN_GENE1", "UNKNOWN_GENE2"]}
        omic_results = {"expr": {"OTHER_GENE": 0.5}}
        results = multi_omic_enrichment(gene_sets, omic_results)
        assert len(results) == 1
        assert results[0]["combined_p"] == 1.0  # no overlap

    def test_invalid_method_raises(self) -> None:
        """Unsupported method should raise ValueError."""
        with pytest.raises(ValueError, match="method must be"):
            multi_omic_enrichment({"p": ["G1"]}, {"o": {"G1": 0.05}}, method="bad")

    def test_empty_gene_sets_raises(self) -> None:
        """Empty gene_sets should raise ValueError."""
        with pytest.raises(ValueError, match="gene_sets must not be empty"):
            multi_omic_enrichment({}, {"o": {"G1": 0.05}})

    def test_empty_omic_results_raises(self) -> None:
        """Empty omic_results should raise ValueError."""
        with pytest.raises(ValueError, match="omic_results must not be empty"):
            multi_omic_enrichment({"p": ["G1"]}, {})

    def test_significant_pathway_has_low_p(self) -> None:
        """Pathways with consistently low p-values should have low combined_p."""
        gene_sets = {"sig_pathway": ["G1", "G2", "G3"]}
        omic_results = {
            "expr": {"G1": 0.001, "G2": 0.002, "G3": 0.005},
            "prot": {"G1": 0.001, "G2": 0.003, "G3": 0.01},
        }
        results = multi_omic_enrichment(gene_sets, omic_results)
        assert results[0]["combined_p"] < 0.05


# ===================================================================
# Active Module Detection Tests
# ===================================================================


class TestActiveModuleDetection:
    """Tests for active_module_detection subnetwork finding."""

    def test_basic_detection(self) -> None:
        """Should detect modules in a PPI network with significant node scores."""
        network = _make_network_adjacency()
        scores = {
            "TP53": 0.001,
            "MDM2": 0.01,
            "BAX": 0.005,
            "BCL2": 0.002,
            "CDKN1A": 0.03,
            "CDK2": 0.5,
            "CDK1": 0.8,
            "CCNB1": 0.9,
            "CCNE1": 0.7,
            "BRCA1": 0.6,
            "RAD51": 0.7,
        }
        modules = active_module_detection(network, scores, alpha=0.05, n_permutations=100)
        assert isinstance(modules, list)
        assert len(modules) >= 1

        for mod in modules:
            assert "module_genes" in mod
            assert "module_score" in mod
            assert "p_value" in mod
            assert len(mod["module_genes"]) >= 2

    def test_module_genes_from_network(self) -> None:
        """Module genes should be nodes in the network."""
        network = _make_network_adjacency()
        all_nodes = set(network.keys())
        for adj_list in network.values():
            all_nodes.update(adj_list)

        scores = {node: 0.001 for node in list(all_nodes)[:5]}
        scores.update({node: 0.8 for node in list(all_nodes)[5:]})

        modules = active_module_detection(network, scores, alpha=0.05, n_permutations=50)
        for mod in modules:
            for gene in mod["module_genes"]:
                assert gene in all_nodes

    def test_no_significant_nodes_returns_empty(self) -> None:
        """All non-significant nodes should yield no modules."""
        network = _make_network_adjacency()
        scores = {node: 0.9 for node in network}
        modules = active_module_detection(network, scores, alpha=0.01, n_permutations=50)
        assert modules == []

    def test_module_p_values_valid(self) -> None:
        """Module p-values should be in (0, 1]."""
        network = _make_network_adjacency()
        scores = {node: 0.001 for node in network}
        modules = active_module_detection(network, scores, alpha=0.05, n_permutations=100)
        for mod in modules:
            assert 0.0 < mod["p_value"] <= 1.0

    def test_modules_sorted_by_p(self) -> None:
        """Modules should be sorted by ascending p-value."""
        network = _make_network_adjacency()
        scores = {
            "TP53": 0.001,
            "MDM2": 0.005,
            "BAX": 0.002,
            "BCL2": 0.01,
            "CDKN1A": 0.02,
            "CDK2": 0.003,
            "CDK1": 0.004,
            "CCNB1": 0.015,
            "CCNE1": 0.025,
            "BRCA1": 0.03,
            "RAD51": 0.035,
        }
        modules = active_module_detection(network, scores, alpha=0.05, n_permutations=100)
        if len(modules) > 1:
            pvals = [m["p_value"] for m in modules]
            assert pvals == sorted(pvals)

    def test_empty_network_raises(self) -> None:
        """Empty network should raise ValueError."""
        with pytest.raises(ValueError, match="network must not be empty"):
            active_module_detection({}, {"G1": 0.01})

    def test_empty_scores_raises(self) -> None:
        """Empty scores should raise ValueError."""
        with pytest.raises(ValueError, match="scores must not be empty"):
            active_module_detection({"G1": ["G2"]}, {})

    def test_alpha_affects_seeds(self) -> None:
        """Stricter alpha should produce fewer or smaller modules."""
        network = _make_network_adjacency()
        scores = {
            "TP53": 0.001,
            "MDM2": 0.04,
            "BAX": 0.03,
            "BCL2": 0.8,
            "CDKN1A": 0.9,
            "CDK2": 0.7,
            "CDK1": 0.6,
            "CCNB1": 0.5,
            "CCNE1": 0.8,
            "BRCA1": 0.9,
            "RAD51": 0.7,
        }
        strict = active_module_detection(network, scores, alpha=0.001, n_permutations=50)
        lenient = active_module_detection(network, scores, alpha=0.05, n_permutations=50)
        # Strict should find fewer or equal modules
        strict_genes = sum(len(m["module_genes"]) for m in strict)
        lenient_genes = sum(len(m["module_genes"]) for m in lenient)
        assert strict_genes <= lenient_genes


# ===================================================================
# Pathway Topology Analysis Tests
# ===================================================================


class TestPathwayTopologyAnalysis:
    """Tests for pathway_topology_analysis."""

    def test_basic_topology(self) -> None:
        """Should compute impact factor and perturbation details."""
        graph = _make_pathway_graph()
        gene_scores = {
            "PIK3CA": 0.001,
            "AKT1": 0.005,
            "MTOR": 0.01,
            "PTEN": 0.8,
            "TSC1": 0.5,
            "TSC2": 0.6,
            "GSK3B": 0.3,
        }
        result = pathway_topology_analysis(graph, gene_scores)

        assert "impact_factor" in result
        assert "p_value" in result
        assert "perturbed_genes" in result
        assert "pathway_perturbation" in result
        assert result["impact_factor"] >= 0.0
        assert 0.0 <= result["p_value"] <= 1.0

    def test_perturbed_genes_significant(self) -> None:
        """Perturbed genes should be those with p < 0.05."""
        graph = _make_pathway_graph()
        gene_scores = {
            "PIK3CA": 0.001,
            "AKT1": 0.01,
            "MTOR": 0.03,
            "PTEN": 0.8,
            "TSC1": 0.9,
            "TSC2": 0.7,
            "GSK3B": 0.6,
        }
        result = pathway_topology_analysis(graph, gene_scores)
        for gene in result["perturbed_genes"]:
            assert gene_scores.get(gene, 1.0) < 0.05

    def test_perturbation_details_present(self) -> None:
        """Each gene should have score, centrality, weighted_score."""
        graph = _make_pathway_graph()
        gene_scores = {"PIK3CA": 0.01, "AKT1": 0.05}
        result = pathway_topology_analysis(graph, gene_scores)

        for gene, details in result["pathway_perturbation"].items():
            assert "score" in details
            assert "centrality" in details
            assert "weighted_score" in details
            assert details["score"] >= 0.0
            assert details["centrality"] > 0.0

    def test_hub_gene_higher_centrality(self) -> None:
        """Hub genes (high degree) should have higher centrality."""
        graph = {
            "HUB": ["A", "B", "C", "D"],
            "A": ["HUB"],
            "B": ["HUB"],
            "C": ["HUB"],
            "D": ["HUB"],
        }
        gene_scores = {"HUB": 0.01, "A": 0.01, "B": 0.01, "C": 0.01, "D": 0.01}
        result = pathway_topology_analysis(graph, gene_scores)
        hub_centrality = result["pathway_perturbation"]["HUB"]["centrality"]
        leaf_centrality = result["pathway_perturbation"]["A"]["centrality"]
        assert hub_centrality > leaf_centrality

    def test_empty_graph_raises(self) -> None:
        """Empty graph should raise ValueError."""
        with pytest.raises(ValueError, match="pathway_graph must not be empty"):
            pathway_topology_analysis({}, {"G1": 0.05})

    def test_empty_scores_raises(self) -> None:
        """Empty gene_scores should raise ValueError."""
        with pytest.raises(ValueError, match="gene_scores must not be empty"):
            pathway_topology_analysis({"G1": ["G2"]}, {})

    def test_all_significant_genes(self) -> None:
        """When all genes are significant, all should be in perturbed_genes."""
        graph = {"A": ["B"], "B": ["C"], "C": []}
        gene_scores = {"A": 0.001, "B": 0.01, "C": 0.02}
        result = pathway_topology_analysis(graph, gene_scores)
        assert set(result["perturbed_genes"]) == {"A", "B", "C"}


# ===================================================================
# Cross-Omic Pathway Concordance Tests
# ===================================================================


class TestCrossOmicConcordance:
    """Tests for cross_omic_pathway_concordance."""

    def test_basic_concordance(self) -> None:
        """Should classify pathways as concordant or discordant."""
        pathway_results = {
            "transcriptomics": {
                "KEGG_CELL_CYCLE": 0.001,
                "KEGG_APOPTOSIS": 0.01,
                "KEGG_WNT": 0.8,
                "KEGG_MAPK": 0.05,
            },
            "proteomics": {
                "KEGG_CELL_CYCLE": 0.002,
                "KEGG_APOPTOSIS": 0.5,  # discordant with transcriptomics
                "KEGG_WNT": 0.7,
                "KEGG_MAPK": 0.06,
            },
        }
        result = cross_omic_pathway_concordance(pathway_results)

        assert "concordant_pathways" in result
        assert "discordant_pathways" in result
        assert "concordance_score" in result
        assert "heatmap_data" in result

        all_pathways = set(result["concordant_pathways"]) | set(result["discordant_pathways"])
        assert len(all_pathways) == 4

    def test_concordance_score_bounded(self) -> None:
        """Concordance score should be in [0, 1]."""
        pathway_results = {
            "omic_a": {"P1": 0.01, "P2": 0.5},
            "omic_b": {"P1": 0.02, "P2": 0.6},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        assert 0.0 <= result["concordance_score"] <= 1.0

    def test_perfectly_concordant(self) -> None:
        """Identical p-values across omics should yield high concordance."""
        pathway_results = {
            "omic_a": {"P1": 0.001, "P2": 0.5, "P3": 0.9},
            "omic_b": {"P1": 0.001, "P2": 0.5, "P3": 0.9},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        assert result["concordance_score"] >= 0.9

    def test_discordant_detected(self) -> None:
        """Pathway significant in one omic but not another should be discordant."""
        pathway_results = {
            "omic_a": {"SIG_ONLY_A": 0.001, "BOTH_SIG": 0.001, "NEITHER": 0.9},
            "omic_b": {"SIG_ONLY_A": 0.9, "BOTH_SIG": 0.001, "NEITHER": 0.8},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        # SIG_ONLY_A has very different -log10(p) across omics
        assert "SIG_ONLY_A" in result["discordant_pathways"]

    def test_heatmap_data_structure(self) -> None:
        """Heatmap data should have -log10(p) per pathway per omic."""
        pathway_results = {
            "expr": {"P1": 0.01, "P2": 0.5},
            "prot": {"P1": 0.02, "P2": 0.6},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        for pathway_id, omic_vals in result["heatmap_data"].items():
            assert "expr" in omic_vals
            assert "prot" in omic_vals
            for val in omic_vals.values():
                assert val >= 0.0  # -log10(p) should be non-negative

    def test_pathways_only_in_one_omic(self) -> None:
        """Pathways missing from one omic should get default p=1.0."""
        pathway_results = {
            "omic_a": {"P1": 0.01, "P_ONLY_A": 0.001},
            "omic_b": {"P1": 0.02, "P_ONLY_B": 0.001},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        # All pathways should appear
        all_pathways = set(result["concordant_pathways"]) | set(result["discordant_pathways"])
        assert "P_ONLY_A" in all_pathways
        assert "P_ONLY_B" in all_pathways

    def test_fewer_than_two_omics_raises(self) -> None:
        """Fewer than 2 omic layers should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2"):
            cross_omic_pathway_concordance({"only_one": {"P1": 0.05}})

    def test_empty_raises(self) -> None:
        """Empty pathway_results should raise ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            cross_omic_pathway_concordance({})

    def test_three_omics(self) -> None:
        """Should work with three omic layers."""
        pathway_results = {
            "expr": {"P1": 0.01, "P2": 0.5},
            "prot": {"P1": 0.02, "P2": 0.6},
            "meth": {"P1": 0.015, "P2": 0.55},
        }
        result = cross_omic_pathway_concordance(pathway_results)
        assert len(result["heatmap_data"]) == 2
        for pw_data in result["heatmap_data"].values():
            assert len(pw_data) == 3
