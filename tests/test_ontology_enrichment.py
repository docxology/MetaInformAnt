"""Tests for pathway enrichment: ORA, GSEA, enrichment scores, pathway networks.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import pytest

from metainformant.ontology.pathway_enrichment.enrichment import (
    compare_enrichments,
    compute_enrichment_score,
    gsea,
    over_representation_analysis,
    pathway_network,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_gene_sets():
    return {
        "pathway_A": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
        "pathway_B": ["GENE3", "GENE4", "GENE5", "GENE6", "GENE7"],
        "pathway_C": ["GENE8", "GENE9", "GENE10", "GENE11", "GENE12"],
    }


def _make_background():
    return [f"GENE{i}" for i in range(1, 51)]


# ---------------------------------------------------------------------------
# over_representation_analysis
# ---------------------------------------------------------------------------


class TestOverRepresentationAnalysis:
    def test_basic_ora(self):
        gene_list = ["GENE1", "GENE2", "GENE3"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, background=_make_background())
        assert isinstance(results, list)
        assert len(results) > 0
        assert "term_id" in results[0]
        assert "p_value" in results[0]

    def test_overlap_genes_returned(self):
        gene_list = ["GENE1", "GENE2", "GENE3"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, background=_make_background())
        pathway_a = next(r for r in results if r["term_id"] == "pathway_A")
        assert pathway_a["n_overlap"] == 3
        assert set(pathway_a["overlap_genes"]) == {"GENE1", "GENE2", "GENE3"}

    def test_no_overlap(self):
        gene_list = ["GENE50", "GENE49"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, background=_make_background())
        # All should have n_overlap = 0 or not appear
        for r in results:
            assert r["n_overlap"] == 0 or r["p_value"] >= 1.0

    def test_fdr_correction(self):
        gene_list = ["GENE1", "GENE2", "GENE3", "GENE4"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, correction="fdr_bh", background=_make_background())
        for r in results:
            assert r["adjusted_p"] >= r["p_value"] or r["adjusted_p"] == pytest.approx(r["p_value"])

    def test_bonferroni_correction(self):
        gene_list = ["GENE1", "GENE2", "GENE3"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, correction="bonferroni", background=_make_background())
        for r in results:
            assert r["adjusted_p"] <= 1.0

    def test_no_correction(self):
        gene_list = ["GENE1", "GENE2"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, correction="none", background=_make_background())
        for r in results:
            assert r["adjusted_p"] == pytest.approx(r["p_value"])

    def test_auto_background(self):
        gene_list = ["GENE1", "GENE3"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets)
        assert isinstance(results, list)

    def test_sorted_by_adjusted_p(self):
        gene_list = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, background=_make_background())
        for i in range(len(results) - 1):
            assert results[i]["adjusted_p"] <= results[i + 1]["adjusted_p"]

    def test_odds_ratio_positive(self):
        gene_list = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        gene_sets = _make_gene_sets()
        results = over_representation_analysis(gene_list, gene_sets, background=_make_background())
        pathway_a = next(r for r in results if r["term_id"] == "pathway_A")
        assert pathway_a["odds_ratio"] > 0


# ---------------------------------------------------------------------------
# compute_enrichment_score
# ---------------------------------------------------------------------------


class TestComputeEnrichmentScore:
    def test_basic_es(self):
        ranked_list = [f"GENE{i}" for i in range(1, 21)]
        gene_set = {"GENE1", "GENE2", "GENE3"}
        result = compute_enrichment_score(ranked_list, gene_set)
        assert "es" in result
        assert "running_es" in result
        assert "hit_indices" in result

    def test_no_overlap_zero_es(self):
        ranked_list = [f"GENE{i}" for i in range(1, 21)]
        gene_set = {"NOTEXIST1", "NOTEXIST2"}
        result = compute_enrichment_score(ranked_list, gene_set)
        assert result["es"] == 0.0

    def test_hit_indices_correct(self):
        ranked_list = ["A", "B", "C", "D", "E"]
        gene_set = {"B", "D"}
        result = compute_enrichment_score(ranked_list, gene_set)
        assert 1 in result["hit_indices"]
        assert 3 in result["hit_indices"]

    def test_running_es_length(self):
        ranked_list = [f"G{i}" for i in range(30)]
        gene_set = {"G5", "G10", "G20"}
        result = compute_enrichment_score(ranked_list, gene_set)
        assert len(result["running_es"]) == 30

    def test_weighted_vs_unweighted(self):
        ranked_list = [f"GENE{i}" for i in range(1, 21)]
        gene_set = {"GENE1", "GENE2", "GENE3"}
        weights = [20 - i for i in range(20)]
        weighted = compute_enrichment_score(ranked_list, gene_set, weighted=True, weights=weights)
        unweighted = compute_enrichment_score(ranked_list, gene_set, weighted=False)
        # Both should produce valid results
        assert isinstance(weighted["es"], float)
        assert isinstance(unweighted["es"], float)


# ---------------------------------------------------------------------------
# gsea
# ---------------------------------------------------------------------------


class TestGSEA:
    def test_basic_gsea(self):
        ranked_genes = [(f"GENE{i}", float(50 - i)) for i in range(1, 51)]
        gene_sets = {
            "pathway_top": [f"GENE{i}" for i in range(1, 20)],
            "pathway_bottom": [f"GENE{i}" for i in range(30, 50)],
        }
        results = gsea(ranked_genes, gene_sets, n_permutations=50, min_size=10)
        assert isinstance(results, list)
        for r in results:
            assert "es" in r
            assert "nes" in r
            assert "p_value" in r
            assert "fdr" in r

    def test_gene_set_too_small_filtered(self):
        ranked_genes = [(f"GENE{i}", float(50 - i)) for i in range(1, 51)]
        gene_sets = {"tiny": ["GENE1", "GENE2"]}
        results = gsea(ranked_genes, gene_sets, n_permutations=10, min_size=15)
        assert len(results) == 0  # Filtered out


# ---------------------------------------------------------------------------
# pathway_network
# ---------------------------------------------------------------------------


class TestPathwayNetwork:
    def test_basic_network(self):
        enrichment_results = [
            {"term_id": "pathway_A", "p_value": 0.01, "adjusted_p": 0.03},
            {"term_id": "pathway_B", "p_value": 0.02, "adjusted_p": 0.05},
            {"term_id": "pathway_C", "p_value": 0.03, "adjusted_p": 0.06},
        ]
        gene_sets = _make_gene_sets()
        result = pathway_network(enrichment_results, gene_sets, similarity_threshold=0.1)
        assert "nodes" in result
        assert "edges" in result
        assert "clusters" in result

    def test_high_threshold_no_edges(self):
        enrichment_results = [
            {"term_id": "pathway_A", "p_value": 0.01, "adjusted_p": 0.03},
            {"term_id": "pathway_C", "p_value": 0.03, "adjusted_p": 0.06},
        ]
        gene_sets = _make_gene_sets()
        result = pathway_network(enrichment_results, gene_sets, similarity_threshold=0.99)
        assert len(result["edges"]) == 0


# ---------------------------------------------------------------------------
# compare_enrichments
# ---------------------------------------------------------------------------


class TestCompareEnrichments:
    def test_basic_comparison(self):
        results_a = [
            {"term_id": "pathway_A", "adjusted_p": 0.01, "es": 2.0},
            {"term_id": "pathway_B", "adjusted_p": 0.02, "es": 1.5},
        ]
        results_b = [
            {"term_id": "pathway_A", "adjusted_p": 0.03, "es": 1.8},
            {"term_id": "pathway_C", "adjusted_p": 0.01, "es": 2.5},
        ]
        result = compare_enrichments(results_a, results_b)
        assert "shared_terms" in result
        assert "unique_a" in result
        assert "unique_b" in result
        assert "concordance" in result

    def test_shared_terms_identified(self):
        results_a = [
            {"term_id": "T1", "adjusted_p": 0.01, "es": 2.0},
            {"term_id": "T2", "adjusted_p": 0.04, "es": 1.0},
        ]
        results_b = [
            {"term_id": "T1", "adjusted_p": 0.02, "es": 1.5},
            {"term_id": "T2", "adjusted_p": 0.04, "es": 0.8},
        ]
        result = compare_enrichments(results_a, results_b)
        assert "T1" in result["shared_terms"]

    def test_concordance_direction(self):
        results_a = [
            {"term_id": "T1", "adjusted_p": 0.01, "es": 2.0},
        ]
        results_b = [
            {"term_id": "T1", "adjusted_p": 0.01, "es": 3.0},  # Same direction
        ]
        result = compare_enrichments(results_a, results_b)
        assert result["concordance"] == 1.0

    def test_empty_results(self):
        result = compare_enrichments([], [])
        assert len(result["shared_terms"]) == 0
