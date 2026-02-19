"""Tests for metabolomics pathway enrichment analysis."""
from __future__ import annotations

import pytest

from metainformant.metabolomics.pathways.enrichment import (
    EnrichmentResult,
    PathwayActivityScore,
    benjamini_hochberg,
    enrichment_with_fdr,
    metabolite_set_enrichment,
    pathway_activity_scoring,
)


class TestMetaboliteSetEnrichment:
    """Tests for over-representation analysis."""

    def setup_method(self) -> None:
        self.pathway_db = {
            "glycolysis": ["glucose", "fructose", "pyruvate", "lactate", "ATP"],
            "tca_cycle": ["citrate", "succinate", "fumarate", "malate", "oxaloacetate"],
            "amino_acid": ["alanine", "glycine", "serine", "threonine", "valine"],
        }

    def test_enriched_pathway(self) -> None:
        """Query containing many glycolysis metabolites should show enrichment."""
        query = ["glucose", "fructose", "pyruvate", "lactate", "random1"]
        results = metabolite_set_enrichment(query, self.pathway_db)
        assert len(results) > 0
        # Glycolysis should be the top hit
        assert results[0].pathway_name == "glycolysis"
        assert results[0].overlap == 4
        assert results[0].fold_enrichment > 1.0

    def test_no_enrichment(self) -> None:
        """Query with no pathway overlap should have p-value = 1."""
        query = ["unknown1", "unknown2", "unknown3"]
        results = metabolite_set_enrichment(query, self.pathway_db)
        assert all(r.overlap == 0 for r in results)
        assert all(r.p_value == pytest.approx(1.0) for r in results)

    def test_sorted_by_pvalue(self) -> None:
        query = ["glucose", "citrate", "alanine"]
        results = metabolite_set_enrichment(query, self.pathway_db)
        p_vals = [r.p_value for r in results]
        assert p_vals == sorted(p_vals)

    def test_result_fields(self) -> None:
        query = ["glucose", "fructose"]
        results = metabolite_set_enrichment(query, self.pathway_db)
        for r in results:
            assert isinstance(r, EnrichmentResult)
            assert isinstance(r.pathway_name, str)
            assert r.pathway_size > 0
            assert 0 <= r.overlap <= r.pathway_size
            assert 0.0 <= r.p_value <= 1.0
            assert isinstance(r.metabolites_in_pathway, list)

    def test_custom_background_size(self) -> None:
        """Background size affects enrichment fold and p-value calculation."""
        query = ["glucose", "fructose"]
        small_bg = metabolite_set_enrichment(query, self.pathway_db, background_size=20)
        large_bg = metabolite_set_enrichment(query, self.pathway_db, background_size=2000)
        glycolysis_small = next(r for r in small_bg if r.pathway_name == "glycolysis")
        glycolysis_large = next(r for r in large_bg if r.pathway_name == "glycolysis")
        # Larger background means higher fold enrichment (more surprising to see overlap)
        assert glycolysis_large.fold_enrichment > glycolysis_small.fold_enrichment

    def test_empty_query(self) -> None:
        results = metabolite_set_enrichment([], self.pathway_db)
        assert all(r.overlap == 0 for r in results)


class TestBenjaminiHochberg:
    """Tests for FDR correction."""

    def test_monotonic_adjustment(self) -> None:
        """Adjusted values should be >= raw p-values."""
        raw = [0.001, 0.01, 0.05, 0.1, 0.5]
        adjusted = benjamini_hochberg(raw)
        for raw_p, adj_p in zip(raw, adjusted):
            assert adj_p >= raw_p

    def test_upper_bound(self) -> None:
        """Adjusted p-values should never exceed 1.0."""
        raw = [0.5, 0.8, 0.9, 0.99]
        adjusted = benjamini_hochberg(raw)
        assert all(q <= 1.0 for q in adjusted)

    def test_empty_input(self) -> None:
        assert benjamini_hochberg([]) == []

    def test_single_pvalue(self) -> None:
        assert benjamini_hochberg([0.05]) == [0.05]

    def test_preserves_order(self) -> None:
        """Output should be in the same order as input."""
        raw = [0.05, 0.001, 0.5, 0.01]
        adjusted = benjamini_hochberg(raw)
        assert len(adjusted) == 4
        # The smallest raw p-value (0.001 at index 1) should still have smallest q-value
        assert adjusted[1] == min(adjusted)


class TestEnrichmentWithFDR:
    """Tests for FDR-corrected enrichment."""

    def test_fdr_filtering(self) -> None:
        pathway_db = {
            "pathway_A": [f"met_{i}" for i in range(50)],
            "pathway_B": ["met_0", "met_1"],
            "pathway_C": [f"other_{i}" for i in range(100)],
        }
        query = [f"met_{i}" for i in range(10)]
        results = enrichment_with_fdr(query, pathway_db, fdr_threshold=0.05)
        # Only truly enriched pathways should survive FDR correction
        assert all(isinstance(r, EnrichmentResult) for r in results)

    def test_empty_query(self) -> None:
        results = enrichment_with_fdr([], {"pw": ["a", "b"]})
        assert results == []


class TestPathwayActivityScoring:
    """Tests for pathway activity scoring."""

    def setup_method(self) -> None:
        self.pathway_db = {
            "glycolysis": ["glucose", "fructose", "pyruvate"],
            "tca_cycle": ["citrate", "succinate", "malate"],
        }

    def test_positive_activity(self) -> None:
        scores = {"glucose": 2.0, "fructose": 1.5, "pyruvate": 1.0}
        results = pathway_activity_scoring(scores, self.pathway_db)
        glycolysis = next(r for r in results if r.pathway_name == "glycolysis")
        assert glycolysis.activity_score == pytest.approx(1.5)  # mean of 2.0, 1.5, 1.0
        assert glycolysis.n_measured == 3
        assert glycolysis.n_total == 3

    def test_min_coverage_filter(self) -> None:
        """Pathways below coverage threshold should be excluded."""
        scores = {"glucose": 1.0}
        results = pathway_activity_scoring(scores, self.pathway_db, min_coverage=0.5)
        # glycolysis has 1/3 coverage < 0.5
        assert not any(r.pathway_name == "glycolysis" for r in results)

    def test_sorted_by_abs_score(self) -> None:
        scores = {"glucose": -3.0, "fructose": -2.0, "pyruvate": -1.0, "citrate": 0.5, "succinate": 0.3, "malate": 0.2}
        results = pathway_activity_scoring(scores, self.pathway_db)
        abs_scores = [abs(r.activity_score) for r in results]
        assert abs_scores == sorted(abs_scores, reverse=True)

    def test_result_fields(self) -> None:
        scores = {"glucose": 1.0, "fructose": 2.0, "pyruvate": 3.0}
        results = pathway_activity_scoring(scores, self.pathway_db)
        for r in results:
            assert isinstance(r, PathwayActivityScore)
            assert isinstance(r.contributing_metabolites, dict)
            assert r.n_measured <= r.n_total

    def test_empty_scores(self) -> None:
        results = pathway_activity_scoring({}, self.pathway_db)
        assert results == []
