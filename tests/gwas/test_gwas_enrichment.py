"""Zero-mocks tests for gwas.analysis.enrichment.

Real computation: Fisher's exact test on hand-constructed 2x2 tables with
closed-form expected values (cross-checked against scipy.stats.fisher_exact
and scipy.stats.hypergeom.sf during development), BH FDR on known p-value
vectors, and pathway enrichment over an explicit in-memory gene universe.
"""

from __future__ import annotations

import math

from metainformant.gwas.analysis import enrichment


class TestFisherEnrichmentTest:
    def test_all_hits_in_set_enriched(self) -> None:
        # Background of 100 genes; 10 genes in the set; 5 hits all inside the set.
        background = [f"G{i:03d}" for i in range(100)]
        gene_set = background[:10]
        hits = background[:5]
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        assert result["n_hit"] == 5
        assert result["n_set"] == 10
        assert result["n_background"] == 100
        assert result["n_overlap"] == 5
        # One-sided greater test: 5/5 hits in a 10% set is significant.
        assert 0.0 < result["p_value"] < 1e-3
        assert result["odds_ratio"] > 1.0
        # Fold enrichment = (5/5)/(10/100) = 10
        assert abs(result["fold_enrichment"] - 10.0) < 1e-9
        assert sorted(result["genes_in_set"]) == sorted(hits)

    def test_no_enrichment_when_zero_overlap(self) -> None:
        background = [f"G{i:03d}" for i in range(200)]
        gene_set = background[:20]
        hits = [f"G{i:03d}" for i in range(20, 40)]  # zero overlap
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        assert result["n_overlap"] == 0
        assert result["p_value"] == 1.0
        assert result["fold_enrichment"] == 0.0

    def test_partial_overlap_enriched(self) -> None:
        # 30 hits ALL drawn from within the 20-gene set → maximally enriched
        background = [f"G{i:03d}" for i in range(100)]
        gene_set = background[:20]
        hits = background[:20]
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        assert result["n_overlap"] == 20
        assert 0.0 < result["p_value"] < 0.05
        assert result["fold_enrichment"] > 1.0

    def test_depleted_overlap_not_significant(self) -> None:
        # Hits drawn from outside the gene set: p must be large (one-sided greater).
        background = [f"G{i:03d}" for i in range(100)]
        gene_set = background[:20]
        hits = background[20:50]  # 30 hits, zero in the set
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        assert result["n_overlap"] == 0
        assert result["p_value"] == 1.0

    def test_matches_scipy_fisher_exact(self) -> None:
        pytest = __import__("pytest")
        pytest.importorskip("scipy")
        import scipy.stats as st

        background = [f"G{i:03d}" for i in range(100)]
        gene_set = background[:20]
        hits = background[:20]
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        a, b = 20, 0
        c = 0  # all 20 set genes are hits
        d = 80
        _, ref_p = st.fisher_exact([[a, b], [c, d]], alternative="greater")
        assert abs(result["p_value"] - ref_p) < 1e-12

    def test_hits_outside_background_constrained(self) -> None:
        background = [f"G{i:03d}" for i in range(100)]
        gene_set = background[:10]
        hits = ["ZZZ1", "ZZZ2"]  # not in background at all
        result = enrichment.fisher_enrichment_test(hits, gene_set, background)
        assert result["n_hit"] == 0
        assert result["n_overlap"] == 0

    def test_empty_hits(self) -> None:
        background = [f"G{i:03d}" for i in range(50)]
        result = enrichment.fisher_enrichment_test([], background[:10], background)
        assert result["n_overlap"] == 0
        assert result["p_value"] == 1.0


class TestBHCorrection:
    def test_known_bh_values(self) -> None:
        # BH step-up: sorted p=[0.005,0.01,0.03,0.04], n=4
        # raw q = p*n/rank = [0.02,0.02,0.04,0.04]; cumulative min from largest rank
        p = [0.01, 0.04, 0.03, 0.005]
        q, sig = enrichment.bh_correction(p, alpha=0.05)
        expected = {0.01: 0.02, 0.04: 0.04, 0.03: 0.04, 0.005: 0.02}
        for pi, qi in zip(p, q):
            assert abs(qi - expected[pi]) < 1e-9
        assert sig == [True, True, True, True]

    def test_bh_table_example(self) -> None:
        # Classic 10-p BH table: only the first two are significant at FDR 0.05.
        p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
        q, sig = enrichment.bh_correction(p, alpha=0.05)
        assert sig == [True, True] + [False] * 8
        # q for largest p = 0.216 (raw q = 0.216)
        assert abs(q[-1] - 0.216) < 1e-9
        assert abs(q[0] - 0.01) < 1e-9

    def test_monotonic_q_values_in_p_order(self) -> None:
        p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
        q, _sig = enrichment.bh_correction(p, alpha=0.05)
        order = sorted(range(len(p)), key=lambda i: p[i])
        qs_in_order = [q[i] for i in order]
        for a, b in zip(qs_in_order, qs_in_order[1:]):
            assert b >= a - 1e-12

    def test_all_nonsignificant(self) -> None:
        q, sig = enrichment.bh_correction([0.5, 0.9, 0.7], alpha=0.05)
        assert sig == [False, False, False]
        assert all(v <= 1.0 for v in q)

    def test_q_values_bounded_by_one(self) -> None:
        q, _ = enrichment.bh_correction([0.9, 0.95, 0.99], alpha=0.05)
        assert all(v <= 1.0 for v in q)

    def test_empty(self) -> None:
        assert enrichment.bh_correction([], alpha=0.05) == ([], [])


class TestPathwayEnrichmentFromAnnotations:
    @staticmethod
    def _hits(genes: list) -> list:
        return [{"nearby_genes": [{"gene_name": g} for g in genes]}]

    def test_explicit_gene_sets_ranked(self) -> None:
        background = [f"G{i:03d}" for i in range(100)] + ["PATH1A", "PATH1B", "PATH1C", "PATH1D"]
        hits = ["PATH1A", "PATH1B", "PATH1C", "PATH1D", "G001"]
        gene_sets = {"pathway_one": ["PATH1A", "PATH1B", "PATH1C", "PATH1D"]}
        results = enrichment.pathway_enrichment_from_annotations(self._hits(hits), background, gene_sets=gene_sets)
        assert len(results) == 1
        row = results[0]
        assert row["pathway"] == "pathway_one"
        assert row["n_overlap"] == 4
        assert 0.0 < row["p_value"] < 1e-3
        assert "q_value" in row
        assert "significant_fdr05" in row

    def test_no_hits_returns_empty_list(self) -> None:
        """Documented behavior: with zero hit genes the function short-circuits to []."""
        background = [f"G{i:03d}" for i in range(50)]
        gene_sets = {"empty_case": background[:5]}
        results = enrichment.pathway_enrichment_from_annotations([], background, gene_sets=gene_sets)
        assert results == []

    def test_default_gene_sets_available_and_sorted(self) -> None:
        background = [f"G{i:03d}" for i in range(100)] + ["HSP70", "HSP90", "SOD1"]
        hits = ["HSP70", "HSP90", "SOD1"]
        results = enrichment.pathway_enrichment_from_annotations(self._hits(hits), background)
        names = {r["pathway"] for r in results}
        assert "Stress_response" in names
        assert results == sorted(results, key=lambda r: r["p_value"])


class TestHypergeometricHelpers:
    def test_log_comb_identity(self) -> None:
        for n, k in [(10, 3), (52, 5), (100, 50)]:
            expected = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            assert abs(enrichment._log_comb(n, k) - expected) < 1e-9

    def test_sf_matches_scipy_across_grid(self) -> None:
        pytest = __import__("pytest")
        pytest.importorskip("scipy")
        import scipy.stats as st

        cases = [
            (5, 5, 10, 100),
            (10, 30, 20, 100),
            (0, 5, 10, 100),
            (3, 5, 10, 50),
            (2, 10, 8, 60),
            (7, 25, 15, 120),
        ]
        for a, n_hits, n_set, n_bg in cases:
            ours = enrichment._hypergeometric_sf(a, n_hits, n_set, n_bg)
            ref = st.hypergeom.sf(a - 1, n_bg, n_set, n_hits)
            assert abs(ours - ref) < 1e-9, (a, n_hits, n_set, n_bg, ours, ref)

    def test_sf_decreases_with_k(self) -> None:
        prev = enrichment._hypergeometric_sf(1, 10, 10, 100)
        for k in range(2, 8):
            cur = enrichment._hypergeometric_sf(k, 10, 10, 100)
            assert cur <= prev + 1e-12
            prev = cur
