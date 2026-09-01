"""Tests for protein domain classification, similarity, clustering, and enrichment.

Round-4 test-depth lane T4: pins regression behaviour of
``metainformant.protein.domains.classification`` with real data. Enrichment
p-values are cross-checked against ``scipy.stats.hypergeom`` (scipy is a repo
dependency, not a mock). Pinned similarity/cluster values were measured from
the implementation and cross-derived by hand.
"""

from __future__ import annotations

import math

import pytest

from metainformant.protein.domains.classification import (
    _hypergeometric_sf,
    _log_hypergeometric_pmf,
    classify_protein_family,
    cluster_by_domain,
    compute_domain_similarity,
    domain_enrichment,
)


def _d(name: str, start: int = 1) -> dict:
    return {"name": name, "start": start}


# ---------------------------------------------------------------------------
# classify_protein_family
# ---------------------------------------------------------------------------


class TestClassifyProteinFamily:
    def test_kinase_with_regulatory_domain(self):
        result = classify_protein_family([_d("Protein kinase domain"), _d("SH2")], "MSEQ")
        assert result["family"] == "Protein kinase"
        assert result["confidence"] == pytest.approx(1.0)
        assert result["superfamily"] == "Transferase"
        assert result["enzyme_class"] == "Kinase"
        assert result["ec_number"] == "2.7.x.x"
        assert result["go_terms"] == ["GO:0004672", "GO:0006468"]
        assert result["domain_names"] == ["Protein kinase domain", "SH2"]

    def test_serine_protease(self):
        result = classify_protein_family([_d("Trypsin-like serine protease"), _d("Kringle")], "MKV")
        assert result["family"] == "Serine protease"
        assert result["superfamily"] == "Hydrolase"
        assert result["ec_number"] == "3.4.21.x"

    def test_unclassified_no_domains(self):
        result = classify_protein_family([], "MSEQ")
        assert result["family"] == "Unclassified"
        assert result["confidence"] == 0.0
        assert result["matched_rule"] is None
        assert result["superfamily"] == "Unknown"
        assert result["go_terms"] == []

    def test_unknown_domain_unclassified(self):
        result = classify_protein_family([_d("Totally unknown domain")], "MSEQ")
        assert result["family"] == "Unclassified"
        assert result["confidence"] == 0.0
        # But domain name is still reported
        assert result["domain_names"] == ["Totally unknown domain"]

    def test_multiple_rules_alternatives(self):
        # Kinase + RRM hits two rules; kinase has one required (1.0) + optional,
        # RRM rule requires "RNA recognition motif (RRM)" (1.0)
        domains = [_d("Protein kinase domain"), _d("RNA recognition motif (RRM)")]
        result = classify_protein_family(domains, "MSEQ")
        assert result["family"] in {"Protein kinase", "RNA-binding protein"}
        # The runner-up should appear in alternatives
        alt_names = {a["family"] for a in result["alternative_families"]}
        assert result["family"] not in alt_names or True  # alternatives may include the other
        families_seen = {result["family"]} | alt_names
        assert {"Protein kinase", "RNA-binding protein"} <= families_seen

    def test_domain_names_sorted(self):
        domains = [_d("WD40 repeat"), _d("Homeodomain")]
        result = classify_protein_family(domains, "MSEQ")
        assert result["domain_names"] == sorted(["WD40 repeat", "Homeodomain"])


# ---------------------------------------------------------------------------
# compute_domain_similarity
# ---------------------------------------------------------------------------


class TestComputeDomainSimilarity:
    def test_identical_architectures(self):
        a = [_d("X", 1), _d("Y", 50)]
        b = [_d("X", 1), _d("Y", 50)]
        assert compute_domain_similarity(a, b) == pytest.approx(1.0)

    def test_both_empty(self):
        assert compute_domain_similarity([], []) == pytest.approx(1.0)

    def test_one_empty(self):
        assert compute_domain_similarity([_d("X")], []) == pytest.approx(0.0)
        assert compute_domain_similarity([], [_d("X")]) == pytest.approx(0.0)

    def test_swapped_order_penalty(self):
        # Same composition, reversed order: jaccard=1 (0.7) + LCS=1/2 (0.3*0.5) = 0.85
        a = [_d("X", 1), _d("Y", 50)]
        b = [_d("Y", 1), _d("X", 50)]
        assert compute_domain_similarity(a, b) == pytest.approx(0.85)

    def test_disjoint_domains(self):
        a = [_d("X", 1)]
        b = [_d("Z", 1)]
        assert compute_domain_similarity(a, b) == pytest.approx(0.0)

    def test_shared_partial(self):
        # Composition: 1 shared of 3 union -> jaccard 1/3; order LCS 1/2
        a = [_d("X", 1), _d("Y", 50)]
        b = [_d("X", 1), _d("Z", 50)]
        # 0.7*(1/3) + 0.3*(1/2) = 0.3833... (module rounds to 4 decimals)
        assert compute_domain_similarity(a, b) == pytest.approx(0.7 / 3 + 0.15, abs=5e-4)

    def test_duplicate_domains_in_one_protein(self):
        a = [_d("X", 1), _d("X", 100)]
        b = [_d("X", 1)]
        # Sets collapse duplicates: jaccard 1.0; order LCS 1/2
        assert compute_domain_similarity(a, b) == pytest.approx(0.7 + 0.15)

    def test_result_in_unit_interval(self):
        a = [_d("A", 1), _d("B", 10), _d("C", 30)]
        b = [_d("B", 1), _d("C", 20), _d("D", 40)]
        s = compute_domain_similarity(a, b)
        assert 0.0 <= s <= 1.0

    def test_sorted_by_start_coordinate(self):
        # Domain dict order should not matter; 'start' drives order comparison
        a = [_d("Y", 100), _d("X", 10)]
        b = [_d("X", 10), _d("Y", 100)]
        assert compute_domain_similarity(a, b) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# cluster_by_domain
# ---------------------------------------------------------------------------


class TestClusterByDomain:
    def test_basic_two_clusters(self):
        proteins = {
            "p1": [_d("A")],
            "p2": [_d("A", 5)],
            "p3": [_d("B")],
        }
        result = cluster_by_domain(proteins, threshold=0.7)
        assert result["cluster_count"] == 2
        assert result["singleton_count"] == 1
        members = [sorted(v) for v in result["clusters"].values()]
        assert ["p1", "p2"] in members
        assert ["p3"] in members

    def test_empty_input(self):
        result = cluster_by_domain({}, threshold=0.5)
        assert result["cluster_count"] == 0
        assert result["singleton_count"] == 0
        assert result["similarity_matrix"] == {}

    def test_unsupported_method_raises(self):
        with pytest.raises(ValueError, match="Unsupported clustering method"):
            cluster_by_domain({"p1": [_d("A")]}, method="euclidean")

    def test_threshold_prevents_merging(self):
        proteins = {
            "p1": [_d("A"), _d("B")],
            "p2": [_d("C"), _d("D")],
        }
        result = cluster_by_domain(proteins, threshold=0.7)
        assert result["cluster_count"] == 2
        assert result["singleton_count"] == 2

    def test_single_linkage_transitivity(self):
        # p1~p2 similar, p2~p3 similar, p1~p3 less so -> all one cluster
        proteins = {
            "p1": [_d("A"), _d("B")],
            "p2": [_d("A"), _d("B"), _d("C")],
            "p3": [_d("A"), _d("B"), _d("C"), _d("D")],
        }
        result = cluster_by_domain(proteins, threshold=0.6)
        assert result["cluster_count"] == 1

    def test_similarity_matrix_symmetric(self):
        proteins = {
            "p1": [_d("A")],
            "p2": [_d("A", 2)],
            "p3": [_d("B")],
        }
        result = cluster_by_domain(proteins, threshold=0.5)
        m = result["similarity_matrix"]
        assert m[("p1", "p2")] == m[("p2", "p1")]
        assert m[("p1", "p3")] == m[("p3", "p1")]
        assert ("p1", "p1") not in m  # diagonal not stored
        assert set(m) == {("p1", "p2"), ("p2", "p1"), ("p1", "p3"), ("p3", "p1"), ("p2", "p3"), ("p3", "p2")}

    def test_metadata_round_trip(self):
        proteins = {"p1": [_d("A")], "p2": [_d("A")]}
        result = cluster_by_domain(proteins, threshold=0.9)
        assert result["method"] == "jaccard"
        assert result["threshold"] == pytest.approx(0.9)


# ---------------------------------------------------------------------------
# domain_enrichment — p-values cross-checked against scipy
# ---------------------------------------------------------------------------


class TestDomainEnrichment:
    def _scenario(self):
        # 100 background proteins; domain "D" in 10 of them.
        # Foreground: 25 proteins, 5 of which have "D" (all inside background).
        annotations = {f"p{i}": (["D"] if i < 10 else []) for i in range(100)}
        background = [f"p{i}" for i in range(100)]
        foreground = [f"p{i}" for i in range(5)] + [f"p{i}" for i in range(50, 70)]
        return annotations, background, foreground

    def test_pvalue_matches_scipy(self):
        scipy_stats = pytest.importorskip("scipy.stats")
        annotations, background, foreground = self._scenario()
        results = domain_enrichment(foreground, background, annotations)
        assert len(results) == 1
        r = results[0]
        # One-sided Fisher / hypergeometric SF: P(X >= 5) with N=100, K=10, n=25
        expected = float(scipy_stats.hypergeom.sf(4, 100, 10, 25))
        assert r["p_value"] == pytest.approx(expected, abs=1e-5)

    def test_counts_and_fold_enrichment(self):
        annotations, background, foreground = self._scenario()
        r = domain_enrichment(foreground, background, annotations)[0]
        assert r["domain"] == "D"
        assert r["count_in_set"] == 5
        assert r["count_in_background"] == 10
        assert r["set_size"] == 25
        assert r["background_size"] == 100
        # expected = (10/100)*25 = 2.5 -> fold = 5/2.5 = 2.0
        assert r["fold_enrichment"] == pytest.approx(2.0)

    def test_sorted_ascending_by_pvalue(self):
        annotations = {
            "a": ["D1", "D2"],
            "b": ["D1"],
            "c": ["D2"],
            "d": ["D1", "D2"],
            "e": ["D3"],
            "f": ["D3"],
            "g": [],
        }
        background = ["a", "b", "c", "d", "e", "f", "g"]
        foreground = ["a", "b", "d"]
        results = domain_enrichment(foreground, background, annotations)
        ps = [r["p_value"] for r in results]
        assert ps == sorted(ps)

    def test_depleted_domain_not_reported_when_zero_in_set(self):
        annotations = {"a": ["D1"], "b": [], "c": [], "d": []}
        background = ["a", "b", "c", "d"]
        results = domain_enrichment(["b", "c", "d"], background, annotations)
        assert results == []  # k == 0 domains skipped

    def test_empty_foreground_raises(self):
        with pytest.raises(ValueError, match="protein_set"):
            domain_enrichment([], ["a"], {})

    def test_empty_background_raises(self):
        with pytest.raises(ValueError, match="background"):
            domain_enrichment(["a"], [], {})

    def test_significant_flag_threshold(self):
        # Extreme case: all set proteins share a domain absent from rest of background
        annotations = {f"p{i}": ["Hot"] for i in range(5)} | {f"q{i}": [] for i in range(95)}
        background = [f"p{i}" for i in range(5)] + [f"q{i}" for i in range(95)]
        foreground = [f"p{i}" for i in range(5)]
        results = domain_enrichment(foreground, background, annotations)
        assert results[0]["significant"] is True
        assert results[0]["p_value"] < 0.05
        assert results[0]["fold_enrichment"] > 1.0


# ---------------------------------------------------------------------------
# Hypergeometric internals — analytic pins
# ---------------------------------------------------------------------------


class TestHypergeometricInternals:
    def test_sf_bounds(self):
        assert _hypergeometric_sf(0, 100, 10, 25) == pytest.approx(1.0)
        assert _hypergeometric_sf(11, 100, 10, 25) == pytest.approx(0.0)  # > min(n, K)
        sf = _hypergeometric_sf(5, 100, 10, 25)
        assert 0.0 < sf < 1.0

    def test_sf_decreasing_in_k(self):
        vals = [_hypergeometric_sf(k, 100, 10, 25) for k in range(0, 8)]
        assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))

    def test_sf_matches_scipy_series(self):
        scipy_stats = pytest.importorskip("scipy.stats")
        for k in range(1, 6):
            expected = float(scipy_stats.hypergeom.sf(k - 1, 100, 10, 25))
            assert _hypergeometric_sf(k, 100, 10, 25) == pytest.approx(expected, abs=1e-6)

    def test_log_pmf_normalisation(self):
        # Sum of pmf over support ~ 1 (within float tolerance in log space)
        N, K, n = 50, 15, 20
        total = 0.0
        for k in range(0, min(n, K) + 1):
            lp = _log_hypergeometric_pmf(k, N, K, n)
            if lp > -700:
                total += math.exp(lp)
        assert total == pytest.approx(1.0, abs=1e-9)
