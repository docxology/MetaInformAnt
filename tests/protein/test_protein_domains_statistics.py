"""Tests: domain statistics count proteins, and a regular-expression hit does not carry an invented E-value.

Defects:

- ``domain_enrichment`` counted *occurrences* of a domain where its docstring says proteins (a protein with three zinc
  fingers counted three times), and grew the background by one for every domain of a foreground protein that was
  missing from the background instead of by one protein, so its hypergeometric p-values were computed on wrong counts.
- ``scan_domains`` reported an ``e_value`` obtained by pushing a heuristic score (match length plus bonuses) through the
  BLAST formula ``K * m * n * exp(-lambda * score)``. The score is not an alignment score, so the number is not an
  E-value of anything. A hit now says ``None``.

Oracle: scipy's hypergeometric distribution, ``hypergeom.sf(k - 1, N, K, n)`` (the one-sided Fisher exact p-value of
over-representation), computed from counts worked out by hand below.
"""

from __future__ import annotations

import pytest
from scipy.stats import hypergeom

from metainformant.protein.domains.classification import domain_enrichment
from metainformant.protein.domains.detection import scan_domains

BACKGROUND = [f"p{i}" for i in range(1, 21)]
ANNOTATIONS = {
    "p1": ["D1", "D1", "D1"],  # three copies of D1: one protein
    "p2": ["D1"],
    "p3": ["D1", "D2"],
    "p6": ["D1", "D1"],
    "p7": ["D1"],
    "p10": ["D2"],
}


def by_domain(results):
    return {row["domain"]: row for row in results}


def test_a_protein_counts_once_however_often_it_repeats_a_domain() -> None:
    results = by_domain(domain_enrichment(["p1", "p2", "p3", "p4", "p5"], BACKGROUND, ANNOTATIONS))

    # D1 is in p1, p2, p3 of the five foreground proteins and in p1, p2, p3, p6, p7 of the twenty proteins
    d1 = results["D1"]
    assert (d1["count_in_set"], d1["count_in_background"], d1["set_size"], d1["background_size"]) == (3, 5, 5, 20)
    assert d1["p_value"] == pytest.approx(hypergeom.sf(2, 20, 5, 5), abs=1e-6)
    assert d1["fold_enrichment"] == pytest.approx(3 / (5 / 20 * 5))  # 2.4
    # D2 is in p3 of the foreground and in p3 and p10 of the background
    d2 = results["D2"]
    assert (d2["count_in_set"], d2["count_in_background"]) == (1, 2)
    assert d2["p_value"] == pytest.approx(hypergeom.sf(0, 20, 2, 5), abs=1e-6)


def test_the_background_grows_by_one_protein_for_a_foreground_protein_that_is_missing_from_it() -> None:
    annotations = dict(ANNOTATIONS, q1=["D1", "D2", "D3"])  # q1 has three domains and is not in the background

    results = by_domain(domain_enrichment(["p1", "p2", "p3", "p4", "p5", "q1"], BACKGROUND, annotations))

    d1 = results["D1"]
    assert d1["background_size"] == 21  # twenty proteins plus q1, not plus three
    assert (d1["count_in_set"], d1["count_in_background"], d1["set_size"]) == (4, 6, 6)
    assert d1["p_value"] == pytest.approx(hypergeom.sf(3, 21, 6, 6), abs=1e-6)
    assert results["D3"]["count_in_set"] == results["D3"]["count_in_background"] == 1


def test_results_are_sorted_by_p_value_and_flag_the_significant_domains() -> None:
    background = [f"p{i}" for i in range(1, 41)]
    annotations = {f"p{i}": ["strong"] for i in range(1, 9)}  # all of the foreground and nothing else
    annotations.update({f"p{i}": ["weak"] for i in (1, 20, 21, 22, 23, 24, 25, 26)})

    results = domain_enrichment([f"p{i}" for i in range(1, 9)], background, annotations)

    assert [row["domain"] for row in results] == ["strong", "weak"]
    assert results[0]["p_value"] == pytest.approx(hypergeom.sf(7, 40, 8, 8), abs=1e-6)
    assert results[0]["significant"] is True and results[1]["significant"] is False
    assert results[0]["p_value"] <= results[1]["p_value"]


def test_a_foreground_listed_twice_is_one_set() -> None:
    once = domain_enrichment(["p1", "p2"], BACKGROUND, ANNOTATIONS)
    twice = domain_enrichment(["p1", "p2", "p1", "p2"], BACKGROUND, ANNOTATIONS)

    assert once == twice


def test_empty_sets_are_still_refused() -> None:
    with pytest.raises(ValueError):
        domain_enrichment([], BACKGROUND, ANNOTATIONS)
    with pytest.raises(ValueError):
        domain_enrichment(["p1"], [], ANNOTATIONS)


# ---------------------------------------------------------------------------
# No E-value for a pattern hit
# ---------------------------------------------------------------------------

# matches the C2H2 zinc finger pattern C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H of the built-in profiles
ZINC_FINGER = "C" + "AA" + "C" + "AAA" + "L" + "A" * 8 + "H" + "AAA" + "H"


def test_a_hit_carries_a_score_and_no_made_up_statistic() -> None:
    hits = scan_domains("M" * 10 + ZINC_FINGER + "G" * 10)

    assert [hit["name"] for hit in hits] == ["C2H2 zinc finger"]
    for hit in hits:
        assert hit["e_value"] is None  # the former value applied the BLAST formula to this heuristic score
        assert hit["score"] > 0
        assert hit["matched_sequence"] == ZINC_FINGER


def test_custom_profiles_also_get_no_e_value() -> None:
    profiles = {"X1": {"name": "tiny", "pattern": "WWW", "description": "three tryptophans"}}

    hits = scan_domains("AAWWWAA", profiles=profiles)

    assert len(hits) == 1 and hits[0]["e_value"] is None and (hits[0]["start"], hits[0]["end"]) == (2, 5)


def test_the_documentation_of_a_hit_says_what_the_score_and_the_e_value_are() -> None:
    documentation = scan_domains.__doc__ or ""

    assert "e_value: always None" in documentation
    assert "not an alignment" in documentation
