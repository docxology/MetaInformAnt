"""Regression tests for the hypergeometric enrichment branch and the corrected p-value of a lone pathway.

Three defects are covered:

* ``pathway_enrichment_analysis(method="hypergeometric")`` computed its odds ratio as
  ``(a / max(b, 1)) / (c / max(d, 1))``. ``c`` is the number of background genes outside the query that
  belong to the pathway, so it is zero whenever every pathway gene is in the query (or the pathway has no
  gene in the background at all), and the division raised ``ZeroDivisionError``. The ratio is now the
  sample odds ratio ``(a * d) / (b * c)``, which is ``inf`` there, exactly as ``scipy.stats.fisher_exact``
  reports it.
* The same branch took the success count ``n`` from the whole pathway instead of from the pathway genes in
  the background (the urn the genes are drawn from). A pathway listing genes the background does not
  contain therefore got a larger ``n`` than the population can hold, which understated the enrichment.
  Pathways, and the query, are now restricted to the background, as ``pathway_enrichment`` already did.
* ``pathway_enrichment`` only attached ``corrected_p_value`` when more than one pathway passed
  ``min_overlap``, so ``results[pid]["corrected_p_value"]`` raised ``KeyError`` for a lone survivor.

The oracle is ``scipy.stats.hypergeom.sf`` (upper tail) and ``scipy.stats.fisher_exact``, independent of
the implementation, fed with set sizes computed here. All gene sets are real data run through the real
enrichment functions; no test doubles.
"""

from __future__ import annotations

import math
import random
from typing import Any

import pytest
from scipy import stats

from metainformant.networks.analysis.pathway import PathwayNetwork, pathway_enrichment, pathway_enrichment_analysis

BACKGROUND = [f"g{i:02d}" for i in range(40)]


def _reference(query: list[str], background: list[str], pathway: list[str]) -> dict[str, Any]:
    """Reference values for one pathway, from set sizes and scipy alone."""
    universe = set(background)
    drawn = set(query) & universe
    members = set(pathway) & universe
    m, n, big_n, k = len(universe), len(members), len(drawn), len(drawn & members)
    table = [[k, big_n - k], [n - k, m - n - (big_n - k)]]
    return {
        "p_value": float(stats.hypergeom.sf(k - 1, m, n, big_n)),
        "odds_ratio": float(stats.fisher_exact(table, alternative="greater").statistic),
        "genes_in_pathway": k,
        "pathway_size": n,
        "genes_tested": big_n,
        "expected": n / m * big_n,
    }


def _by_name(results: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {result["pathway"]: result for result in results}


ALL_PATHWAY_GENES_IN_QUERY = [
    pytest.param(BACKGROUND[:10], BACKGROUND[:4], id="pathway-inside-a-larger-query"),
    pytest.param(BACKGROUND[:4], BACKGROUND[:4], id="query-equals-pathway"),
    pytest.param(BACKGROUND, BACKGROUND[:6], id="query-is-the-whole-background"),
    pytest.param(BACKGROUND[:10], BACKGROUND[:4] + ["x1", "x2"], id="rest-of-pathway-is-outside-the-background"),
]


@pytest.mark.parametrize(("query", "pathway"), ALL_PATHWAY_GENES_IN_QUERY)
def test_pathway_fully_inside_the_query_no_longer_divides_by_zero(query: list[str], pathway: list[str]) -> None:
    results = pathway_enrichment_analysis(query, BACKGROUND, {"covered": pathway}, method="hypergeometric")

    reference = _reference(query, BACKGROUND, pathway)
    (result,) = results
    assert result["p_value"] == pytest.approx(reference["p_value"], rel=1e-12)
    assert result["odds_ratio"] == pytest.approx(reference["odds_ratio"], nan_ok=True)
    assert result["genes_in_pathway"] == reference["genes_in_pathway"] == reference["pathway_size"]


def test_odds_ratio_is_infinite_when_no_background_gene_outside_the_query_is_in_the_pathway() -> None:
    (result,) = pathway_enrichment_analysis(
        BACKGROUND[:10], BACKGROUND, {"covered": BACKGROUND[:4]}, method="hypergeometric"
    )

    assert math.isinf(result["odds_ratio"])
    assert result["odds_ratio"] > 0


def test_pathway_genes_outside_the_background_are_not_counted() -> None:
    query = BACKGROUND[:3]
    pathways = {
        "padded": BACKGROUND[:5] + ["x1", "x2", "x3"],  # 5 of its 8 genes are in the background
        "elsewhere": ["x4", "x5", "x6", "x7"],  # none are
    }

    results = _by_name(pathway_enrichment_analysis(query, BACKGROUND, pathways, method="hypergeometric"))

    padded = _reference(query, BACKGROUND, pathways["padded"])
    assert padded["pathway_size"] == 5
    assert results["padded"]["pathway_size"] == 5
    assert results["padded"]["p_value"] == pytest.approx(padded["p_value"], rel=1e-12)
    assert results["padded"]["enrichment_ratio"] == pytest.approx(3 / padded["expected"])
    assert results["elsewhere"]["pathway_size"] == 0
    assert results["elsewhere"]["genes_in_pathway"] == 0
    assert results["elsewhere"]["p_value"] == pytest.approx(1.0)


@pytest.mark.parametrize("method", ["fisher", "hypergeometric"])
def test_query_genes_outside_the_background_are_ignored(method: str) -> None:
    pathways = {"P": BACKGROUND[:3] + BACKGROUND[10:13], "Q": BACKGROUND[20:26]}

    with_noise = pathway_enrichment_analysis(BACKGROUND[:5] + ["zz1", "zz2"], BACKGROUND, pathways, method=method)
    without = pathway_enrichment_analysis(BACKGROUND[:5], BACKGROUND, pathways, method=method)

    assert with_noise == without
    assert {result["genes_tested"] for result in with_noise} == {5}


def test_hypergeometric_results_match_scipy_on_random_gene_sets() -> None:
    rng = random.Random(20261026)
    background = [f"g{i:03d}" for i in range(80)]
    covered = 0
    for _ in range(40):
        query = rng.sample(background, rng.randint(1, 30))
        outside = [f"x{i}" for i in range(rng.randint(0, 4))]
        pathways = {}
        for index in range(rng.randint(1, 5)):
            if rng.random() < 0.3:  # every background gene of this pathway is in the query
                members = rng.sample(query, rng.randint(1, min(len(query), 8)))
            else:
                members = rng.sample(background, rng.randint(1, 15))
            pathways[f"P{index}"] = members + outside

        results = _by_name(pathway_enrichment_analysis(query, background, pathways, method="hypergeometric"))

        assert set(results) == set(pathways)
        for name, genes in pathways.items():
            reference = _reference(query, background, genes)
            result = results[name]
            covered += reference["genes_in_pathway"] == reference["pathway_size"]
            assert result["p_value"] == pytest.approx(reference["p_value"], rel=1e-12, abs=1e-300), name
            assert result["odds_ratio"] == pytest.approx(reference["odds_ratio"], nan_ok=True), name
            assert result["genes_in_pathway"] == reference["genes_in_pathway"], name
            assert result["pathway_size"] == reference["pathway_size"], name
            assert result["genes_tested"] == reference["genes_tested"], name
            assert result["enrichment_ratio"] == pytest.approx(
                reference["genes_in_pathway"] / reference["expected"], rel=1e-12
            ), name
    assert covered > 5, "the scenarios must include pathways fully inside the query"


def _network(pathways: dict[str, list[str]]) -> PathwayNetwork:
    return PathwayNetwork(name="lone-survivor", pathways=pathways)


@pytest.mark.parametrize("correction", ["bonferroni", "fdr"])
@pytest.mark.parametrize("method", ["fisher", "hypergeometric"])
def test_lone_pathway_above_min_overlap_still_gets_a_corrected_p_value(method: str, correction: str) -> None:
    network = _network({"hit": BACKGROUND[:6], "missed": BACKGROUND[20:26]})
    query = BACKGROUND[:4]

    results = pathway_enrichment(query, network, BACKGROUND, method=method, correction=correction, min_overlap=1)

    assert list(results) == ["hit"]
    assert results["hit"]["p_value"] == pytest.approx(_reference(query, BACKGROUND, BACKGROUND[:6])["p_value"])
    assert results["hit"]["corrected_p_value"] == pytest.approx(results["hit"]["p_value"])


@pytest.mark.parametrize("correction", ["bonferroni", "fdr"])
def test_min_overlap_can_leave_exactly_one_pathway(correction: str) -> None:
    network = _network({"strong": BACKGROUND[:6], "weak": BACKGROUND[3:5] + BACKGROUND[20:24]})
    query = BACKGROUND[:4]  # overlap with the query: strong 4 genes, weak 1

    results = pathway_enrichment(query, network, BACKGROUND, correction=correction, min_overlap=2)

    assert list(results) == ["strong"]
    assert results["strong"]["corrected_p_value"] == pytest.approx(results["strong"]["p_value"])


def test_no_pathway_above_min_overlap_returns_no_results() -> None:
    network = _network({"far": BACKGROUND[20:26]})

    assert pathway_enrichment(BACKGROUND[:4], network, BACKGROUND, min_overlap=1) == {}


def test_correction_of_several_pathways_is_unchanged() -> None:
    network = _network({"a": BACKGROUND[:6], "b": BACKGROUND[2:8], "c": BACKGROUND[20:26]})

    results = pathway_enrichment(BACKGROUND[:4], network, BACKGROUND, correction="bonferroni", min_overlap=1)

    assert set(results) == {"a", "b"}
    for entry in results.values():
        assert entry["corrected_p_value"] == pytest.approx(min(entry["p_value"] * 2, 1.0))
