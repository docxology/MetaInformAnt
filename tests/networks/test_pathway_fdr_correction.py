"""Regression tests for the Benjamini-Hochberg FDR correction in pathway enrichment.

``_benjamini_hochberg`` enforced monotonicity in the wrong direction: each
adjusted p-value was clamped by the adjusted value of the next *smaller* raw
p-value instead of the next *larger* one. Every adjusted value therefore
collapsed to the smallest ``p * m / rank`` of the whole list, so a single
enriched pathway dragged every null pathway to a tiny q-value and a
``significant=True`` flag.

The oracle is ``scipy.stats.false_discovery_control(..., method="bh")``, an
independent implementation. All gene sets are real data run through the real
enrichment functions; no test doubles.
"""

from __future__ import annotations

import random

import pytest
from scipy import stats

from metainformant.networks.analysis.pathway import (
    PathwayNetwork,
    _benjamini_hochberg,
    pathway_enrichment,
    pathway_enrichment_analysis,
)


def _reference_bh(p_values: list[float]) -> list[float]:
    return [float(q) for q in stats.false_discovery_control(p_values, method="bh")]


P_VALUE_VECTORS = [
    pytest.param([0.001, 0.04, 0.9], id="small-mid-large"),
    pytest.param([0.01, 0.02, 0.03, 0.2, 0.5, 0.8], id="mixed"),
    pytest.param([0.02, 0.02, 0.9], id="tied-small-values"),
    pytest.param([0.5, 0.01, 0.2, 0.04], id="unsorted-input"),
    pytest.param([0.3], id="single-value"),
    pytest.param([1.0, 1.0, 1.0], id="all-ones"),
    pytest.param([0.0, 0.5, 1.0], id="includes-zero"),
]


@pytest.mark.parametrize("p_values", P_VALUE_VECTORS)
def test_bh_matches_the_scipy_reference(p_values: list[float]) -> None:
    assert _benjamini_hochberg(p_values) == pytest.approx(_reference_bh(p_values), rel=1e-12, abs=1e-15)


def test_bh_large_raw_p_values_are_not_pulled_down_to_the_smallest_adjusted_value() -> None:
    adjusted = _benjamini_hochberg([0.001, 0.04, 0.9])

    assert adjusted == pytest.approx([0.003, 0.06, 0.9])


def test_bh_of_empty_input_is_empty() -> None:
    assert _benjamini_hochberg([]) == []


def test_bh_is_monotone_bounded_and_matches_reference_on_random_vectors() -> None:
    rng = random.Random(20261026)
    for _ in range(60):
        p_values = [rng.random() ** rng.choice([1, 2, 4]) for _ in range(rng.randint(2, 40))]

        adjusted = _benjamini_hochberg(p_values)

        by_raw_p = [adjusted[i] for i in sorted(range(len(p_values)), key=p_values.__getitem__)]
        assert by_raw_p == sorted(by_raw_p), "adjusted values must not decrease as the raw p-value grows"
        # q >= p holds up to one unit of floating-point rounding in ``p * m / m``.
        assert all(p - 1e-12 <= q <= 1.0 for p, q in zip(p_values, adjusted))
        assert adjusted == pytest.approx(_reference_bh(p_values), rel=1e-12, abs=1e-15)


def _null_pathway_scenario() -> tuple[list[str], list[str], dict[str, list[str]]]:
    """One enriched pathway among nine pathways that share no genes with the query."""
    background = [f"g{i:02d}" for i in range(60)]
    query = background[:5]
    pathways = {"P0": background[:6]}
    for k in range(1, 10):
        start = 10 + (k - 1) * 5
        pathways[f"N{k}"] = background[start : start + 6]
    return query, background, pathways


def test_pathway_enrichment_analysis_fdr_flags_only_the_enriched_pathway() -> None:
    query, background, pathways = _null_pathway_scenario()

    results = pathway_enrichment_analysis(query, background, pathways, correction="fdr")

    assert [r["pathway"] for r in results if r["significant"]] == ["P0"]
    raw = [r["p_value"] for r in results]
    assert [r["p_value_corrected"] for r in results] == pytest.approx(_reference_bh(raw), rel=1e-12, abs=1e-15)
    assert all(r["p_value_corrected"] == pytest.approx(1.0) for r in results if r["pathway"] != "P0")


def test_pathway_enrichment_analysis_fdr_on_overlapping_pathways_matches_reference() -> None:
    pathways = {
        "glycolysis": ["G1", "G2", "G3", "G4"],
        "tca": ["G3", "G4", "G5", "G6"],
        "oxphos": ["G5", "G6", "G7", "G8"],
        "unrelated": ["U1", "U2", "U3", "U4"],
    }
    background = sorted({gene for genes in pathways.values() for gene in genes})

    results = pathway_enrichment_analysis(["G1", "G2", "G3"], background, pathways, correction="fdr")

    raw = [r["p_value"] for r in results]
    assert [r["p_value_corrected"] for r in results] == pytest.approx(_reference_bh(raw), rel=1e-12, abs=1e-15)
    by_name = {r["pathway"]: r for r in results}
    assert by_name["tca"]["p_value_corrected"] == pytest.approx(1.0)


def test_network_pathway_enrichment_fdr_matches_reference() -> None:
    background = [f"g{i:02d}" for i in range(60)]
    query = background[:8]
    pathways = {"P0": background[:6] + background[40:42]}
    for k in range(1, 6):
        start = 10 + (k - 1) * 6
        pathways[f"N{k}"] = [background[k - 1]] + background[start : start + 5]
    network = PathwayNetwork(name="fdr-regression", pathways=pathways)

    results = pathway_enrichment(query, network, background, correction="fdr")

    raw = [entry["p_value"] for entry in results.values()]
    corrected = [entry["corrected_p_value"] for entry in results.values()]
    assert corrected == pytest.approx(_reference_bh(raw), rel=1e-12, abs=1e-15)
    assert results["P0"]["corrected_p_value"] < 0.05
    assert all(results[f"N{k}"]["corrected_p_value"] > 0.5 for k in range(1, 6))
