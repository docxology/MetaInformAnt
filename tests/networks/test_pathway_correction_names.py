"""An unknown multiple-testing correction is an error; it no longer leaves p-values uncorrected.

``pathway_enrichment_analysis`` and ``pathway_enrichment`` applied Bonferroni for ``"bonferroni"``,
Benjamini-Hochberg for ``"fdr"`` and, for any other string (``"holm"``, ``"BH"``, a typo), nothing at all, then
published the unchanged p-values as ``p_value_corrected`` / ``corrected_p_value``. They now raise ``ValueError`` for
any name but ``'bonferroni'``, ``'fdr'`` and the no-correction spellings (``None``, ``'none'``). The oracles are the
Bonferroni and Benjamini-Hochberg formulas written out in the test. No test doubles.
"""

from __future__ import annotations

import pytest

pytest.importorskip("scipy")

from metainformant.networks.analysis.pathway import (  # noqa: E402
    PathwayNetwork,
    pathway_enrichment,
    pathway_enrichment_analysis,
)

BACKGROUND = [f"g{i}" for i in range(40)]
PATHWAYS = {"P1": BACKGROUND[:6], "P2": BACKGROUND[3:9], "P3": BACKGROUND[20:26], "P4": BACKGROUND[2:5]}
QUERY = BACKGROUND[:5]


def _bh(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=lambda i: p_values[i])
    adjusted = [0.0] * len(p_values)
    running = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * len(p_values) / rank)
        adjusted[index] = running
    return adjusted


@pytest.mark.parametrize("name", ["holm", "BH", "fdr_bh", "Bonferroni", "bogus"])
def test_unknown_names_raise_in_both_functions(name: str) -> None:
    with pytest.raises(ValueError, match="Unknown correction"):
        pathway_enrichment_analysis(QUERY, BACKGROUND, PATHWAYS, correction=name)
    with pytest.raises(ValueError, match="Unknown correction"):
        pathway_enrichment(QUERY, PathwayNetwork("n", PATHWAYS), background_genes=BACKGROUND, correction=name)


def test_analysis_corrections_are_the_textbook_ones() -> None:
    raw = {
        r["pathway"]: r["p_value"] for r in pathway_enrichment_analysis(QUERY, BACKGROUND, PATHWAYS, correction=None)
    }
    order = list(raw)

    bonferroni = {
        r["pathway"]: r["p_value_corrected"] for r in pathway_enrichment_analysis(QUERY, BACKGROUND, PATHWAYS)
    }
    fdr = {
        r["pathway"]: r["p_value_corrected"]
        for r in pathway_enrichment_analysis(QUERY, BACKGROUND, PATHWAYS, correction="fdr")
    }

    assert bonferroni == pytest.approx({p: min(raw[p] * len(raw), 1.0) for p in order})
    assert fdr == pytest.approx(dict(zip(order, _bh([raw[p] for p in order]))))


@pytest.mark.parametrize("spelling", [None, "none", ""])
def test_no_correction_spellings_leave_the_p_values_alone(spelling) -> None:
    rows = pathway_enrichment_analysis(QUERY, BACKGROUND, PATHWAYS, correction=spelling)
    results = pathway_enrichment(QUERY, PathwayNetwork("n", PATHWAYS), background_genes=BACKGROUND, correction=spelling)

    assert all(row["p_value_corrected"] == row["p_value"] for row in rows)
    assert all("corrected_p_value" not in result for result in results.values())


def test_the_network_function_applies_bonferroni_and_fdr() -> None:
    network = PathwayNetwork("n", PATHWAYS)
    raw = {
        pid: r["p_value"]
        for pid, r in pathway_enrichment(QUERY, network, background_genes=BACKGROUND, correction=None).items()
    }

    bonferroni = pathway_enrichment(QUERY, network, background_genes=BACKGROUND, correction="bonferroni")
    fdr = pathway_enrichment(QUERY, network, background_genes=BACKGROUND, correction="fdr")

    assert {p: r["corrected_p_value"] for p, r in bonferroni.items()} == pytest.approx(
        {p: min(v * len(raw), 1.0) for p, v in raw.items()}
    )
    assert {p: r["corrected_p_value"] for p, r in fdr.items()} == pytest.approx(dict(zip(raw, _bh(list(raw.values())))))
