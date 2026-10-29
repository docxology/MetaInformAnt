"""``over_representation_analysis`` rejects correction names it does not implement.

The documented corrections are ``"fdr_bh"``, ``"bonferroni"`` and ``"none"``, but every other string (``"fdr"``, which
is how the networks pathway functions spell Benjamini-Hochberg, ``"BH"``, ``"holm"``, a typo) fell into the final
``else`` and left the p-values unchanged under the key ``adjusted_p``. Unknown names now raise ``ValueError``. The
oracles are the Bonferroni and Benjamini-Hochberg formulas written out in the test. No test doubles.
"""

from __future__ import annotations

import pytest

from metainformant.ontology.pathway_enrichment.enrichment import over_representation_analysis

BACKGROUND = [f"g{i}" for i in range(50)]
GENE_SETS = {"A": BACKGROUND[:8], "B": BACKGROUND[4:12], "C": BACKGROUND[30:38], "D": BACKGROUND[1:4]}
QUERY = BACKGROUND[:6]


def _by_term(rows):
    return {row["term_id"]: row for row in rows}


@pytest.mark.parametrize("name", ["fdr", "BH", "holm", "FDR_BH", "bogus", ""])
def test_unknown_corrections_raise(name: str) -> None:
    with pytest.raises(ValueError, match="Unknown correction"):
        over_representation_analysis(QUERY, GENE_SETS, background=BACKGROUND, correction=name)


def test_the_implemented_corrections_match_their_formulas() -> None:
    raw = {
        t: r["p_value"] for t, r in _by_term(over_representation_analysis(QUERY, GENE_SETS, BACKGROUND, "none")).items()
    }
    order = sorted(raw, key=raw.get)
    m = len(raw)
    bh, running = {}, 1.0
    for rank in range(m, 0, -1):
        term = order[rank - 1]
        running = min(running, raw[term] * m / rank)
        bh[term] = running

    bonferroni = {
        t: r["adjusted_p"]
        for t, r in _by_term(over_representation_analysis(QUERY, GENE_SETS, BACKGROUND, "bonferroni")).items()
    }
    fdr = {t: r["adjusted_p"] for t, r in _by_term(over_representation_analysis(QUERY, GENE_SETS, BACKGROUND)).items()}
    none = {
        t: r["adjusted_p"]
        for t, r in _by_term(over_representation_analysis(QUERY, GENE_SETS, BACKGROUND, "none")).items()
    }

    assert bonferroni == pytest.approx({t: min(p * m, 1.0) for t, p in raw.items()})
    assert fdr == pytest.approx(bh)
    assert none == pytest.approx(raw)
    assert any(abs(fdr[t] - raw[t]) > 1e-12 for t in raw)  # the correction is not a no-op on this data
