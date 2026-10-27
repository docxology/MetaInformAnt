"""Regression tests: cross-omic concordance must split pathways at the true median CV.

``cross_omic_pathway_concordance`` calls a pathway concordant when its coefficient of variation (CV) of
-log10(p) across omics is at most the median CV. The median was ``sorted(cvs)[len(cvs) // 2]``, the upper
of the two middle values for an even number of pathways, so with four pathways three were called
concordant and only one discordant instead of two and two (and with two pathways both were concordant).

The oracle recomputes the CVs and their median with numpy, which the implementation does not use. All
p-values are real data run through the real function; no test doubles.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest

from metainformant.multiomics.pathways.enrichment import cross_omic_pathway_concordance


def _reference_split(pathway_results: dict[str, dict[str, float]]) -> tuple[list[str], list[str]]:
    """Concordant and discordant pathways: CV <= numpy's median CV, CV from the population std of -log10(p)."""
    pathways = sorted({name for p_values in pathway_results.values() for name in p_values})
    cvs: dict[str, float] = {}
    for name in pathways:
        values = np.array([-math.log10(max(p_values.get(name, 1.0), 1e-300)) for p_values in pathway_results.values()])
        cvs[name] = float(np.std(values) / np.mean(values)) if np.mean(values) > 0 else 0.0
    median_cv = float(np.median(list(cvs.values())))
    return [n for n in pathways if cvs[n] <= median_cv], [n for n in pathways if cvs[n] > median_cv]


# CVs of -log10(p) across the two omics: A 0, B 1/3, C 2/3, D 4/5 (values 2/2, 2/4, 1/5, 1/9).
FOUR_PATHWAYS = {
    "transcriptome": {"A": 1e-2, "B": 1e-2, "C": 1e-1, "D": 1e-1},
    "proteome": {"A": 1e-2, "B": 1e-4, "C": 1e-5, "D": 1e-9},
}


def test_an_even_number_of_pathways_splits_evenly_at_the_mean_of_the_middle_cvs() -> None:
    result = cross_omic_pathway_concordance(FOUR_PATHWAYS)

    # The median CV is (1/3 + 2/3) / 2 = 0.5, so C (CV 2/3) is discordant; the upper median 2/3 kept it concordant.
    assert result["concordant_pathways"] == ["A", "B"]
    assert result["discordant_pathways"] == ["C", "D"]
    assert (result["concordant_pathways"], result["discordant_pathways"]) == _reference_split(FOUR_PATHWAYS)


def test_two_pathways_are_no_longer_both_concordant() -> None:
    pathway_results = {"transcriptome": {"A": 1e-2, "B": 1e-1}, "proteome": {"A": 1e-2, "B": 1e-5}}

    result = cross_omic_pathway_concordance(pathway_results)

    assert result["concordant_pathways"] == ["A"]
    assert result["discordant_pathways"] == ["B"]


def test_an_odd_number_of_pathways_keeps_the_middle_pathway_concordant() -> None:
    pathway_results = {
        "transcriptome": {"A": 1e-2, "B": 1e-2, "C": 1e-1},
        "proteome": {"A": 1e-2, "B": 1e-4, "C": 1e-5},
    }

    result = cross_omic_pathway_concordance(pathway_results)

    assert result["concordant_pathways"] == ["A", "B"]
    assert result["discordant_pathways"] == ["C"]


def test_split_matches_numpy_on_random_p_values() -> None:
    rng = random.Random(20261026)
    even = 0
    for _ in range(40):
        n_pathways = rng.randint(2, 9)
        omics = [f"omic{i}" for i in range(rng.randint(2, 4))]
        pathway_results = {omic: {f"P{j}": 10 ** -rng.uniform(0.1, 12.0) for j in range(n_pathways)} for omic in omics}

        result = cross_omic_pathway_concordance(pathway_results)

        assert (result["concordant_pathways"], result["discordant_pathways"]) == _reference_split(pathway_results)
        even += n_pathways % 2 == 0
    assert even >= 10


@pytest.mark.parametrize("n_pathways", [2, 4, 6, 8])
def test_even_counts_with_distinct_cvs_split_half_and_half(n_pathways: int) -> None:
    # Pathway j has -log10(p) of 2 and 2 + j across the two omics, so its CV grows strictly with j.
    pathway_results = {
        "transcriptome": {f"P{j}": 1e-2 for j in range(n_pathways)},
        "proteome": {f"P{j}": 10.0 ** -(2 + j) for j in range(n_pathways)},
    }

    result = cross_omic_pathway_concordance(pathway_results)

    assert len(result["concordant_pathways"]) == len(result["discordant_pathways"]) == n_pathways // 2
