"""Regression tests for ``enrich_genes``, its no-scipy p-value fallback and the ``semantic_similarity`` documentation.

* ``enrich_genes`` took an ``alpha`` ("significance threshold") that nothing used: it only appeared in a log message,
  and the returned frame carried no significance information. ``alpha`` now sets a boolean ``significant`` column
  (``bonferroni_p <= alpha``); every tested term is still returned.
* Without scipy, ``_hypergeometric_test`` returned the point probability ``P(X = k)`` computed from a handful of
  ``lgamma`` terms ("very simplified approximation"), not the upper tail ``P(X >= k)`` that the scipy branch returns,
  and fell back to 1.0 on any numerical error. ``_fisher_exact_test`` shares it. Both now compute the exact tail.
* ``semantic_similarity`` documented "a score between 0 and 1" for every method, but Resnik returns the information
  content of the most informative common ancestor (bits), which is not bounded by 1. The documentation now says so.

The oracles are ``scipy.stats.fisher_exact`` / ``scipy.stats.hypergeom`` fed with sets built in the test. The no-scipy
cases run the real module in a subprocess in which scipy cannot be imported. No test doubles.
"""

from __future__ import annotations

import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest
from scipy import stats

import metainformant
from metainformant.ontology.core import go as go_module
from metainformant.ontology.core.go import enrich_genes, semantic_similarity

BACKGROUND = [f"G{i:02d}" for i in range(60)]
QUERY = BACKGROUND[:10]
ANNOTATIONS = {
    "T_strong": set(BACKGROUND[:8]) | set(BACKGROUND[30:34]),
    "T_medium": set(BACKGROUND[:5]) | set(BACKGROUND[20:25]),
    "T_weak": {BACKGROUND[0]} | set(BACKGROUND[40:60]),
    "T_none": set(BACKGROUND[50:60]),  # shares no gene with the query: not tested
}


def _fisher_oracle(query: list[str], background: list[str], term_genes: set[str]) -> float:
    universe = set(background)
    drawn = set(query) & universe  # genes outside the background are not tested
    members = term_genes & universe
    k = len(drawn & members)
    table = [[k, len(drawn) - k], [len(members) - k, len(universe) - len(drawn) - (len(members) - k)]]
    return float(stats.fisher_exact(table, alternative="greater")[1])


def _oracle_frame() -> dict[str, tuple[float, float]]:
    """term -> (p, bonferroni p) for the terms that share a gene with the query."""
    p_values = {
        term: _fisher_oracle(QUERY, BACKGROUND, genes) for term, genes in ANNOTATIONS.items() if set(QUERY) & genes
    }
    return {term: (p, min(p * len(p_values), 1.0)) for term, p in p_values.items()}


# ---------------------------------------------------------------------------
# alpha
# ---------------------------------------------------------------------------


def test_alpha_sets_the_significant_flag_against_the_bonferroni_p_values() -> None:
    oracle = _oracle_frame()
    ordered = sorted(oracle.values(), key=lambda pair: pair[1])
    strict_alpha = (ordered[0][1] + ordered[1][1]) / 2  # only the strongest term passes
    loose_alpha = (ordered[1][1] + ordered[2][1]) / 2  # the two enriched terms pass

    strict = enrich_genes(QUERY, BACKGROUND, ANNOTATIONS, alpha=strict_alpha)
    loose = enrich_genes(QUERY, BACKGROUND, ANNOTATIONS, alpha=loose_alpha)

    for frame, alpha in ((strict, strict_alpha), (loose, loose_alpha)):
        assert set(frame["go_term"]) == set(oracle), "every tested term stays in the result"
        assert frame["p_value"].is_monotonic_increasing
        for _, row in frame.iterrows():
            expected_p, expected_bonferroni = oracle[row["go_term"]]
            assert row["p_value"] == pytest.approx(expected_p, rel=1e-9)
            assert row["bonferroni_p"] == pytest.approx(expected_bonferroni, rel=1e-9)
            assert bool(row["significant"]) == (expected_bonferroni <= alpha)
    assert strict["significant"].sum() == 1
    assert loose["significant"].sum() == 2


def test_default_alpha_is_the_conventional_0_05() -> None:
    frame = enrich_genes(QUERY, BACKGROUND, ANNOTATIONS)

    assert dict(zip(frame["go_term"], frame["significant"])) == {"T_strong": True, "T_medium": True, "T_weak": False}


def test_a_very_small_alpha_flags_nothing_but_keeps_every_term() -> None:
    frame = enrich_genes(QUERY, BACKGROUND, ANNOTATIONS, alpha=1e-12)

    assert len(frame) == 3
    assert not frame["significant"].any()


# ---------------------------------------------------------------------------
# p-values without scipy
# ---------------------------------------------------------------------------

NO_SCIPY_SCRIPT = """
import json, sys

sys.modules["scipy"] = None  # any `import scipy...` now raises ImportError, as on a machine without scipy

from metainformant.ontology.core.go import enrich_genes

payload = json.load(sys.stdin)
annotations = {term: set(genes) for term, genes in payload["annotations"].items()}
out = {}
for method in ("fisher", "hypergeometric"):
    frame = enrich_genes(payload["query"], None, annotations, method=method)
    out[method] = dict(zip(frame["go_term"], frame["p_value"]))
print(json.dumps(out))
"""


def _subprocess_env() -> dict[str, str]:
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    return {**os.environ, "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else []))}


def test_enrichment_without_scipy_returns_the_exact_upper_tail() -> None:
    # background=None means "all annotated genes", so the oracle universe is the union of the annotation sets
    # (G08 and G09 are in no annotation set, so they drop out of the query)
    universe = sorted(set().union(*ANNOTATIONS.values()))
    payload = {"query": QUERY, "annotations": {term: sorted(genes) for term, genes in ANNOTATIONS.items()}}

    completed = subprocess.run(
        [sys.executable, "-c", NO_SCIPY_SCRIPT],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=_subprocess_env(),
        timeout=300,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    observed = json.loads(completed.stdout.strip().splitlines()[-1])

    expected = {
        term: _fisher_oracle(QUERY, universe, genes) for term, genes in ANNOTATIONS.items() if set(QUERY) & genes
    }
    assert expected["T_strong"] < 1e-3 and 0.9 < expected["T_weak"], "scenario must span enriched and depleted terms"
    for method in ("fisher", "hypergeometric"):
        assert observed[method] == pytest.approx(expected, rel=1e-9), method


def test_upper_tail_helper_matches_scipy_on_random_tables() -> None:
    rng = random.Random(20261027)
    for _ in range(400):
        population = rng.randint(1, 300)
        successes = rng.randint(0, population)
        draws = rng.randint(0, population)
        k = rng.randint(-1, min(successes, draws) + 2)

        reference = float(stats.hypergeom.sf(k - 1, population, successes, draws))

        assert go_module._hypergeometric_upper_tail(k, population, successes, draws) == pytest.approx(
            reference, rel=1e-12, abs=1e-300
        ), (k, population, successes, draws)


def test_upper_tail_is_not_the_point_probability() -> None:
    # P(X >= 3) = 1/30 for 3 draws from 10 with 4 successes... and P(X = 2) is larger than P(X >= 3)
    assert go_module._hypergeometric_upper_tail(2, 10, 4, 3) == pytest.approx(float(stats.hypergeom.sf(1, 10, 4, 3)))
    assert go_module._hypergeometric_upper_tail(2, 10, 4, 3) > float(stats.hypergeom.pmf(2, 10, 4, 3))


# ---------------------------------------------------------------------------
# semantic_similarity: documented range
# ---------------------------------------------------------------------------

TERM_IC = {"A": 0.5, "B": 3.0, "C": 5.0}  # information content in bits
HIERARCHY = {"B": {"A"}, "C": {"A", "B"}}


def test_resnik_is_the_unbounded_information_content_of_the_common_ancestor() -> None:
    assert semantic_similarity("B", "C", TERM_IC, HIERARCHY, method="resnik") == pytest.approx(3.0)


def test_lin_and_jiang_conrath_stay_between_0_and_1() -> None:
    rng = random.Random(7)
    for _ in range(200):
        # information content never decreases from a term to its descendants: A <- B <- C
        term_ic = dict(zip("ABC", sorted(rng.uniform(0, 8) for _ in range(3))))
        for method in ("lin", "jiang-conrath"):
            assert 0.0 <= semantic_similarity("B", "C", term_ic, HIERARCHY, method=method) <= 1.0


def test_the_documentation_says_resnik_is_not_bounded_by_one() -> None:
    returns = re.sub(r"\s+", " ", (semantic_similarity.__doc__ or "").split("Returns:")[1].split("Raises:")[0])

    assert re.search(r"resnik.*not bounded by 1", returns, re.IGNORECASE), returns
