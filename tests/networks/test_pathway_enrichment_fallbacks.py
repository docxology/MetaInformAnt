"""Regression tests for three defects in ``metainformant.networks.analysis.pathway``.

* ``pathway_enrichment`` computed ``max(0.001, 1 / (overlap / expected))`` as its p-value whenever scipy could not be
  imported. That is not a Fisher or hypergeometric test (it is floored at 0.001, ignores the pathway and background
  sizes beyond the expected count, and returns 1.0 for every overlap at or below the expectation). The fallback now
  computes the exact upper-tail probability, and ``correction="fdr"`` no longer silently returns the raw p-values when
  scipy is missing (Benjamini-Hochberg is pure Python).
* ``pathway_enrichment`` ran any ``method`` other than ``"fisher"`` as the hypergeometric test, so a misspelt method
  produced results instead of an error.
* ``create_pathway_network`` recovered the two pathway names from the ``"{a}_{b}"`` similarity key with
  ``key.split("_", 1)``, which cuts every pathway ID that itself contains an underscore (``KEGG_GLYCOLYSIS``,
  ``HALLMARK_APOPTOSIS`` ...) in the wrong place, so the edges landed on nodes that do not exist.

The oracles are ``scipy.stats.hypergeom.sf`` and ``scipy.stats.fisher_exact`` (an implementation independent of the
code under test) and ``scipy.stats.false_discovery_control``; the graph assertions are written from the Jaccard
definition. The no-scipy cases run the real module in a subprocess in which scipy cannot be imported. No test doubles.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest
from scipy import stats

import metainformant
from metainformant.networks.analysis import pathway as pathway_module
from metainformant.networks.analysis.pathway import (
    PathwayNetwork,
    create_pathway_network,
    network_enrichment_analysis,
    pathway_enrichment,
)

BACKGROUND = [f"g{i:02d}" for i in range(40)]


def _subprocess_env() -> dict[str, str]:
    """Environment that makes the child import the same ``metainformant`` as this process."""
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    return {**os.environ, "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else []))}


NO_SCIPY_SCRIPT = """
import json, sys

sys.modules["scipy"] = None  # any `import scipy...` now raises ImportError, as on a machine without scipy

from metainformant.networks.analysis import pathway

payload = json.load(sys.stdin)
assert pathway.HAS_SCIPY is False, "scipy was importable"
network = pathway.PathwayNetwork("n", payload["pathways"])
out = {}
for method in ("fisher", "hypergeometric"):
    for correction in (None, "bonferroni", "fdr"):
        results = pathway.pathway_enrichment(
            payload["query"], network, background_genes=payload["background"],
            method=method, correction=correction,
        )
        out[f"{method}|{correction}"] = {
            pid: [r["p_value"], r.get("corrected_p_value")] for pid, r in results.items()
        }
print(json.dumps(out))
"""


def _hypergeometric_oracle(query: list[str], background: list[str], pathway: list[str]) -> float:
    universe = set(background)
    drawn = set(query) & universe
    members = set(pathway) & universe
    overlap = len(drawn & members)
    return float(stats.hypergeom.sf(overlap - 1, len(universe), len(members), len(drawn)))


def _random_scenario(seed: int) -> tuple[list[str], dict[str, list[str]]]:
    rng = random.Random(seed)
    query = rng.sample(BACKGROUND, rng.randint(4, 14))
    pathways = {}
    for index in range(8):
        # half of the pathways are seeded with query genes so some overlaps are enriched and some are not
        base = rng.sample(query, rng.randint(1, min(len(query), 5))) if index % 2 == 0 else []
        pathways[f"P{index}"] = sorted(set(base + rng.sample(BACKGROUND, rng.randint(3, 12))))
    return query, pathways


# ---------------------------------------------------------------------------
# p-value without scipy
# ---------------------------------------------------------------------------


def test_upper_tail_matches_scipy_on_random_tables() -> None:
    rng = random.Random(20261027)
    for _ in range(400):
        population = rng.randint(1, 300)
        successes = rng.randint(0, population)
        draws = rng.randint(0, population)
        k = rng.randint(-1, min(successes, draws) + 2)

        reference = float(stats.hypergeom.sf(k - 1, population, successes, draws))

        assert pathway_module._hypergeometric_upper_tail(k, population, successes, draws) == pytest.approx(
            reference, rel=1e-12, abs=1e-300
        ), (k, population, successes, draws)


def test_upper_tail_is_the_one_sided_fisher_exact_p_value() -> None:
    # 2x2 table [[k, draws - k], [successes - k, population - successes - draws + k]]
    for k, population, successes, draws in [(3, 10, 4, 3), (5, 40, 12, 9), (0, 20, 6, 6), (7, 100, 30, 15)]:
        table = [[k, draws - k], [successes - k, population - successes - draws + k]]
        _, fisher_p = stats.fisher_exact(table, alternative="greater")

        assert pathway_module._hypergeometric_upper_tail(k, population, successes, draws) == pytest.approx(
            fisher_p, rel=1e-12
        )


def test_upper_tail_handles_large_populations_without_underflow_or_overflow() -> None:
    # an enriched overlap (tiny p) and a depleted one (p close to 1), both far outside float-only territory
    assert pathway_module._hypergeometric_upper_tail(300, 20000, 1000, 2000) == pytest.approx(
        float(stats.hypergeom.sf(299, 20000, 1000, 2000)), rel=1e-9
    )
    assert pathway_module._hypergeometric_upper_tail(1, 20000, 10000, 10000) == pytest.approx(1.0)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_enrichment_without_scipy_matches_the_scipy_tests(seed: int) -> None:
    query, pathways = _random_scenario(seed)
    payload = {"query": query, "background": BACKGROUND, "pathways": pathways}

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

    expected_p = {pid: _hypergeometric_oracle(query, BACKGROUND, genes) for pid, genes in pathways.items()}
    # only pathways that share a gene with the query are reported (min_overlap=1)
    reported = {pid for pid, genes in pathways.items() if set(query) & set(genes)}
    assert reported, "scenario must report at least one pathway"

    for method in ("fisher", "hypergeometric"):
        for correction in (None, "bonferroni", "fdr"):
            results = observed[f"{method}|{correction}"]
            assert set(results) == reported
            raw = {pid: pair[0] for pid, pair in results.items()}
            assert raw == pytest.approx({pid: expected_p[pid] for pid in reported}, rel=1e-12, abs=1e-300)

            ordered = list(results)
            corrected = [results[pid][1] for pid in ordered]
            p_values = [expected_p[pid] for pid in ordered]
            if correction is None:
                assert corrected == [None] * len(ordered)
            elif correction == "bonferroni":
                assert corrected == pytest.approx([min(p * len(ordered), 1.0) for p in p_values], rel=1e-12)
            else:
                # the old code returned the raw p-values here when scipy was missing
                expected_q = [float(q) for q in stats.false_discovery_control(p_values, method="bh")]
                assert corrected == pytest.approx(expected_q, rel=1e-12, abs=1e-15)


def test_enrichment_p_values_with_scipy_match_the_oracle_for_both_methods() -> None:
    query, pathways = _random_scenario(11)
    network = PathwayNetwork("n", pathways)

    for method in ("fisher", "hypergeometric"):
        results = pathway_enrichment(query, network, background_genes=BACKGROUND, method=method, correction=None)

        assert results
        for pid, result in results.items():
            assert result["p_value"] == pytest.approx(
                _hypergeometric_oracle(query, BACKGROUND, pathways[pid]), rel=1e-12, abs=1e-300
            )


# ---------------------------------------------------------------------------
# unknown method
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["chi2", "FISHER", "hypergeom", ""])
def test_unknown_method_raises_instead_of_running_as_hypergeometric(method: str) -> None:
    network = PathwayNetwork("n", {"P1": BACKGROUND[:6], "P2": BACKGROUND[3:9]})

    with pytest.raises(ValueError, match="Unknown method"):
        pathway_enrichment(BACKGROUND[:5], network, background_genes=BACKGROUND, method=method)


def test_network_enrichment_alias_rejects_an_unknown_method_too() -> None:
    network = PathwayNetwork("n", {"P1": BACKGROUND[:6]})

    with pytest.raises(ValueError, match="Unknown method"):
        network_enrichment_analysis(BACKGROUND[:5], network, background_genes=BACKGROUND, method="bogus")


def test_unknown_method_is_rejected_even_for_an_empty_network() -> None:
    with pytest.raises(ValueError, match="Unknown method"):
        pathway_enrichment(BACKGROUND[:5], PathwayNetwork("empty"), method="bogus")


# ---------------------------------------------------------------------------
# pathway IDs that contain underscores
# ---------------------------------------------------------------------------

UNDERSCORE_PATHWAYS = {
    "KEGG_GLYCOLYSIS": ["G1", "G2", "G3", "G4"],
    "KEGG_TCA_CYCLE": ["G3", "G4", "G5", "G6"],
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": ["G5", "G6", "G7", "G8"],
    "UNRELATED": ["U1", "U2", "U3", "U4"],
}


def _jaccard(a: list[str], b: list[str]) -> float:
    return len(set(a) & set(b)) / len(set(a) | set(b))


def _overlap_coefficient(a: list[str], b: list[str]) -> float:
    return len(set(a) & set(b)) / min(len(set(a)), len(set(b)))


@pytest.mark.parametrize(
    ("method", "measure"), [("jaccard", _jaccard), ("overlap", _overlap_coefficient)], ids=["jaccard", "overlap"]
)
def test_pathway_network_edges_join_the_real_pathway_ids(method: str, measure) -> None:
    threshold = 0.2
    names = list(UNDERSCORE_PATHWAYS)
    expected = {}
    for i, first in enumerate(names):
        for second in names[i + 1 :]:
            similarity = measure(UNDERSCORE_PATHWAYS[first], UNDERSCORE_PATHWAYS[second])
            if similarity >= threshold:
                expected[frozenset((first, second))] = similarity
    assert len(expected) == 2, "scenario must contain edges whose ids contain underscores"

    network = create_pathway_network(UNDERSCORE_PATHWAYS, similarity_threshold=threshold, method=method)

    assert set(network.nodes()) == set(UNDERSCORE_PATHWAYS), "no phantom nodes from a mis-split key"
    assert {frozenset(edge) for edge in network.edges()} == set(expected)
    for edge, similarity in expected.items():
        first, second = tuple(edge)
        assert network[first][second]["similarity"] == pytest.approx(similarity)


def test_pathway_network_edges_are_not_confused_by_ids_that_prefix_each_other() -> None:
    pathways = {"A": ["g1", "g2"], "A_B": ["g1", "g2", "g3"], "A_B_C": ["g7", "g8"], "C": ["g7", "g8", "g9"]}

    network = create_pathway_network(pathways, similarity_threshold=0.5)

    assert {frozenset(edge) for edge in network.edges()} == {frozenset(("A", "A_B")), frozenset(("A_B_C", "C"))}
