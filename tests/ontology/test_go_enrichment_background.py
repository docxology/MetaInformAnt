"""Regression tests: ``enrich_genes`` counts only the annotated genes that are in the background, for both methods.

With an explicit ``background`` that does not contain every annotated gene, the Fisher branch built its 2x2 table from
the background-restricted sets, but the term size and the expected overlap used the *whole* annotation set, and the
hypergeometric branch took its success count ``K`` from that unrestricted size. A term listing genes the background
does not contain therefore got a ``K`` larger than the urn can hold, which understated its enrichment, and
``method="hypergeometric"`` disagreed with ``method="fisher"`` although both test the same upper tail.

Both methods now work on the annotations intersected with the background. The oracle is ``scipy.stats.fisher_exact``
(on the 2x2 table of the restricted sets) and ``scipy.stats.hypergeom.sf``, both fed with sets built in the test. No
test doubles.
"""

from __future__ import annotations

import random

import pytest
from scipy import stats

from metainformant.ontology.core.go import enrich_genes

UNIVERSE = [f"G{i:02d}" for i in range(60)]
BACKGROUND = UNIVERSE[:40]  # the annotations below also name genes G40..G59, which the background lacks


def _scenario(seed: int) -> tuple[list[str], dict[str, set[str]]]:
    rng = random.Random(seed)
    query = rng.sample(BACKGROUND, 12)
    annotations: dict[str, set[str]] = {}
    for index in range(8):
        inside = rng.sample(query, rng.randint(1, 6)) + rng.sample(BACKGROUND, rng.randint(2, 10))
        outside = rng.sample(UNIVERSE[40:], rng.randint(3, 15))  # annotated, but not in the background
        annotations[f"T{index}"] = set(inside) | set(outside)
    return query, annotations


def _oracle(query: list[str], background: list[str], annotations: dict[str, set[str]]) -> dict[str, dict[str, float]]:
    universe = set(background)
    drawn = set(query) & universe
    expected = {}
    for term, genes in annotations.items():
        members = genes & universe
        overlap = len(drawn & members)
        if overlap == 0:
            continue
        table = [
            [overlap, len(drawn) - overlap],
            [len(members) - overlap, len(universe) - len(drawn) - len(members) + overlap],
        ]
        expected[term] = {
            "fisher": float(stats.fisher_exact(table, alternative="greater")[1]),
            "hypergeom": float(stats.hypergeom.sf(overlap - 1, len(universe), len(members), len(drawn))),
            "term_size": len(members),
            "expected": len(members) * len(drawn) / len(universe),
        }
    return expected


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_both_methods_agree_with_the_oracle_for_an_explicit_background(seed: int) -> None:
    query, annotations = _scenario(seed)
    oracle = _oracle(query, BACKGROUND, annotations)
    assert oracle, "scenario must test at least one term"

    for method in ("fisher", "hypergeometric"):
        frame = enrich_genes(query, BACKGROUND, annotations, method=method)

        assert set(frame["go_term"]) == set(oracle), method
        for row in frame.itertuples():
            assert row.p_value == pytest.approx(oracle[row.go_term]["fisher"], rel=1e-9, abs=1e-300), (
                method,
                row.go_term,
            )
            assert row.p_value == pytest.approx(oracle[row.go_term]["hypergeom"], rel=1e-9, abs=1e-300), (
                method,
                row.go_term,
            )


def test_the_two_methods_return_the_same_frame_for_an_explicit_background() -> None:
    query, annotations = _scenario(4)

    fisher = enrich_genes(query, BACKGROUND, annotations, method="fisher").sort_values("go_term").reset_index(drop=True)
    hyper = (
        enrich_genes(query, BACKGROUND, annotations, method="hypergeometric")
        .sort_values("go_term")
        .reset_index(drop=True)
    )

    assert list(fisher["go_term"]) == list(hyper["go_term"])
    assert list(hyper["p_value"]) == pytest.approx(list(fisher["p_value"]), rel=1e-9, abs=1e-300)
    assert list(hyper["bonferroni_p"]) == pytest.approx(list(fisher["bonferroni_p"]), rel=1e-9, abs=1e-300)


@pytest.mark.parametrize("method", ["fisher", "hypergeometric"])
def test_term_size_and_expected_overlap_count_only_background_genes(method: str) -> None:
    query, annotations = _scenario(5)
    oracle = _oracle(query, BACKGROUND, annotations)

    frame = enrich_genes(query, BACKGROUND, annotations, method=method)

    for row in frame.itertuples():
        assert row.term_size == oracle[row.go_term]["term_size"]
        assert row.term_size < len(annotations[row.go_term]), "scenario annotates genes outside the background"
        assert row.expected == pytest.approx(oracle[row.go_term]["expected"])
        assert row.enrichment == pytest.approx(row.observed / oracle[row.go_term]["expected"])


@pytest.mark.parametrize("method", ["fisher", "hypergeometric"])
def test_a_term_annotated_only_outside_the_background_is_not_tested(method: str) -> None:
    query, annotations = _scenario(6)
    annotations["OUTSIDE"] = set(UNIVERSE[40:50])

    frame = enrich_genes(query, BACKGROUND, annotations, method=method)

    assert "OUTSIDE" not in set(frame["go_term"])


def test_without_an_explicit_background_the_methods_still_agree() -> None:
    query, annotations = _scenario(7)
    default_universe = sorted(set().union(*annotations.values()))
    oracle = _oracle(query, default_universe, annotations)

    for method in ("fisher", "hypergeometric"):
        frame = enrich_genes(query, None, annotations, method=method)

        assert {row.go_term: row.p_value for row in frame.itertuples()} == pytest.approx(
            {term: values["fisher"] for term, values in oracle.items()}, rel=1e-9, abs=1e-300
        )
