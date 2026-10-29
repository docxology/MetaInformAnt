"""``functional_enrichment_ppi`` tests the hypergeometric upper tail of a consistent population and honours filters.

The function summed ``math.comb`` terms in a loop guarded by ``except (ValueError, OverflowError): p_value = 1e-10``.
``math.comb`` works on arbitrary integers and ``int / int`` only overflows for a quotient above the float range, so
that branch (a made-up "very significant" p-value) could never run; it is gone. The same function had defects that the
dead branch hid:

* the draw ``n`` was ``len(protein_list)`` (duplicates, and proteins outside the background, counted) and the
  population ``N`` was ``len(background_proteins)`` (duplicates counted), while the overlap only counted distinct
  proteins inside the background, so the 2x2 table did not add up and the p-value was wrong;
* when ``n`` came out larger than ``N`` the annotation was skipped without a word;
* ``min_overlap`` was compared with the number of test proteins carrying the annotation, not with the overlap that is
  tested;
* ``p_value <= max_p_value or enrichment_ratio >= 1.0`` is true for every non-depleted annotation (a p-value below the
  threshold implies a ratio above 1), so ``max_p_value`` filtered nothing;
* a ``None`` annotation raised ``TypeError``.

The oracle is ``scipy.stats.hypergeom.sf`` / ``scipy.stats.fisher_exact`` fed with sets built in the test. No test
doubles.
"""

from __future__ import annotations

import math
import random

import pytest

stats = pytest.importorskip("scipy.stats")

from metainformant.networks.interaction.ppi import ProteinNetwork, functional_enrichment_ppi  # noqa: E402

ANNOTATED: list[tuple[str, object]] = [
    ("A1", "DNA repair"),
    ("A2", "DNA repair"),
    ("A3", "DNA repair"),
    ("B1", "cell cycle"),
    ("B2", "cell cycle"),
    ("C1", "metabolism"),
]


def _network(annotated: list[tuple[str, object]] = ANNOTATED) -> ProteinNetwork:
    """A chain network over exactly the given proteins, each carrying its ``function`` annotation."""
    network = ProteinNetwork()
    for protein, function in annotated:
        network.add_protein_metadata(protein, function=function)
    for (first, _), (second, _) in zip(annotated, annotated[1:]):
        network.add_interaction(first, second, 0.9, ["experimental"])
    return network


def _upper_tail(overlap: int, population: int, successes: int, draws: int) -> float:
    return float(stats.hypergeom.sf(overlap - 1, population, successes, draws))


class TestPopulation:
    def test_duplicates_and_outsiders_do_not_enter_the_population_or_the_draw(self) -> None:
        network = _network()
        background = ["A1", "A2", "A3", "A3", "B1", "B2", "C1"]  # A3 listed twice: six distinct proteins
        test_list = ["A1", "A2", "A3", "X1", "X2", "A1"]  # three distinct proteins in the background, two outside

        result = functional_enrichment_ppi(test_list, network, background_proteins=background, min_overlap=2)

        # N=6, K=3, n=3, k=3: C(3,3) / C(6,3) = 0.05; the old N=7, n=6 gave 4/7 = 0.571
        assert set(result) == {"DNA repair"}
        row = result["DNA repair"]
        assert row["p_value"] == pytest.approx(0.05, rel=1e-12, abs=1e-300)
        assert row["p_value"] == pytest.approx(_upper_tail(3, 6, 3, 3), rel=1e-12, abs=1e-300)
        assert row["count"] == 3
        assert row["background_count"] == 3
        assert row["enrichment_ratio"] == pytest.approx((3 / 3) / (3 / 6))
        assert row["proteins"] == ["A1", "A2", "A3"]

    def test_p_values_match_the_hypergeometric_oracle_on_random_networks(self) -> None:
        rng = random.Random(20261028)
        functions = ["f1", "f2", "f3", "f4", "f5"]
        checked = 0
        for _ in range(15):
            annotated = [(f"P{i}", rng.choice(functions)) for i in range(rng.randint(40, 90))]
            proteins = [protein for protein, _ in annotated]
            network = _network(annotated)
            background = rng.sample(proteins, rng.randint(25, len(proteins))) + [f"Q{i}" for i in range(3)]
            test_list = rng.sample(proteins, rng.randint(8, 30)) + ["OUTSIDE1", "OUTSIDE2"]
            test_list += test_list[:3]  # duplicates

            result = functional_enrichment_ppi(
                test_list, network, background_proteins=background + background[:4], min_overlap=1, max_p_value=1.0
            )

            population = set(background)
            drawn = set(test_list) & population
            for function, row in result.items():
                members = {protein for protein, label in annotated if label == function}
                successes = len(members & population)
                overlap = len(members & drawn)
                assert row["count"] == overlap
                assert row["background_count"] == successes
                expected = _upper_tail(overlap, len(population), successes, len(drawn))
                assert row["p_value"] == pytest.approx(expected, rel=1e-10, abs=1e-300)
                table = [
                    [overlap, len(drawn) - overlap],
                    [successes - overlap, len(population) - successes - len(drawn) + overlap],
                ]
                fisher_p = float(stats.fisher_exact(table, alternative="greater")[1])
                assert row["p_value"] == pytest.approx(fisher_p, rel=1e-9, abs=1e-300)
                checked += 1
        assert checked >= 40

    def test_a_test_list_longer_than_the_background_list_is_still_tested(self) -> None:
        network = _network()
        test_list = ["A1", "A2", "A3", "B1", "X1", "X2", "X3", "X4", "X5", "X6", "X7"]  # 11 names, N = 4

        result = functional_enrichment_ppi(
            test_list, network, background_proteins=["A1", "A2", "B1", "C1"], min_overlap=2, max_p_value=1.0
        )

        # the old guard (n <= N) skipped every annotation; the tested draw is {A1, A2, B1} from N=4, K=2 repairs
        assert set(result) == {"DNA repair"}
        assert result["DNA repair"]["count"] == 2
        assert result["DNA repair"]["p_value"] == pytest.approx(_upper_tail(2, 4, 2, 3), rel=1e-12, abs=1e-300)


class TestFilters:
    def test_min_overlap_counts_the_proteins_that_are_actually_tested(self) -> None:
        network = _network()
        # three DNA-repair proteins in the list but only A1 in the background: the overlap is 1
        kwargs = {"background_proteins": ["A1", "B1", "B2", "C1"], "max_p_value": 1.0}

        assert functional_enrichment_ppi(["A1", "A2", "A3"], network, min_overlap=2, **kwargs) == {}
        assert set(functional_enrichment_ppi(["A1", "A2", "A3"], network, min_overlap=1, **kwargs)) == {"DNA repair"}

    def test_max_p_value_is_honoured(self) -> None:
        network = _network()
        drawn = ["A1", "A2", "B1", "B2"]
        # DNA repair: 2 of its 3 proteins in a draw of 4 from 6, P(X >= 2) = 12/15 = 0.8 (ratio 1.0);
        # cell cycle: both of its 2 proteins, P(X >= 2) = 6/15 = 0.4 (ratio 1.5). Both have a ratio of at least 1,
        # which used to keep them whatever the threshold.
        assert functional_enrichment_ppi(drawn, network, max_p_value=0.05) == {}
        assert set(functional_enrichment_ppi(drawn, network, max_p_value=0.5)) == {"cell cycle"}

        everything = functional_enrichment_ppi(drawn, network, max_p_value=1.0)
        assert set(everything) == {"DNA repair", "cell cycle"}
        assert everything["DNA repair"]["p_value"] == pytest.approx(0.8)
        assert everything["cell cycle"]["p_value"] == pytest.approx(0.4)

    def test_a_significant_annotation_survives_the_default_threshold(self) -> None:
        network = _network([(f"P{i}", "target" if i < 5 else "other") for i in range(60)])

        result = functional_enrichment_ppi([f"P{i}" for i in range(5)], network, min_overlap=2)

        # all 5 annotated proteins in a draw of 5 from 60: 1 / C(60, 5)
        assert set(result) == {"target"}
        assert result["target"]["p_value"] == pytest.approx(1 / math.comb(60, 5), rel=1e-12, abs=1e-300)

    def test_depleted_annotations_are_not_reported_but_enriched_ones_are(self) -> None:
        network = _network([(f"P{i}", "common" if i < 40 else "rare") for i in range(60)])
        draw = [f"P{i}" for i in range(40, 60)] + ["P0", "P1"]  # mostly rare, only two common proteins

        result = functional_enrichment_ppi(draw, network, min_overlap=2, max_p_value=0.05)

        assert set(result) == {"rare"}
        assert result["rare"]["p_value"] == pytest.approx(_upper_tail(20, 60, 20, 22), rel=1e-10, abs=1e-300)


class TestBigPopulationsAndAnnotations:
    def test_a_large_population_gets_the_exact_tiny_p_value_not_a_made_up_one(self) -> None:
        network = _network([(f"P{i}", "marked" if i < 1500 else "other") for i in range(3000)])

        result = functional_enrichment_ppi([f"P{i}" for i in range(400)], network, min_overlap=2)

        # 400 of 400 draws are marked although only half the background is: about 1e-121, never 1e-10
        expected = _upper_tail(400, 3000, 1500, 400)
        assert 0.0 < expected < 1e-100
        assert result["marked"]["p_value"] == pytest.approx(expected, rel=1e-9, abs=1e-300)

    def test_a_missing_or_none_annotation_is_no_annotation(self) -> None:
        network = _network([("A1", "x"), ("A2", "x"), ("A3", "x"), ("B1", None), ("B2", None), ("C1", "y")])

        result = functional_enrichment_ppi(["A1", "A2", "A3", "B1"], network, min_overlap=2, max_p_value=1.0)

        assert set(result) == {"x"}  # the old code raised "TypeError: 'NoneType' object is not iterable"

    def test_proteins_with_several_annotations_count_once_per_annotation(self) -> None:
        network = _network([("A1", ["x", "y", "x"]), ("A2", ["x", "y"]), ("A3", ["x"]), ("B1", ["y"]), ("B2", "z")])

        result = functional_enrichment_ppi(["A1", "A2", "A3"], network, min_overlap=2, max_p_value=1.0)

        assert result["x"]["count"] == 3 and result["x"]["background_count"] == 3
        assert result["y"]["count"] == 2 and result["y"]["background_count"] == 3
        assert result["x"]["p_value"] == pytest.approx(_upper_tail(3, 5, 3, 3), rel=1e-12, abs=1e-300)
        assert result["y"]["p_value"] == pytest.approx(_upper_tail(2, 5, 3, 3), rel=1e-12, abs=1e-300)

    def test_the_documented_result_keys_are_unchanged(self) -> None:
        result = functional_enrichment_ppi(["A1", "A2", "A3"], _network(), min_overlap=2, max_p_value=1.0)

        assert set(result["DNA repair"]) == {
            "count",
            "test_count",
            "background_count",
            "enrichment_ratio",
            "p_value",
            "proteins",
        }
