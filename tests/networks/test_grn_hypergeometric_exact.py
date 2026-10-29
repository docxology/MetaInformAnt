"""The GRN hypergeometric p-value is the exact upper tail, with or without SciPy, on a valid population.

``grn_inference._hypergeometric_pvalue`` called ``scipy.stats.hypergeom.sf`` when SciPy was importable and otherwise
fell back to a normal approximation, which is far off exactly where enrichment is tested (small counts, a small gene
set in a large population: 0.5 for 0.593, 7e-05 for 0.0066, and 0.0 for 3e-20). It also returned 1.0 for any call with a
zero in the population, including an observed overlap that cannot occur, and for a gene set larger than the target
universe (SciPy returns NaN there and the clamp ``max(0, min(1, nan))`` turned it into 1.0).

``score_regulators`` made the last case routine: it passed ``len(gene_set)`` as the number of successes in the
population. Gene-set members the GRN never targets (and duplicate entries) inflated it, so every p-value came out too
large, and it rounded the result to 8 decimals, so every p-value below 5e-9 was reported as exactly 0.0.

The p-value is now summed in log space, exactly, with no SciPy branch. The oracles are ``scipy.stats.hypergeom.sf``,
``scipy.stats.fisher_exact`` and a ``math.comb`` integer sum, none of which the implementation uses. The no-SciPy
cases run the real module in a subprocess in which SciPy cannot be imported. No test doubles.
"""

from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

stats = pytest.importorskip("scipy.stats")

import metainformant  # noqa: E402
from metainformant.networks.analysis.pathway import _hypergeometric_upper_tail  # noqa: E402
from metainformant.networks.regulatory.grn_inference import _hypergeometric_pvalue, score_regulators  # noqa: E402


def _scipy_upper_tail(k: int, n: int, K: int, N: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws), by scipy."""
    return float(stats.hypergeom.sf(k - 1, N, K, n))


def _integer_upper_tail(k: int, n: int, K: int, N: int) -> float:
    """The same tail from exact integer arithmetic, correctly rounded once."""
    lowest, highest = max(k, 0, n + K - N), min(K, n)
    numerator = sum(math.comb(K, x) * math.comb(N - K, n - x) for x in range(lowest, highest + 1))
    return float(Fraction(numerator, math.comb(N, n)))


def _subprocess_env() -> dict[str, str]:
    """Environment that makes the child import the same ``metainformant`` as this process."""
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    return {**os.environ, "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else []))}


class TestExactUpperTail:
    def test_hand_computed_tail(self) -> None:
        # N=10, K=4, n=3: P(X >= 2) = (C(4,2)*C(6,1) + C(4,3)*C(6,0)) / C(10,3) = (36 + 4) / 120 = 1/3
        assert _hypergeometric_pvalue(k=2, n=3, K=4, N=10) == pytest.approx(1 / 3, rel=1e-12, abs=1e-300)

    def test_matches_scipy_over_a_grid_of_parameters(self) -> None:
        checked = 0
        for N in (12, 40, 200):
            for K in sorted({0, 1, 5, N // 2, N}):
                for n in sorted({0, 1, 7, N // 3, N}):
                    for k in range(0, min(K, n) + 1):
                        assert _hypergeometric_pvalue(k, n, K, N) == pytest.approx(
                            _scipy_upper_tail(k, n, K, N), rel=1e-9, abs=1e-300
                        ), (k, n, K, N)
                        checked += 1
        assert checked > 300

    def test_matches_scipy_on_random_tables(self) -> None:
        rng = random.Random(20261028)
        for _ in range(300):
            population = rng.randint(1, 400)
            successes = rng.randint(0, population)
            draws = rng.randint(0, population)
            k = rng.randint(-1, min(successes, draws) + 2)

            assert _hypergeometric_pvalue(k, draws, successes, population) == pytest.approx(
                _scipy_upper_tail(k, draws, successes, population), rel=1e-9, abs=1e-300
            ), (k, draws, successes, population)

    def test_is_the_one_sided_fisher_exact_p_value(self) -> None:
        for k, n, K, N in [(4, 15, 20, 100), (0, 15, 20, 100), (9, 12, 10, 60), (5, 5, 30, 31), (3, 10, 4, 12)]:
            table = [[k, n - k], [K - k, N - K - n + k]]
            _, fisher_p = stats.fisher_exact(table, alternative="greater")

            assert _hypergeometric_pvalue(k, n, K, N) == pytest.approx(float(fisher_p), rel=1e-9, abs=1e-300)

    @pytest.mark.parametrize(
        ("k", "n", "K", "N"),
        [
            (300, 2000, 1000, 20000),  # 4e-74
            (40, 60, 45, 20000),  # 4e-103
            (12, 40, 30, 2000),  # 5e-14
            (8, 10, 15, 5000),  # 3e-20
            (150, 3000, 900, 12000),  # a depleted-looking but non-trivial tail
        ],
    )
    def test_stays_accurate_for_large_populations_and_tiny_tails(self, k: int, n: int, K: int, N: int) -> None:
        got = _hypergeometric_pvalue(k, n, K, N)

        assert got == pytest.approx(_scipy_upper_tail(k, n, K, N), rel=1e-9, abs=1e-300)
        assert got == pytest.approx(_integer_upper_tail(k, n, K, N), rel=1e-9, abs=1e-300)

    def test_agrees_with_the_exact_integer_helper_used_by_pathway_enrichment(self) -> None:
        for k, population, successes, draws in [(3, 50, 12, 9), (7, 300, 40, 25), (0, 20, 6, 6), (2, 30, 30, 4)]:
            assert _hypergeometric_pvalue(k, draws, successes, population) == pytest.approx(
                _hypergeometric_upper_tail(k, population, successes, draws), rel=1e-9
            )

    def test_a_zero_or_negative_observed_overlap_is_certain(self) -> None:
        assert _hypergeometric_pvalue(0, 5, 3, 20) == 1.0
        assert _hypergeometric_pvalue(-2, 5, 3, 20) == 1.0
        assert _hypergeometric_pvalue(0, 0, 10, 100) == 1.0
        assert _hypergeometric_pvalue(0, 10, 0, 100) == 1.0

    def test_an_overlap_forced_by_the_pigeonhole_principle_is_certain(self) -> None:
        # 8 draws from 10 items of which 5 are successes: at least 3 successes must be drawn
        assert _hypergeometric_pvalue(3, 8, 5, 10) == 1.0
        assert _hypergeometric_pvalue(2, 8, 5, 10) == 1.0

    def test_an_overlap_that_cannot_occur_has_probability_zero_not_one(self) -> None:
        # the old code answered 1.0 whenever n, K or N was zero, whatever k was
        assert _hypergeometric_pvalue(k=5, n=0, K=10, N=100) == 0.0
        assert _hypergeometric_pvalue(k=5, n=10, K=0, N=100) == 0.0
        assert _hypergeometric_pvalue(k=6, n=8, K=5, N=20) == 0.0
        assert _hypergeometric_pvalue(k=1, n=0, K=0, N=0) == 0.0

    @pytest.mark.parametrize(
        ("k", "n", "K", "N"),
        [
            (3, 10, 50, 40),  # more successes than the population
            (1, 12, 2, 10),  # more draws than the population
            (0, 2, -1, 5),
            (0, -1, 2, 5),
            (0, 0, 0, -1),
        ],
    )
    def test_parameters_that_are_not_a_draw_from_a_population_raise(self, k: int, n: int, K: int, N: int) -> None:
        with pytest.raises(ValueError, match="inconsistent hypergeometric parameters"):
            _hypergeometric_pvalue(k, n, K, N)

    def test_result_is_always_a_probability(self) -> None:
        for k in range(-1, 25):
            assert 0.0 <= _hypergeometric_pvalue(k, 20, 30, 100) <= 1.0


NO_SCIPY_SCRIPT = """
import json, sys

sys.modules["scipy"] = None  # any `import scipy...` now raises ImportError, as on a machine without scipy

from metainformant.networks.regulatory import grn_inference

assert grn_inference.HAS_SCIPY is False, "scipy was importable"
cases = json.load(sys.stdin)
print(json.dumps([grn_inference._hypergeometric_pvalue(*case) for case in cases]))
"""


class TestWithoutScipy:
    CASES = [
        (5, 20, 50, 200),
        (30, 30, 50, 200),
        (3, 5, 10, 100),
        (2, 3, 4, 10),
        (12, 40, 30, 2000),
        (8, 10, 15, 5000),
        (1, 10, 5, 1000),
        (300, 2000, 1000, 20000),
        (0, 10, 5, 1000),
        (7, 6, 20, 100),
    ]

    def test_the_exact_tail_does_not_depend_on_scipy(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-c", NO_SCIPY_SCRIPT],
            input=json.dumps(self.CASES),
            capture_output=True,
            text=True,
            env=_subprocess_env(),
            timeout=300,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
        observed = json.loads(completed.stdout.strip().splitlines()[-1])

        # the normal approximation answered 0.5 for 0.593, 7e-05 for 0.0066, 0.0 for 3e-20 and 0.14 for 1/3
        assert observed == pytest.approx(
            [_integer_upper_tail(k, n, K, N) for k, n, K, N in self.CASES], rel=1e-9, abs=1e-300
        )


class TestScoreRegulators:
    @staticmethod
    def _grn(edges: list[tuple[str, str]]) -> dict:
        return {"edges": [{"source": tf, "target": target, "weight": 1.0} for tf, target in edges]}

    @staticmethod
    def _fisher(grn: dict, tf: str, gene_set: list[str]) -> float:
        """One-sided Fisher p-value on the universe of distinct targets, built independently of the code under test."""
        targets_of = {e["target"] for e in grn["edges"] if e["source"] == tf}
        universe = {e["target"] for e in grn["edges"]}
        wanted = {g.lower() for g in gene_set}
        successes = {t for t in universe if t.lower() in wanted}
        k, n, K, N = len(targets_of & successes), len(targets_of), len(successes), len(universe)
        return float(stats.fisher_exact([[k, n - k], [K - k, N - K - n + k]], alternative="greater")[1])

    def _three_tf_grn(self) -> dict:
        targets = [f"T{i:02d}" for i in range(30)]
        edges = (
            [("TF_A", t) for t in targets[:10]] + [("TF_B", t) for t in targets] + [("TF_C", t) for t in targets[:3]]
        )
        return self._grn(edges)

    def test_gene_set_members_outside_the_grn_and_duplicates_do_not_change_the_p_values(self) -> None:
        grn = self._three_tf_grn()
        inside = [f"T{i:02d}" for i in range(12)]
        padded = inside + [f"OUTSIDE{i}" for i in range(8)] + inside[:2] + [g.lower() for g in inside[:3]]

        clean = {row["tf"]: row["enrichment_p"] for row in score_regulators(grn, inside)}
        noisy = {row["tf"]: row["enrichment_p"] for row in score_regulators(grn, padded)}

        assert noisy == pytest.approx(clean, rel=1e-12, abs=1e-300)
        for tf, p in noisy.items():
            assert p == pytest.approx(self._fisher(grn, tf, inside), rel=1e-9, abs=1e-300)
        # the old len(gene_set) made these 0.0215 and 0.379
        assert noisy["TF_A"] == pytest.approx(2.1967038458792583e-06, rel=1e-9, abs=1e-300)
        assert noisy["TF_C"] == pytest.approx(0.0541871921182266, rel=1e-9, abs=1e-300)

    def test_a_gene_set_much_larger_than_the_target_universe_still_gives_the_fisher_p_value(self) -> None:
        grn = self._three_tf_grn()
        huge = [f"T{i:02d}" for i in range(12)] + [f"ELSEWHERE{i}" for i in range(500)]

        got = {row["tf"]: row["enrichment_p"] for row in score_regulators(grn, huge)}

        for tf, p in got.items():
            assert p == pytest.approx(self._fisher(grn, tf, huge), rel=1e-9, abs=1e-300)
        assert got["TF_A"] < 1e-5  # the old successes count of 512 > 30 targets gave SciPy NaN, clamped to 1.0

    def test_a_strongly_enriched_tf_keeps_its_tiny_p_value(self) -> None:
        genes = [f"G{i:03d}" for i in range(400)]
        grn = self._grn([("TF_X", g) for g in genes[:30]] + [("TF_Y", g) for g in genes[30:]])

        rows = {row["tf"]: row for row in score_regulators(grn, genes[:30])}

        # every one of TF_X's 30 targets is in the gene set: C(30, 30) / C(400, 30), 7e-46, not round(p, 8) == 0.0
        assert rows["TF_X"]["enrichment_p"] == pytest.approx(
            float(stats.hypergeom.sf(29, 400, 30, 30)), rel=1e-9, abs=1e-300
        )
        assert rows["TF_X"]["enrichment_p"] > 0.0
        assert rows["TF_Y"]["enrichment_p"] == 1.0

    def test_p_values_match_fisher_on_random_networks(self) -> None:
        rng = random.Random(7)
        genes = [f"g{i}" for i in range(80)]
        checked = 0
        for _ in range(12):
            edges = [(f"TF{rng.randint(0, 5)}", rng.choice(genes)) for _ in range(rng.randint(30, 120))]
            grn = self._grn(edges)
            gene_set = rng.sample(genes, rng.randint(3, 30)) + ["not_in_grn_1", "not_in_grn_2"]

            for row in score_regulators(grn, gene_set):
                assert row["enrichment_p"] == pytest.approx(
                    self._fisher(grn, row["tf"], gene_set), rel=1e-9, abs=1e-300
                )
                checked += 1
        assert checked >= 40

    def test_results_are_sorted_by_p_value_and_an_empty_gene_set_is_never_significant(self) -> None:
        grn = self._three_tf_grn()
        ranked = [row["enrichment_p"] for row in score_regulators(grn, [f"T{i:02d}" for i in range(12)])]

        assert ranked == sorted(ranked)
        assert [row["enrichment_p"] for row in score_regulators(grn, [])] == [1.0, 1.0, 1.0]
