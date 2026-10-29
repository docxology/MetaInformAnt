"""The permutation seed of ``active_module_detection`` and ``pathway_topology_analysis`` is an argument.

Regression tests for ``metainformant.multiomics.pathways.enrichment``. ``active_module_detection`` built
``random.Random(42)`` inside its module loop, so every module's permutation test used the same hard-coded seed and
a caller could neither vary it nor check that a p-value was stable across seeds; ``pathway_topology_analysis``
hard-coded the same seed and a fixed 1000 permutations. Both now take ``random_state`` (default 42, which keeps the
previous results exactly) and ``active_module_detection`` also refuses arguments that made its answer meaningless:
``n_permutations=0`` reported a p-value of 1.0 for a test that never ran, a negative p-value was scored as the most
significant node, and an ``alpha`` outside (0, 1] was clamped silently.

The golden p-values below were computed with the implementation before the change (fixed seed 42). Real
computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import inspect
import math
import random
from typing import Any

import numpy as np
import pytest

from metainformant.multiomics.pathways.enrichment import active_module_detection, pathway_topology_analysis

N_PERMUTATIONS = 300
OTHER_SEEDS = (0, 1, 7, 2024)


def _network_and_scores(seed: int = 11, n_background: int = 60) -> tuple[dict[str, list[str]], dict[str, float]]:
    """Three significant chains in a background where a good share of the nodes is significant as well.

    Random node sets therefore often score as high as a chain, which keeps the permutation p-values away from
    their minimum of 1 / (n_permutations + 1) and makes them depend on the random stream.
    """
    rng = random.Random(seed)
    network: dict[str, list[str]] = {}
    scores: dict[str, float] = {}
    for prefix, p_values in {"A": [0.03, 0.02, 0.04, 0.01], "B": [0.045, 0.03, 0.02], "C": [0.008, 0.04]}.items():
        nodes = [f"{prefix}{i}" for i in range(len(p_values))]
        for u, v in zip(nodes, nodes[1:]):
            network.setdefault(u, []).append(v)
            network.setdefault(v, []).append(u)
        scores.update(zip(nodes, p_values))
    for i in range(n_background):
        name = f"G{i:02d}"
        network[name] = []
        scores[name] = round(rng.choice([rng.uniform(0.001, 0.05), rng.uniform(0.05, 1.0)]), 4)
    return network, scores


def _p_by_module(modules: list[dict[str, Any]]) -> dict[tuple[str, ...], float]:
    return {tuple(m["module_genes"]): m["p_value"] for m in modules}


GRAPH = {"S": ["T1", "T2"], "T1": ["U"], "T2": ["U"], "U": ["V"], "V": []}
GENE_SCORES = {"S": 0.001, "T1": 0.02, "T2": 0.4, "U": 0.003, "V": 0.9, "W": 0.2, "X": 0.7}


class TestActiveModuleRandomState:
    def test_the_parameter_exists_and_defaults_to_the_previously_fixed_seed(self) -> None:
        assert inspect.signature(active_module_detection).parameters["random_state"].default == 42
        assert inspect.signature(pathway_topology_analysis).parameters["random_state"].default == 42

    def test_the_default_reproduces_the_results_of_the_fixed_seed_exactly(self) -> None:
        network, scores = _network_and_scores()

        modules = active_module_detection(network, scores, alpha=0.05, n_permutations=N_PERMUTATIONS)

        # p = (permutations scoring at least as high + 1) / (n_permutations + 1), computed before the change
        assert _p_by_module(modules) == {
            ("A0", "A1", "A2", "A3"): 34 / 301,
            ("C0", "C1"): 49 / 301,
            ("B0", "B1", "B2"): 74 / 301,
        }
        assert [m["module_genes"] for m in modules] == [["A0", "A1", "A2", "A3"], ["C0", "C1"], ["B0", "B1", "B2"]]

    def test_an_explicit_42_is_the_default(self) -> None:
        network, scores = _network_and_scores()

        default = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS)
        explicit = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS, random_state=42)

        assert default == explicit

    @pytest.mark.parametrize("seed", OTHER_SEEDS)
    def test_another_seed_changes_the_p_values_but_not_the_modules_or_their_scores(self, seed: int) -> None:
        network, scores = _network_and_scores()
        reference = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS)

        modules = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS, random_state=seed)

        assert [m["module_genes"] for m in modules] == [m["module_genes"] for m in reference]
        assert [m["module_score"] for m in modules] == [m["module_score"] for m in reference]
        assert _p_by_module(modules) != _p_by_module(reference)  # the seed reached the permutation test
        for module in modules:
            assert 1 / (N_PERMUTATIONS + 1) <= module["p_value"] <= 1.0
            assert (module["p_value"] * (N_PERMUTATIONS + 1)) == pytest.approx(
                round(module["p_value"] * (N_PERMUTATIONS + 1))
            )

    def test_one_seed_always_gives_the_same_answer(self) -> None:
        network, scores = _network_and_scores()

        first = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS, random_state=7)
        second = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS, random_state=7)

        assert first == second

    def test_the_seed_drives_python_random_through_a_generator_started_for_every_module(self) -> None:
        # independent oracle: redo each module's permutation test with random.Random(seed) and the node scores
        network, scores = _network_and_scores()
        alpha = 0.05
        threshold = -math.log(alpha)
        universe = list(scores)

        for seed in (42, 7):
            modules = active_module_detection(
                network, scores, alpha=alpha, n_permutations=N_PERMUTATIONS, random_state=seed
            )
            assert len(modules) == 3
            for module in modules:
                size = len(module["module_genes"])
                rng = random.Random(seed)
                at_least_as_high = 0
                for _ in range(N_PERMUTATIONS):
                    drawn = rng.sample(universe, size)
                    if (
                        sum(-math.log(max(scores[n], 1e-300)) - threshold for n in drawn)
                        >= module["module_score"] - 1e-9
                    ):
                        at_least_as_high += 1
                assert module["p_value"] == (at_least_as_high + 1) / (N_PERMUTATIONS + 1), (
                    seed,
                    module["module_genes"],
                )

    def test_none_gives_valid_p_values_and_the_same_modules(self) -> None:
        network, scores = _network_and_scores()
        reference = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS)

        modules = active_module_detection(network, scores, n_permutations=N_PERMUTATIONS, random_state=None)

        # the same modules; their order follows the (now unseeded) p-values, so it is not compared
        assert sorted(m["module_genes"] for m in modules) == sorted(m["module_genes"] for m in reference)
        assert all(1 / (N_PERMUTATIONS + 1) <= m["p_value"] <= 1.0 for m in modules)

    def test_numpy_integers_are_accepted_as_seeds(self) -> None:
        network, scores = _network_and_scores()

        with_numpy = active_module_detection(network, scores, n_permutations=50, random_state=np.int64(7))
        with_int = active_module_detection(network, scores, n_permutations=50, random_state=7)

        assert with_numpy == with_int

    @pytest.mark.parametrize("bad", [1.5, "42", [42], b"42"])
    def test_a_seed_that_is_not_an_integer_is_a_type_error(self, bad: Any) -> None:
        network, scores = _network_and_scores()

        with pytest.raises(TypeError, match=r"random_state must be an int or None"):
            active_module_detection(network, scores, random_state=bad)


class TestActiveModuleArgumentValidation:
    @pytest.mark.parametrize("bad", [0, -1, -1000, 2.5, "10", None, float("nan")])
    def test_n_permutations_must_be_a_positive_integer(self, bad: Any) -> None:
        network, scores = _network_and_scores()

        with pytest.raises(ValueError, match=r"n_permutations must be an integer >= 1"):
            active_module_detection(network, scores, n_permutations=bad)

    def test_zero_permutations_no_longer_report_a_p_value_of_one(self) -> None:
        network, scores = _network_and_scores()

        with pytest.raises(ValueError):
            active_module_detection(network, scores, n_permutations=0)

    @pytest.mark.parametrize("bad", [0.0, -0.05, 1.5, float("nan"), float("inf")])
    def test_alpha_must_lie_in_zero_one(self, bad: float) -> None:
        network, scores = _network_and_scores()

        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\]"):
            active_module_detection(network, scores, alpha=bad)

    def test_alpha_one_is_allowed(self) -> None:
        network, scores = _network_and_scores()

        modules = active_module_detection(network, scores, alpha=1.0, n_permutations=20)

        assert isinstance(modules, list)

    @pytest.mark.parametrize("bad", [-0.001, 1.2, float("nan"), float("inf")])
    def test_a_node_p_value_outside_zero_one_names_the_node(self, bad: float) -> None:
        network, scores = _network_and_scores()
        scores["B1"] = bad

        with pytest.raises(ValueError, match=r"scores must be p-values in \[0, 1\]; 'B1' has"):
            active_module_detection(network, scores)

    def test_a_negative_p_value_is_not_scored_as_the_most_significant_node(self) -> None:
        network = {"X": ["Y"], "Y": ["X"]}

        with pytest.raises(ValueError):
            active_module_detection(network, {"X": -0.5, "Y": 0.9}, n_permutations=10)

    def test_the_p_values_of_every_omic_are_checked_too(self) -> None:
        network, scores = _network_and_scores()

        with pytest.raises(ValueError, match=r"omic_scores\['rna'\] must be p-values in \[0, 1\]; 'A0' has"):
            active_module_detection(network, scores, omic_scores={"rna": {"A0": 3.0}})

    def test_a_p_value_of_exactly_zero_or_one_is_valid(self) -> None:
        network = {"X": ["Y"], "Y": ["X"], "Z": []}

        modules = active_module_detection(network, {"X": 0.0, "Y": 0.001, "Z": 1.0}, n_permutations=20)

        assert [m["module_genes"] for m in modules] == [["X", "Y"]]


class TestPathwayTopologyPermutations:
    def test_the_default_reproduces_the_fixed_seed_result(self) -> None:
        result = pathway_topology_analysis(GRAPH, GENE_SCORES)

        assert result["p_value"] == 54 / 1001  # 53 of 1000 permutations scored as high, computed before the change
        assert result["impact_factor"] == pytest.approx(1.6982506667511246, rel=1e-12)
        assert result["perturbed_genes"] == ["S", "T1", "U"]

    @pytest.mark.parametrize(("seed", "expected"), [(0, 50), (1, 56), (7, 51), (2024, 57)])
    def test_another_seed_changes_only_the_permutation_p_value(self, seed: int, expected: int) -> None:
        reference = pathway_topology_analysis(GRAPH, GENE_SCORES)

        result = pathway_topology_analysis(GRAPH, GENE_SCORES, random_state=seed)

        assert result["p_value"] == expected / 1001
        assert result["impact_factor"] == reference["impact_factor"]
        assert result["pathway_perturbation"] == reference["pathway_perturbation"]

    def test_the_number_of_permutations_is_an_argument(self) -> None:
        result = pathway_topology_analysis(GRAPH, GENE_SCORES, n_permutations=200)

        assert result["p_value"] == 9 / 201

    @pytest.mark.parametrize("bad", [0, -5, 1.5, None, "100"])
    def test_zero_or_fewer_permutations_raise_instead_of_reporting_a_p_value(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"n_permutations must be an integer >= 1"):
            pathway_topology_analysis(GRAPH, GENE_SCORES, n_permutations=bad)

    @pytest.mark.parametrize("bad", [1.5, "42"])
    def test_a_seed_that_is_not_an_integer_is_a_type_error(self, bad: Any) -> None:
        with pytest.raises(TypeError, match=r"random_state must be an int or None"):
            pathway_topology_analysis(GRAPH, GENE_SCORES, random_state=bad)

    def test_none_gives_a_valid_p_value(self) -> None:
        result = pathway_topology_analysis(GRAPH, GENE_SCORES, n_permutations=50, random_state=None)

        assert 1 / 51 <= result["p_value"] <= 1.0
