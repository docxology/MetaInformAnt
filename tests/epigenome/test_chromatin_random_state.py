"""``learn_chromatin_states`` takes its random seed from the caller.

Defect: the k-means++ start of the EM fit drew from a ``random.Random(42)`` hard-coded inside the function, so the
model for a given matrix could not be restarted from another seed to escape a poor local optimum, and nothing in the
signature said that a seed was involved at all. ``random_state`` (default 42, so results with the default are the same
as before) now sets that seed.

Invalid arguments also failed late and obscurely: ``n_states=0`` raised ``ZeroDivisionError``, ``max_iter=0`` raised
``UnboundLocalError``, and a NaN in the signal turned every responsibility into NaN and put every bin in state 0
without any error.

The reference values for the default seed were computed with the code before ``random_state`` existed.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

import numpy as np
import pytest

from metainformant.epigenome.chromatin_state.state_learning import _kmeans_init, learn_chromatin_states


def _matrix(n_bins: int, n_marks: int, n_clusters: int, seed: int) -> List[List[float]]:
    """Blocks of ten bins drawn around ``n_clusters`` centres, so there is structure to learn."""
    rng = random.Random(seed)
    centres = [[rng.uniform(0, 10) for _ in range(n_marks)] for _ in range(n_clusters)]
    rows = []
    for i in range(n_bins):
        centre = centres[(i // 10) % n_clusters]
        rows.append([max(0.0, rng.gauss(centre[j], 0.8)) for j in range(n_marks)])
    return rows


# (n_bins, n_marks, n_clusters, data seed, n_states) -> (log-likelihood, EM iterations, first twelve assignments)
_REFERENCE: Dict[Tuple[int, int, int, int, int], Tuple[float, int, List[int]]] = {
    (120, 4, 3, 7, 3): (-659.1892, 4, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2]),
    (150, 5, 4, 11, 6): (-1055.2314, 69, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1]),
    (90, 3, 3, 5, 5): (-382.2703, 43, [4, 1, 1, 3, 4, 1, 1, 3, 1, 3, 2, 2]),
}


class TestTheDefaultSeedIsTheSeedItAlwaysHad:
    @pytest.mark.parametrize("shape", sorted(_REFERENCE))
    def test_results_with_the_default_are_those_of_the_fixed_seed(self, shape: Tuple[int, int, int, int, int]) -> None:
        n_bins, n_marks, n_clusters, seed, n_states = shape
        log_likelihood, iterations, assignments = _REFERENCE[shape]

        result = learn_chromatin_states(_matrix(n_bins, n_marks, n_clusters, seed), n_states=n_states)

        # The last digits of the likelihood and the stopping iteration can move with the platform's libm; the model
        # (assignments) and its quality cannot
        assert result["log_likelihood"] == pytest.approx(log_likelihood, abs=1e-2)
        assert abs(result["n_iterations"] - iterations) <= 2
        assert result["converged"] is True
        assert result["assignments"][:12] == assignments

    def test_an_explicit_42_is_the_default(self) -> None:
        data = _matrix(150, 5, 4, 11)

        assert learn_chromatin_states(data, n_states=6, random_state=42) == learn_chromatin_states(data, n_states=6)


class TestRandomState:
    def setup_method(self) -> None:
        self.data = _matrix(150, 5, 4, 11)

    def test_the_same_seed_gives_the_same_model(self) -> None:
        first = learn_chromatin_states(self.data, n_states=6, random_state=3)
        second = learn_chromatin_states(self.data, n_states=6, random_state=3)

        assert first == second

    def test_another_seed_can_reach_another_local_optimum(self) -> None:
        # Six states for four underlying groups: the starting centroids decide where EM ends up
        likelihoods = {
            seed: learn_chromatin_states(self.data, n_states=6, random_state=seed)["log_likelihood"]
            for seed in range(6)
        }

        assert len(set(likelihoods.values())) >= 3
        assert max(likelihoods.values()) - min(likelihoods.values()) > 10.0

    def test_the_seed_only_changes_the_start(self) -> None:
        # Every seed still returns a complete, valid model
        for seed in (0, 1, 2):
            result = learn_chromatin_states(self.data, n_states=6, random_state=seed)

            assert result["states"] == 6 and len(result["assignments"]) == 150
            assert set(result["assignments"]) <= set(range(6))
            assert sum(result["weights"]) == pytest.approx(1.0, abs=1e-4)

    def test_none_seeds_from_the_system_and_still_returns_a_model(self) -> None:
        result = learn_chromatin_states(self.data, n_states=3, random_state=None)

        assert result["states"] == 3 and len(result["assignments"]) == 150

    def test_numpy_integers_are_seeds(self) -> None:
        assert learn_chromatin_states(self.data, n_states=6, random_state=np.int64(5)) == learn_chromatin_states(
            self.data, n_states=6, random_state=5
        )

    @pytest.mark.parametrize("bad", ["42", 1.5, [1], object()])
    def test_a_seed_that_is_not_an_integer_is_rejected(self, bad: object) -> None:
        with pytest.raises(ValueError, match="random_state"):
            learn_chromatin_states(self.data, n_states=3, random_state=bad)  # type: ignore[arg-type]

    def test_the_global_random_state_is_neither_read_nor_changed(self) -> None:
        saved = random.getstate()  # put back afterwards: other tests may share the global generator
        try:
            random.seed(1234)
            before = random.getstate()

            learn_chromatin_states(self.data, n_states=3, random_state=8)
            learn_chromatin_states(self.data, n_states=3)
            learn_chromatin_states(self.data, n_states=3, random_state=None)

            assert random.getstate() == before
        finally:
            random.setstate(saved)

    def test_the_starting_centroids_follow_the_seed(self) -> None:
        points = [list(row) for row in self.data]

        same = [_kmeans_init(points, 6, random_state=9) for _ in range(2)]
        other = _kmeans_init(points, 6, random_state=10)

        assert same[0] == same[1]
        assert other != same[0]


class TestInvalidArguments:
    def setup_method(self) -> None:
        self.data = _matrix(60, 3, 3, 1)

    @pytest.mark.parametrize("n_states", [0, -2])
    def test_fewer_than_one_state_is_a_value_error(self, n_states: int) -> None:
        with pytest.raises(ValueError, match="n_states must be at least 1"):
            learn_chromatin_states(self.data, n_states=n_states)

    def test_no_iterations_is_a_value_error(self) -> None:
        with pytest.raises(ValueError, match="max_iter must be at least 1"):
            learn_chromatin_states(self.data, n_states=3, max_iter=0)

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
    def test_a_signal_that_is_not_finite_is_refused_not_assigned_to_state_zero(self, bad: float) -> None:
        data = [list(row) for row in self.data]
        data[17][1] = bad

        with pytest.raises(ValueError, match="finite"):
            learn_chromatin_states(data, n_states=3)
        with pytest.raises(ValueError, match="finite"):
            learn_chromatin_states(np.array(data), n_states=3)

    def test_an_empty_matrix_is_a_value_error(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            learn_chromatin_states(np.zeros((0, 3)), n_states=2)
