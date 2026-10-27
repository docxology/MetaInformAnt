"""Regression tests: permutation tests and bootstraps in ``information`` take a seed and leave the global state alone.

* ``transfer_entropy``, ``information_flow_network`` and ``mutual_information_network``
  (``network_info.information_flow``) shuffled with the *global* ``random`` module. Their p-values, and therefore
  which edges a network contained, changed on every call, could not be pinned with a seed, and every call consumed the
  caller's global random stream.
* ``entropy_bootstrap_confidence`` (``metrics.core.estimation``) called ``np.random.seed(random_state)``: with a seed
  it clobbered the caller's global NumPy state, and with ``random_state=None`` it *re*seeded it from system entropy.

All four now take a seed (an int, or a ``random.Random`` for the pure-Python permutation tests): the same seed
reproduces the output, different seeds differ, and the global ``random`` / ``numpy.random`` states are identical
before and after. Quantities that are not random (the TE value, the MI matrix) are checked against independent
computations (``sklearn.metrics.mutual_info_score`` on separately binned data). No test doubles.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest
from sklearn.metrics import mutual_info_score

from metainformant.information.metrics.core.estimation import entropy_bootstrap_confidence
from metainformant.information.network_info.information_flow import (
    information_flow_network,
    mutual_information_network,
    transfer_entropy,
)


class GlobalStateGuard:
    """Seeds both global generators to known states and checks that a block of code left them exactly there."""

    def __enter__(self) -> "GlobalStateGuard":
        random.seed(20261027)
        np.random.seed(20261027)
        self.python_state = random.getstate()
        self.numpy_state = np.random.get_state()
        return self

    def __exit__(self, *exc_info) -> None:
        if exc_info[0] is not None:
            return
        assert random.getstate() == self.python_state, "the global `random` state changed"
        after = np.random.get_state()
        assert after[0] == self.numpy_state[0] and after[2:] == self.numpy_state[2:]
        assert np.array_equal(after[1], self.numpy_state[1]), "the global NumPy random state changed"


def _series(n: int = 120, seed: int = 1) -> tuple[list[float], list[float], list[float]]:
    """x, y driven by x one step earlier, and z independent of both (all from an independent generator)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    y = np.empty(n)
    y[0] = rng.normal()
    y[1:] = 0.8 * x[:-1] + 0.3 * rng.normal(size=n - 1)
    z = rng.normal(size=n)
    return x.tolist(), y.tolist(), z.tolist()


X, Y, Z = _series()
W = (np.asarray(X) + 0.3 * np.random.default_rng(5).normal(size=len(X))).tolist()  # moves with x at the same time


# ---------------------------------------------------------------------------
# transfer_entropy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("random_state", [3, None])
def test_transfer_entropy_leaves_the_global_random_state_alone(random_state: int | None) -> None:
    with GlobalStateGuard():
        transfer_entropy(X, Y, lag=1, n_bins=4, random_state=random_state)


def test_transfer_entropy_is_reproducible_for_a_seed_and_the_te_value_does_not_depend_on_it() -> None:
    first = transfer_entropy(X, Z, lag=1, n_bins=4, random_state=3)
    second = transfer_entropy(X, Z, lag=1, n_bins=4, random_state=3)
    others = [transfer_entropy(X, Z, lag=1, n_bins=4, random_state=seed) for seed in range(4, 10)]

    assert first == second
    assert all(other["te"] == first["te"] for other in others)
    assert len({other["p_value"] for other in others} | {first["p_value"]}) > 1, "the permutation p-value must vary"


def test_transfer_entropy_accepts_a_generator_instance() -> None:
    assert transfer_entropy(X, Z, n_bins=4, random_state=random.Random(5)) == transfer_entropy(
        X, Z, n_bins=4, random_state=5
    )


def test_a_driven_series_is_significant_and_an_independent_one_is_not() -> None:
    driven = transfer_entropy(X, Y, lag=1, n_bins=4, random_state=1)
    independent = transfer_entropy(Z, Y, lag=1, n_bins=4, random_state=1)

    assert driven["p_value"] == pytest.approx(1 / 101)  # no permutation reached the observed TE
    assert independent["p_value"] > 0.05
    assert driven["te"] > independent["te"]


# ---------------------------------------------------------------------------
# information_flow_network / mutual_information_network
# ---------------------------------------------------------------------------


def test_information_flow_network_is_reproducible_and_leaves_global_state_alone() -> None:
    series = {"x": X[:60], "y": Y[:60], "z": Z[:60]}
    with GlobalStateGuard():
        first = information_flow_network(series, method="transfer_entropy", threshold=0.5, random_state=2)
    second = information_flow_network(series, method="transfer_entropy", threshold=0.5, random_state=2)
    others = [
        information_flow_network(series, method="transfer_entropy", threshold=0.5, random_state=seed)
        for seed in range(3, 8)
    ]

    assert first == second
    assert any(other["edges"] != first["edges"] for other in others), "p-values (hence edges) must depend on the seed"


def _binned(values: list[float], n_bins: int) -> np.ndarray:
    array = np.asarray(values)
    width = (array.max() - array.min()) / n_bins
    return np.minimum(((array - array.min()) / width).astype(int), n_bins - 1)


def test_mutual_information_network_is_reproducible_and_matches_sklearn_mutual_information() -> None:
    data = np.column_stack([X, W, Z])
    names = ["x", "w", "z"]
    with GlobalStateGuard():
        first = mutual_information_network(data, names, n_bins=4, threshold=0.03, random_state=5)
    second = mutual_information_network(data, names, n_bins=4, threshold=0.03, random_state=5)
    other = mutual_information_network(data, names, n_bins=4, threshold=0.03, random_state=6)

    assert first == second
    assert other["mi_matrix"] == first["mi_matrix"], "the MI matrix is not random"
    columns = [_binned(column, 4) for column in (X, W, Z)]
    for i in range(3):
        for j in range(i + 1, 3):
            expected_bits = mutual_info_score(columns[i], columns[j]) / math.log(2)
            assert first["mi_matrix"][i][j] == pytest.approx(expected_bits, rel=1e-9, abs=1e-12)
    significant = {frozenset(pair) for pair in first["significant_pairs"]}
    assert frozenset(("x", "w")) in significant, "the coupled pair must be found"
    assert frozenset(("x", "z")) not in significant, "the independent pair must not be"


def test_mutual_information_network_p_values_depend_on_the_seed() -> None:
    data = np.column_stack([X, Z])  # independent columns: p-values are spread out
    p_values = {
        mutual_information_network(data, ["x", "z"], n_bins=4, threshold=1.1, random_state=seed)["edges"][0]["p_value"]
        for seed in range(8)
    }

    assert len(p_values) > 1


# ---------------------------------------------------------------------------
# entropy_bootstrap_confidence
# ---------------------------------------------------------------------------

COUNTS = {"A": 30, "B": 20, "C": 10, "D": 5}


@pytest.mark.parametrize("random_state", [7, None])
def test_entropy_bootstrap_confidence_leaves_the_global_numpy_state_alone(random_state: int | None) -> None:
    with GlobalStateGuard():
        entropy_bootstrap_confidence(COUNTS, n_bootstraps=40, random_state=random_state)


def test_entropy_bootstrap_confidence_is_reproducible_and_seed_dependent() -> None:
    first = entropy_bootstrap_confidence(COUNTS, n_bootstraps=60, random_state=7)
    second = entropy_bootstrap_confidence(COUNTS, n_bootstraps=60, random_state=7)
    other = entropy_bootstrap_confidence(COUNTS, n_bootstraps=60, random_state=8)

    assert first == second
    assert (first["ci_lower"], first["ci_upper"]) != (other["ci_lower"], other["ci_upper"])
    assert first["entropy"] == other["entropy"], "the point estimate does not depend on the seed"
    assert first["ci_lower"] <= first["ci_upper"]
