"""Regression tests: the mutual-information GRN inference takes a ``random_state``.

``infer_grn(method="mutual_info")`` and ``infer_regulatory_network_from_expression(method="mutual_info")`` call
``sklearn.feature_selection.mutual_info_regression``, which adds a tiny random jitter to continuous inputs to break
ties. Neither function passed a ``random_state``, so on data with ties (counts, for example) the weights, and with a
threshold the edges, changed from run to run and consumed the global NumPy random stream; a ``random_state`` given in
``**kwargs`` was silently ignored. Both functions now take ``random_state`` (an int or a ``numpy.random.RandomState``)
and share one generator over all gene pairs.

The oracle is scikit-learn's own ``mutual_info_regression`` called by the test, with a ``RandomState`` consumed in the
order the documentation gives (regulators in order, targets in order). The data are small integer counts, which have
the ties that make the jitter matter. No test doubles.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.feature_selection import mutual_info_regression

from metainformant.networks.interaction.regulatory_analysis import infer_grn, infer_regulatory_network_from_expression

_rng = np.random.default_rng(0)
_BASE = _rng.integers(0, 4, size=40).astype(float)
COUNTS = np.column_stack(
    [
        _BASE,
        np.clip(_BASE + _rng.integers(-1, 2, size=40), 0, None),
        _rng.integers(0, 4, size=40),
        _rng.integers(0, 3, size=40),
    ]
).astype(
    float
)  # samples x genes
GENES = ["G0", "G1", "G2", "G3"]
TFS = ["G0", "G2"]


def _grn_weights(**kwargs) -> dict[tuple[str, str], float]:
    network = infer_grn(COUNTS, GENES, method="mutual_info", threshold=0.0, tf_genes=TFS, **kwargs)
    return {(u, v): data["weight"] for u, v, data in network.graph.edges(data=True)}


def _expression_weights(**kwargs) -> dict[tuple[str, str], float]:
    frame = pd.DataFrame(COUNTS.T, index=GENES)  # genes x samples
    graph = infer_regulatory_network_from_expression(frame, ["G0", "G2"], method="mutual_info", threshold=0.0, **kwargs)
    return {(u, v): data["weight"] for u, v, data in graph.edges(data=True)}


def test_infer_grn_is_reproducible_for_a_seed_and_the_seed_matters() -> None:
    first, second = _grn_weights(random_state=3), _grn_weights(random_state=3)

    assert first == second
    assert len({tuple(sorted(_grn_weights(random_state=seed).items())) for seed in range(6)}) > 1


def test_infer_grn_weights_are_scikit_learn_mutual_information_with_a_shared_generator() -> None:
    rng = np.random.RandomState(3)
    expected = {}
    for tf in TFS:
        for target in GENES:
            if target != tf:
                mi = mutual_info_regression(
                    COUNTS[:, [GENES.index(tf)]], COUNTS[:, GENES.index(target)], random_state=rng
                )
                expected[(tf, target)] = float(mi[0])

    observed = _grn_weights(random_state=3)

    assert observed == pytest.approx(expected, rel=1e-12, abs=1e-15)


def test_infer_grn_accepts_a_random_state_instance() -> None:
    assert _grn_weights(random_state=np.random.RandomState(3)) == _grn_weights(random_state=3)


def test_infer_grn_leaves_the_global_numpy_state_alone() -> None:
    np.random.seed(20261027)
    before = np.random.get_state()

    _grn_weights(random_state=3)

    after = np.random.get_state()
    assert after[0] == before[0] and after[2:] == before[2:] and np.array_equal(after[1], before[1])


def test_expression_inference_is_reproducible_for_a_seed_and_leaves_global_state_alone() -> None:
    np.random.seed(20261027)
    before = np.random.get_state()

    first = _expression_weights(random_state=5)

    after = np.random.get_state()
    assert after[0] == before[0] and after[2:] == before[2:] and np.array_equal(after[1], before[1])
    assert first == _expression_weights(random_state=5)
    assert len({tuple(sorted(_expression_weights(random_state=seed).items())) for seed in range(6)}) > 1


def test_expression_inference_weights_are_scikit_learn_mutual_information_with_a_shared_generator() -> None:
    rng = np.random.RandomState(5)
    expected = {}
    for regulator in ["G0", "G2"]:
        for target in GENES:
            if target not in ("G0", "G2"):  # the targets are the genes that are not regulators
                mi = mutual_info_regression(
                    COUNTS[:, [GENES.index(regulator)]], COUNTS[:, GENES.index(target)], random_state=rng
                )
                expected[(regulator, target)] = float(mi[0])

    assert _expression_weights(random_state=5) == pytest.approx(expected, rel=1e-12, abs=1e-15)
