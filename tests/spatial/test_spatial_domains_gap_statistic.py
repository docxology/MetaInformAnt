"""``spatial_domains`` chooses the number of domains with the gap statistic (Tibshirani et al., 2001).

The docstring promised a gap statistic but the code ran an elbow heuristic on the second difference of the
inertia curve whose index arithmetic is off by one: on well-separated planted blobs it returned 8, 2, 3, 4
and 2 domains for 2, 3, 4, 5 and 6 true ones. Real k-means and real random draws only.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from metainformant.spatial.analysis.clustering import spatial_domains


def _planted_domains(
    k: int, per_domain: int = 40, n_genes: int = 12, seed: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``k`` well-separated expression blobs, each in its own patch of the tissue; returns X, coords, truth.

    Forty spots per domain is about where the gap statistic becomes reliable on such data (with thirty it can
    keep preferring a few extra clusters, as the gap statistic is known to do for small samples).
    """
    rng = np.random.default_rng(seed)
    centres = rng.normal(0.0, 4.0, size=(k, n_genes))
    expression = np.vstack([rng.normal(centre, 0.5, size=(per_domain, n_genes)) for centre in centres])
    angles = np.linspace(0.0, 2.0 * np.pi, k, endpoint=False)
    coords = np.vstack([rng.normal([10.0 * np.cos(a), 10.0 * np.sin(a)], 1.0, size=(per_domain, 2)) for a in angles])
    truth = np.repeat(np.arange(k), per_domain)
    return np.abs(expression) + 0.1, coords, truth


class TestAutomaticNumberOfDomains:
    @pytest.mark.parametrize("n_true", [2, 3, 4, 5, 6])
    def test_the_planted_number_of_domains_is_recovered(self, n_true: int) -> None:
        expression, coords, truth = _planted_domains(n_true)

        result = spatial_domains(expression, coords, n_domains=None, n_pcs=8, seed=1)

        assert result.n_clusters == n_true
        assert adjusted_rand_score(truth, result.labels) == pytest.approx(1.0)
        assert result.metadata["n_domains_selection"] == "gap_statistic"

    def test_a_given_number_is_used_and_recorded(self) -> None:
        expression, coords, _ = _planted_domains(3)

        result = spatial_domains(expression, coords, n_domains=4, n_pcs=8, seed=1)

        assert result.n_clusters == 4
        assert result.metadata["n_domains_selection"] == "user"

    def test_the_choice_is_reproducible_and_leaves_the_global_random_state_alone(self) -> None:
        expression, coords, _ = _planted_domains(4)

        np.random.seed(99)
        before = np.random.get_state()
        first = spatial_domains(expression, coords, n_domains=None, n_pcs=8, seed=3)
        after = np.random.get_state()
        again = spatial_domains(expression, coords, n_domains=None, n_pcs=8, seed=3)

        assert before[0] == after[0] and np.array_equal(before[1], after[1]) and before[2:] == after[2:]
        np.testing.assert_array_equal(first.labels, again.labels)

    def test_data_without_domain_structure_gives_the_minimum_of_two(self) -> None:
        from metainformant.spatial.analysis.clustering import _estimate_n_clusters

        rng = np.random.default_rng(5)
        noise = rng.random((120, 6))
        coords = rng.random((120, 2))

        assert _estimate_n_clusters(noise, coords, max_k=10, seed=1) == 2

    def test_a_tiny_dataset_falls_back_to_two(self) -> None:
        from metainformant.spatial.analysis.clustering import _estimate_n_clusters

        assert _estimate_n_clusters(np.array([[0.0], [1.0]]), np.array([[0.0, 0.0], [1.0, 1.0]])) == 2


class TestGapStatisticDefinition:
    def test_gap_and_spread_follow_the_published_formulas(self) -> None:
        from metainformant.spatial.analysis.clustering import _gap_and_spread

        rng = np.random.default_rng(1)
        observed = np.vstack([rng.normal(0.0, 0.3, (25, 3)), rng.normal(4.0, 0.3, (25, 3))])
        references = [rng.uniform(observed.min(axis=0), observed.max(axis=0), observed.shape) for _ in range(6)]
        k = 2

        gap, spread = _gap_and_spread(observed, references, k, seed=0)

        def log_w(data: np.ndarray) -> float:
            return float(np.log(KMeans(n_clusters=k, random_state=0, n_init=5).fit(data).inertia_))

        reference_logs = np.array([log_w(reference) for reference in references])
        assert gap == pytest.approx(reference_logs.mean() - log_w(observed))
        assert spread == pytest.approx(reference_logs.std() * np.sqrt(1.0 + 1.0 / 6.0))
        assert gap > 0.5  # two clean clusters are far more compact than a uniform cloud

    @pytest.mark.parametrize(
        "gaps, spreads, expected",
        [
            ([0.2, 0.6, 0.9, 0.8, 0.85], [0.05] * 5, 3),  # Gap(3) = 0.9 >= Gap(4) - s_4 = 0.75
            ([0.2, 0.6, 0.9, 1.3], [0.05] * 4, None),  # still rising: no k can be judged
            ([0.1, 0.12, 0.2], [0.3, 0.3, 0.3], 1),  # within one standard error already at k = 1
            ([0.5, 0.1, 0.2], [0.01, 0.01, 0.01], 1),  # a falling gap satisfies the rule at once
            ([0.4], [0.1], None),  # a single value cannot be compared with a next one
        ],
    )
    def test_selection_rule(self, gaps: list[float], spreads: list[float], expected: int | None) -> None:
        from metainformant.spatial.analysis.clustering import _first_k_by_gap_rule

        assert _first_k_by_gap_rule(gaps, spreads) == expected
