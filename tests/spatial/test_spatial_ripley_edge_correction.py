"""Ripley's K applies the documented isotropic edge correction for every number of points.

The estimator used to apply no correction at all for 100 points or fewer, and above that a capped ad hoc
weight that is not Ripley's. The expected values below come from an independent brute force: every ordered
pair is weighted by 1 / (share of the circle about the first point that lies inside the window), with the
share found by sampling the circle densely. Real computation only.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.neighborhood import ripley_k


def _inside_share(x: float, y: float, radius: float, window: tuple[float, float, float, float], n_angles: int) -> float:
    """Share of the circle of ``radius`` about ``(x, y)`` inside the window, by midpoint sampling of the angle."""
    x_min, x_max, y_min, y_max = window
    angle = (np.arange(n_angles) + 0.5) * 2.0 * np.pi / n_angles
    px, py = x + radius * np.cos(angle), y + radius * np.sin(angle)
    return float(np.mean((px >= x_min) & (px <= x_max) & (py >= y_min) & (py <= y_max)))


def _brute_force_k(points: np.ndarray, radii: np.ndarray, area: float, n_angles: int = 200_000) -> np.ndarray:
    """K(r) = area / (n (n - 1)) * sum over ordered pairs within r of 1 / inside-share (window = bounding box)."""
    window = (points[:, 0].min(), points[:, 0].max(), points[:, 1].min(), points[:, 1].max())
    n = len(points)
    total = np.zeros(len(radii))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            distance = float(np.hypot(*(points[i] - points[j])))
            if distance <= radii.max():
                weight = 1.0 / _inside_share(points[i, 0], points[i, 1], distance, window, n_angles)
                total += weight * (distance <= radii)
    return area * total / (n * (n - 1))


def _boundary_heavy_points(n: int, seed: int) -> np.ndarray:
    points = np.random.default_rng(seed).random((n, 2))
    points[:4] = [[0.02, 0.5], [0.5, 0.03], [0.97, 0.9], [0.5, 0.5]]  # some close to the window boundary
    return points


class TestEdgeCorrectionAtEveryN:
    def test_a_small_sample_is_corrected_exactly(self) -> None:
        points = _boundary_heavy_points(12, seed=0)
        radii = np.array([0.1, 0.2, 0.3])
        area = float(np.ptp(points[:, 0]) * np.ptp(points[:, 1]))

        result = ripley_k(points, radii, area, n_simulations=3, seed=1)

        np.testing.assert_allclose(result.k_values, _brute_force_k(points, radii, area), rtol=1e-4)

    def test_a_small_sample_is_larger_than_the_uncorrected_estimate(self) -> None:
        points = _boundary_heavy_points(12, seed=0)
        radii = np.array([0.1, 0.2, 0.3])
        area = float(np.ptp(points[:, 0]) * np.ptp(points[:, 1]))

        corrected = ripley_k(points, radii, area, n_simulations=3, seed=1)
        uncorrected = ripley_k(points, radii, area, n_simulations=3, seed=1, edge_correction="none")

        assert np.all(corrected.k_values > uncorrected.k_values)

    def test_more_than_a_hundred_points_use_the_same_correction(self) -> None:
        points = _boundary_heavy_points(150, seed=2)
        radii = np.array([0.03, 0.06, 0.09])
        area = float(np.ptp(points[:, 0]) * np.ptp(points[:, 1]))

        result = ripley_k(points, radii, area, n_simulations=2, seed=1)

        np.testing.assert_allclose(result.k_values, _brute_force_k(points, radii, area), rtol=1e-4)

    def test_uniform_points_give_k_close_to_pi_r_squared_only_with_the_correction(self) -> None:
        points = np.random.default_rng(3).random((1500, 2))
        radii = np.array([0.1, 0.2])

        corrected = ripley_k(points, radii, 1.0, n_simulations=2, seed=1).k_values
        uncorrected = ripley_k(points, radii, 1.0, n_simulations=2, seed=1, edge_correction="none").k_values

        expected = np.pi * radii**2
        np.testing.assert_allclose(corrected, expected, rtol=0.03)
        assert np.all(uncorrected < 0.95 * expected)  # without it the boundary loses pairs: 0.92x and 0.84x here

    def test_the_csr_envelope_uses_the_same_correction(self) -> None:
        points = np.random.default_rng(4).random((200, 2))
        radii = np.array([0.1, 0.2])

        result = ripley_k(points, radii, 1.0, n_simulations=30, seed=5)

        # for uniform points the pi r^2 line sits inside the 95% envelope of the corrected K of CSR patterns
        expected_l = np.zeros(2)  # L(r) = sqrt(K / pi) - r = 0 under CSR
        assert np.all(result.csr_envelope_lower <= expected_l + 0.01)
        assert np.all(result.csr_envelope_upper >= expected_l - 0.01)


class TestEdgeCorrectionParameter:
    def test_none_counts_each_pair_once_in_both_directions(self) -> None:
        points = np.random.default_rng(6).random((40, 2))
        radii = np.array([0.1, 0.25])
        distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)

        result = ripley_k(points, radii, 1.0, n_simulations=2, seed=1, edge_correction="none")

        pairs_within = [(np.triu(distances <= r, k=1)).sum() for r in radii]
        np.testing.assert_allclose(result.k_values, [2 * count / (40 * 39) for count in pairs_within])

    def test_an_unknown_value_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="edge_correction"):
            ripley_k(
                np.random.default_rng(0).random((10, 2)),
                np.array([0.1]),
                1.0,
                edge_correction="border",  # type: ignore[arg-type]
            )

    def test_a_pattern_without_extent_has_no_window_to_correct_for(self) -> None:
        line = np.column_stack([np.linspace(0.0, 1.0, 20), np.full(20, 0.5)])

        with pytest.raises(ValueError, match="no extent"):
            ripley_k(line, np.array([0.1]), 1.0, n_simulations=2)

        # the uncorrected estimator is still defined for it
        assert ripley_k(line, np.array([0.1]), 1.0, n_simulations=2, edge_correction="none").k_values[0] > 0

    @pytest.mark.parametrize("radii", [np.array([]), np.array([0.1, -0.2])])
    def test_empty_or_negative_radii_are_rejected(self, radii: np.ndarray) -> None:
        with pytest.raises(ValueError, match="radii"):
            ripley_k(np.random.default_rng(0).random((10, 2)), radii, 1.0, n_simulations=2)


class TestWeights:
    """The weight function and the pair sums, checked directly (private helpers, imported lazily)."""

    def test_weights_equal_the_inverse_of_the_integrated_inside_share(self) -> None:
        from metainformant.spatial.analysis.neighborhood import _isotropic_edge_weights

        window = (0.0, 1.0, 0.0, 1.0)
        rng = np.random.default_rng(7)
        centres = rng.random((250, 2))
        centres[:10] = [
            [0, 0],
            [1, 1],
            [0, 1],
            [1, 0],
            [0.5, 0],
            [0, 0.5],
            [1, 0.5],
            [0.5, 1],
            [0.001, 0.001],
            [0.5, 0.5],
        ]
        radii = rng.uniform(0.0, 1.3, 250)

        weights = _isotropic_edge_weights(centres[:, 0], centres[:, 1], radii, window)

        shares = np.array([_inside_share(x, y, r, window, 100_000) for (x, y), r in zip(centres, radii)])
        usable = shares > 1e-3  # a circle almost entirely outside has a huge, numerically touchy weight
        np.testing.assert_allclose(1.0 / weights[usable], shares[usable], atol=2e-5)

    def test_a_circle_inside_the_window_has_weight_one_and_a_boundary_point_has_two(self) -> None:
        from metainformant.spatial.analysis.neighborhood import _isotropic_edge_weights

        window = (0.0, 1.0, 0.0, 1.0)

        inside = _isotropic_edge_weights(np.array([0.5]), np.array([0.5]), np.array([0.3]), window)
        on_edge = _isotropic_edge_weights(np.array([0.5]), np.array([0.0]), np.array([1e-6]), window)
        at_corner = _isotropic_edge_weights(np.array([0.0]), np.array([0.0]), np.array([1e-6]), window)

        assert inside[0] == pytest.approx(1.0)
        assert on_edge[0] == pytest.approx(2.0)  # a small circle centred on a straight edge is half inside
        assert at_corner[0] == pytest.approx(4.0)  # and on a corner a quarter inside

    def test_the_kdtree_and_the_plain_pair_search_agree(self) -> None:
        from metainformant.spatial.analysis.neighborhood import _ripley_pair_sums

        points = _boundary_heavy_points(60, seed=8)
        radii = np.array([0.05, 0.1, 0.2])
        window = (points[:, 0].min(), points[:, 0].max(), points[:, 1].min(), points[:, 1].max())

        with_tree = _ripley_pair_sums(points, radii, window, "isotropic", use_tree=True)
        without = _ripley_pair_sums(points, radii, window, "isotropic", use_tree=False)

        np.testing.assert_allclose(with_tree, without)
        assert np.all(np.diff(with_tree) >= 0)  # more pairs qualify as the radius grows
