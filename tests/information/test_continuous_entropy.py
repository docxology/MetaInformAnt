"""Value-pinned tests for continuous (differential) information measures.

Covers the correctness fixes of the information-metrics slice:
- Kozachenko-Leonenko (KSG-style) k-NN entropy: ``psi(n) - psi(k) + log(c_d)
  + (d/n) * sum log(eps_i)`` in the true joint space (nats),
- histogram differential entropy ``-sum_i p_i log(p_i / V)`` (not the discrete
  entropy of bin masses),
- continuous MI / conditional entropy / transfer entropy operating in the true
  joint space (inputs are never flattened to 1D).

All estimators in this module return nats. Real implementations with small
deterministic data per tests/REAL_IMPLEMENTATION_TESTING_POLICY.md.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.information.metrics.core.continuous import (
    _differential_entropy_histogram,
    _differential_entropy_knn,
    conditional_entropy_continuous,
    conditional_entropy_continuous_3d,
    differential_entropy,
    mutual_information_continuous,
    transfer_entropy_continuous,
)

TRUE_NORMAL_NATS = 0.5 * math.log(2 * math.pi * math.e)  # 1.4189385332046727
TRUE_2D_INDEPENDENT_NATS = math.log(2 * math.pi * math.e)  # 2.8378770664093453


class TestHistogramDifferentialEntropy:
    """The histogram estimator must return -sum p log(p/V), not the discrete
    entropy of the bin masses."""

    def test_normal_matches_analytic_value(self) -> None:
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 10000)
        h = differential_entropy(samples, method="histogram", bins=64)
        # previously this returned the discrete entropy of bin masses
        # (no log-volume term), which is ~0.55 nats too low here
        assert h == pytest.approx(TRUE_NORMAL_NATS, abs=0.02)

    def test_normal_default_sturges_bins(self) -> None:
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 10000)
        h = differential_entropy(samples, method="histogram")
        assert h == pytest.approx(TRUE_NORMAL_NATS, abs=0.03)

    def test_uniform_support_entropy_near_zero(self) -> None:
        # H(Uniform[0,1]) = 0 nats exactly.
        rng = np.random.default_rng(7)
        samples = rng.uniform(0.0, 1.0, 10000)
        h = differential_entropy(samples, method="histogram", bins=50)
        assert abs(h) < 0.02

    def test_equals_discrete_mass_entropy_plus_log_volume(self) -> None:
        # Exact algebraic pin: H_hist = -sum p log p + log(bin width).
        rng = np.random.default_rng(3)
        samples = rng.normal(0, 1, 500)
        bins = 20
        counts, edges = np.histogram(samples, bins=bins)
        masses = counts[counts > 0] / samples.size
        width = edges[1] - edges[0]
        expected = -np.sum(masses * np.log(masses)) + math.log(width)
        h = differential_entropy(samples, method="histogram", bins=bins)
        assert h == pytest.approx(float(expected), abs=1e-12)

    def test_helper_rejects_multivariate_input(self) -> None:
        with pytest.raises(ValueError, match="1D"):
            _differential_entropy_histogram(np.ones((50, 2)))


class TestKozachenkoLeonenkoEntropy:
    """KL/KSG k-NN entropy: psi(n) - psi(k) + log(c_d) + (d/n) sum log(eps_i)."""

    def test_normal_matches_analytic_value(self) -> None:
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 5000)
        h = differential_entropy(samples, method="knn")
        # previously implemented mean(log eps) + log(2) + gamma with the k-th
        # NN distance, which is not the KL/KSG estimator
        assert h == pytest.approx(TRUE_NORMAL_NATS, abs=0.05)

    def test_uniform_support_entropy_near_zero(self) -> None:
        # Pins the log(c_d) constant: with eps taken as the plain k-NN
        # distance, H(uniform[0,1]) -> 0.
        rng = np.random.default_rng(7)
        samples = rng.uniform(0.0, 1.0, 2000)
        h = differential_entropy(samples, method="knn")
        # KL estimator sd here is ~0.01; 0.05 still pins the c_1 = 2 convention
        # (a twice-the-distance or log(2) convention would give ~0.69).
        assert abs(h) < 0.05

    def test_exponential_known_entropy(self) -> None:
        # H(Exp(1)) = 1 nat exactly; the bounded-support boundary biases the
        # k-NN estimator slightly upward.
        rng = np.random.default_rng(11)
        samples = rng.exponential(1.0, 3000)
        h = differential_entropy(samples, method="knn")
        assert h == pytest.approx(1.0, abs=0.1)

    def test_joint_two_dimensional_entropy(self) -> None:
        # Independent bivariate normal: H = ln(2 pi e); the dimension handling
        # (c_2 = pi, d = 2 in the mean-log-distance term) must be explicit.
        rng = np.random.default_rng(5)
        joint = rng.normal(0, 1, (2000, 2))
        h = differential_entropy(joint, method="knn")
        assert h == pytest.approx(TRUE_2D_INDEPENDENT_NATS, abs=0.08)

    def test_helper_rejects_multivariate_input(self) -> None:
        with pytest.raises(ValueError, match="1D"):
            _differential_entropy_knn(np.ones((50, 2)))

    def test_duplicate_samples_use_deterministic_jitter(self) -> None:
        # All-duplicate data is degenerate (true entropy -inf); the estimator
        # must not raise or return nan/inf and must be reproducible.
        samples = np.full(50, 3.14)
        h1 = differential_entropy(samples, method="knn")
        h2 = differential_entropy(samples, method="knn")
        assert math.isfinite(h1)
        assert h1 == h2


class TestJointSpaceEstimation:
    """MI / conditional entropy / transfer entropy must work in the true joint
    space (never flattened to (-1, 1))."""

    def test_histogram_joint_entropy_of_independent_normals(self) -> None:
        rng = np.random.default_rng(9)
        joint = rng.normal(0, 1, (5000, 2))
        h = differential_entropy(joint, method="histogram", bins=32)
        assert h == pytest.approx(TRUE_2D_INDEPENDENT_NATS, abs=0.15)

    def test_joint_knn_not_flattened(self) -> None:
        # H(X,Y) of independent normals must be ~H(X) + H(Y). The previous
        # flattening behaviour computed the entropy of the 2n concatenated
        # values, which is close to H(X) alone (not 2 H(X)).
        rng = np.random.default_rng(17)
        joint = rng.normal(0, 1, (2000, 2))
        h_joint = differential_entropy(joint, method="knn")
        h_single = differential_entropy(joint[:, 0], method="knn")
        assert h_joint == pytest.approx(2 * h_single, abs=0.15)

    def test_mutual_information_knn_correlated_gaussian(self) -> None:
        # Bivariate normal with rho = 0.9: I = -0.5 * ln(1 - rho^2) nats.
        rho = 0.9
        rng = np.random.default_rng(23)
        joint = rng.multivariate_normal([0.0, 0.0], [[1.0, rho], [rho, 1.0]], size=2000)
        mi = mutual_information_continuous(joint[:, 0], joint[:, 1], method="knn")
        # The KL/KSG entropy-difference estimator carries a small positive
        # bias from the marginal terms at n=2000 (~+0.10 nats here).
        assert mi == pytest.approx(-0.5 * math.log(1 - rho**2), abs=0.15)

    def test_mutual_information_knn_independent_near_zero(self) -> None:
        rng = np.random.default_rng(29)
        x = rng.normal(0, 1, 2000)
        y = rng.normal(0, 1, 2000)
        mi = mutual_information_continuous(x, y, method="knn")
        assert 0.0 <= mi < 0.15

    def test_mutual_information_histogram_correlated_positive(self) -> None:
        rng = np.random.default_rng(31)
        x = rng.normal(0, 1, 500)
        y = x + rng.normal(0, 0.1, 500)
        assert mutual_information_continuous(x, y) > 0.0

    def test_conditional_entropy_independent_near_marginal(self) -> None:
        # H(X|Y) = H(X) for independent variables (nats).
        rng = np.random.default_rng(37)
        x = rng.normal(0, 1, 300)
        y = rng.normal(0, 1, 300)
        h = conditional_entropy_continuous(x, y, method="knn")
        assert h == pytest.approx(TRUE_NORMAL_NATS, abs=0.2)

    def test_conditional_entropy_3d_deterministic_dependency(self) -> None:
        # y = 2x exactly and z independent: H(X|Y,Z) must collapse to ~0,
        # which requires the (x, y, z) joint entropy to be estimated in the
        # true 3D joint space.
        rng = np.random.default_rng(41)
        x = rng.normal(0, 1, 400)
        y = 2.0 * x
        z = rng.normal(0, 1, 400)
        h = conditional_entropy_continuous_3d(x, y, z, method="knn")
        assert 0.0 <= h < 0.05

    def test_transfer_entropy_knn_coupling_dominates_independent(self) -> None:
        rng = np.random.default_rng(43)
        x = rng.normal(0, 1, 3000)
        coupled = np.empty(3000)
        coupled[0] = rng.normal()
        for t in range(1, 3000):
            coupled[t] = x[t - 1] + 0.5 * rng.normal()
        independent = 0.5 * rng.normal(size=3000)
        te_coupled = transfer_entropy_continuous(x, coupled, lag=1, method="knn")
        te_independent = transfer_entropy_continuous(
            x, independent, lag=1, method="knn"
        )
        # True TE for y_{t} = x_{t-1} + N(0, 0.25) is 0.5*ln(5) ~ 0.805 nats;
        # independent series must give ~0 (clipped noise).
        assert te_coupled > 0.3
        assert te_coupled > te_independent
        assert te_independent < 0.2
