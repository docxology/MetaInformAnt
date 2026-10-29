"""The logarithm ``base`` of the syntactic measures is honoured throughout, and options are not silently ignored.

* ``mutual_information``, ``conditional_entropy``, ``conditional_mutual_information`` and ``transfer_entropy`` took a
  ``base`` but subtracted entropies of two units: the marginal terms came from ``shannon_entropy_from_counts``, which
  was always in bits, and the joint term was in ``base``. In nats the mutual information of a pair of binary
  sequences came out as 0.745 instead of 0.131, a conditional entropy as 0.255 instead of 0.562, and transfer entropy
  was negative (-0.254 instead of +0.049); in base 10 conditional entropies went negative.
  ``normalized_mutual_information`` and ``information_coefficient`` divided a ``base`` mutual information by bit
  entropies (1.45 and 2.67 for ratios that lie in [0, 1]); ``total_correlation`` ignored ``base`` altogether.
* ``tsallis_entropy(q=0)`` returned the Hartley entropy ``log(support)`` (Renyi's order 0); the Tsallis value is
  ``support - 1``.
* ``kl_divergence_estimator`` was advertised "with bias correction" (``bias_correction=True`` by default) and a
  ``method`` argument, and ignored both. Unknown methods now raise, and the unavailable correction is an error rather
  than a silent no-op.
* ``rate_distortion_function`` documented ``max_rate`` as the bound of the rate axis and never read it.

The oracles recompute every quantity from first principles in the test (entropies from ``collections.Counter`` and
``math.log``). No test doubles.
"""

from __future__ import annotations

import math
import random
from collections import Counter

import numpy as np
import pytest

from metainformant.information.metrics.advanced.information_projection import rate_distortion_function
from metainformant.information.metrics.core.estimation import kl_divergence_estimator
from metainformant.information.metrics.core.syntactic import (
    conditional_entropy,
    conditional_mutual_information,
    information_coefficient,
    mutual_information,
    normalized_mutual_information,
    shannon_entropy_from_counts,
    total_correlation,
    transfer_entropy,
    tsallis_entropy,
)

BASES = [2.0, math.e, 10.0]


def _sequences(seed: int = 7, length: int = 400) -> tuple[list[int], list[int], list[int]]:
    rng = random.Random(seed)
    z = [rng.randint(0, 2) for _ in range(length)]
    x = [(value + rng.randint(0, 1)) % 3 for value in z]
    y = [(a + rng.randint(0, 2)) % 3 for a in x]
    return x, y, z


def _entropy(items: list, base: float) -> float:
    counts = Counter(items)
    total = len(items)
    return -sum(c / total * math.log(c / total, base) for c in counts.values())


@pytest.mark.parametrize("base", BASES)
class TestBase:
    def test_mutual_information(self, base: float) -> None:
        x, y, _ = _sequences()
        expected = _entropy(x, base) + _entropy(y, base) - _entropy(list(zip(x, y)), base)

        assert mutual_information(x, y, base=base) == pytest.approx(expected, abs=1e-12)

    def test_conditional_entropy(self, base: float) -> None:
        x, y, _ = _sequences()
        expected = _entropy(list(zip(x, y)), base) - _entropy(y, base)

        assert conditional_entropy(x, y, base=base) == pytest.approx(expected, abs=1e-12)
        assert conditional_entropy(x, y, base=base) >= 0.0

    def test_conditional_mutual_information(self, base: float) -> None:
        x, y, z = _sequences()
        expected = (
            _entropy(list(zip(x, z)), base)
            + _entropy(list(zip(y, z)), base)
            - _entropy(z, base)
            - _entropy(list(zip(x, y, z)), base)
        )

        assert conditional_mutual_information(x, y, z, base=base) == pytest.approx(expected, abs=1e-12)
        assert conditional_mutual_information(x, y, z, base=base) >= -1e-12

    def test_total_correlation(self, base: float) -> None:
        x, y, z = _sequences()
        expected = _entropy(x, base) + _entropy(y, base) + _entropy(z, base) - _entropy(list(zip(x, y, z)), base)

        assert total_correlation([x, y, z], base=base) == pytest.approx(expected, abs=1e-12)

    def test_transfer_entropy(self, base: float) -> None:
        rng = random.Random(3)
        source = [rng.randint(0, 1) for _ in range(500)]
        target = [0] + [source[i] if rng.random() < 0.8 else rng.randint(0, 1) for i in range(499)]
        future, past, driver = target[1:], target[:-1], source[:-1]
        expected = (
            _entropy(list(zip(future, past)), base)
            + _entropy(list(zip(driver, past)), base)
            - _entropy(past, base)
            - _entropy(list(zip(future, driver, past)), base)
        )

        got = transfer_entropy(source, target, lag=1, base=base)

        assert got == pytest.approx(expected, abs=1e-12)
        assert got > 0.1 * (1 / math.log(base, 2))  # the source really drives the target

    def test_entropy_from_counts(self, base: float) -> None:
        counts = [5, 3, 2, 10]
        total = sum(counts)

        assert shannon_entropy_from_counts(counts, base=base) == pytest.approx(
            -sum(c / total * math.log(c / total, base) for c in counts), abs=1e-12
        )

    def test_the_normalised_measures_are_ratios_whatever_the_base(self, base: float) -> None:
        x, y, _ = _sequences()

        for method in ("arithmetic", "geometric", "max", "min"):
            nmi = normalized_mutual_information(x, y, method=method, base=base)
            assert nmi == pytest.approx(normalized_mutual_information(x, y, method=method, base=2.0), abs=1e-12)
            assert 0.0 <= nmi <= 1.0 + 1e-12
        ic = information_coefficient(x, y, base=base)
        assert ic == pytest.approx(information_coefficient(x, y, base=2.0), abs=1e-12)
        assert 0.0 <= ic <= 1.0 + 1e-12


def test_a_unit_change_scales_every_measure_by_the_same_factor() -> None:
    x, y, z = _sequences()
    ln2 = math.log(2)

    assert mutual_information(x, y, base=math.e) == pytest.approx(mutual_information(x, y) * ln2, rel=1e-12)
    assert conditional_entropy(x, y, base=math.e) == pytest.approx(conditional_entropy(x, y) * ln2, rel=1e-12)
    assert conditional_mutual_information(x, y, z, base=math.e) == pytest.approx(
        conditional_mutual_information(x, y, z) * ln2, rel=1e-12
    )
    assert total_correlation([x, y], base=math.e) == pytest.approx(total_correlation([x, y]) * ln2, rel=1e-12)


class TestTsallis:
    def test_order_zero_is_the_support_size_minus_one(self) -> None:
        assert tsallis_entropy([0.25] * 4, q=0.0) == pytest.approx(3.0, abs=1e-12)
        assert tsallis_entropy([0.7, 0.3, 0.0], q=0.0) == pytest.approx(
            1.0, abs=1e-12
        )  # zero-probability states drop out

    def test_general_order_and_the_base_has_no_effect(self) -> None:
        probs = [0.5, 0.3, 0.2]
        expected = (1.0 - sum(p**2 for p in probs)) / (2.0 - 1.0)

        assert tsallis_entropy(probs, q=2.0) == pytest.approx(expected, abs=1e-12)
        assert tsallis_entropy(probs, q=2.0, base=10.0) == pytest.approx(expected, abs=1e-12)


class TestKlDivergenceEstimator:
    P = [0, 0, 1, 1, 1, 2, 2, 0, 1, 2]
    Q = [0, 1, 1, 1, 2, 2, 2, 0, 1, 0]

    def test_the_default_is_the_plug_in_estimate(self) -> None:
        p, q = Counter(self.P), Counter(self.Q)
        n = len(self.P)
        expected = sum(p[k] / n * math.log2((p[k] / n) / (q[k] / n)) for k in p)

        assert kl_divergence_estimator(self.P, self.Q) == pytest.approx(expected, abs=1e-12)
        assert kl_divergence_estimator(self.P, self.Q, method="plugin", bias_correction=False) == pytest.approx(
            expected
        )

    @pytest.mark.parametrize("method", ["miller_madow", "chao_shen", "jackknife", "Plugin", ""])
    def test_an_unknown_method_raises_instead_of_running_the_plug_in(self, method: str) -> None:
        with pytest.raises(ValueError, match="method"):
            kl_divergence_estimator(self.P, self.Q, method=method)

    def test_a_bias_correction_that_does_not_exist_is_an_error_not_a_no_op(self) -> None:
        with pytest.raises(NotImplementedError, match="bias_correction"):
            kl_divergence_estimator(self.P, self.Q, bias_correction=True)

    def test_unequal_samples_still_raise(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            kl_divergence_estimator([0, 1], [0])


class TestRateDistortionMaxRate:
    SOURCE = [0.5, 0.5]
    DISTORTION = np.array([[0.0, 1.0], [1.0, 0.0]])

    def test_every_returned_point_respects_max_rate(self) -> None:
        full = rate_distortion_function(self.SOURCE, self.DISTORTION, n_points=40)
        limit = 0.3
        assert max(full["rates"]) > limit  # the default sweep does reach beyond the requested axis

        cut = rate_distortion_function(self.SOURCE, self.DISTORTION, max_rate=limit, n_points=40)

        assert cut["rates"] and max(cut["rates"]) <= limit
        assert len(cut["rates"]) < len(full["rates"])
        # the same sweep, stopped early: the surviving points are those of the full curve with R <= max_rate
        kept = [(r, d) for r, d in zip(full["rates"], full["distortions"]) if r <= limit]
        assert cut["rates"] == pytest.approx([r for r, _ in kept], abs=1e-12)
        assert cut["distortions"] == pytest.approx([d for _, d in kept], abs=1e-12)
        assert cut["critical_distortion"] == full["critical_distortion"]

    def test_the_curve_still_follows_the_closed_form(self) -> None:
        # binary symmetric source, Hamming distortion: R(D) = ln 2 - H_b(D) nats
        cut = rate_distortion_function(self.SOURCE, self.DISTORTION, max_rate=0.4, n_points=40)

        for rate, distortion in zip(cut["rates"], cut["distortions"]):
            binary_entropy = -(distortion * math.log(distortion) + (1 - distortion) * math.log(1 - distortion))
            assert rate == pytest.approx(math.log(2) - binary_entropy, abs=1e-6)

    def test_a_generous_max_rate_changes_nothing(self) -> None:
        default = rate_distortion_function(self.SOURCE, self.DISTORTION, n_points=25)
        generous = rate_distortion_function(self.SOURCE, self.DISTORTION, max_rate=50.0, n_points=25)

        assert generous == default

    @pytest.mark.parametrize("max_rate", [0.0, -1.0, float("nan")])
    def test_a_non_positive_max_rate_raises(self, max_rate: float) -> None:
        with pytest.raises(ValueError, match="max_rate"):
            rate_distortion_function(self.SOURCE, self.DISTORTION, max_rate=max_rate)
