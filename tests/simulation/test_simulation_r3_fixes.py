"""Regression tests for simulation defects fixed in the round-3 audit.

REAL IMPLEMENTATION POLICY: every test runs the real simulators and checks them against independent
oracles (closed forms and scipy); nothing is mocked or patched.

Covered here:
- ``simulate_differential_expression`` used the magnitude of a negative fold change, so a fold change
  of -4 raised the group-1 mean fourfold instead of lowering it.
- ``generate_linkage_disequilibrium_data`` used ``r_squared_target`` as a copy probability (the realised
  r-squared was its square), did not keep the requested allele frequencies, and accepted values that
  are not valid probabilities.
- ``simulate_time_series_expression`` redrew phase, amplitude and baseline for every gene at every time
  point, so successive time points were independent (lag-1 autocorrelation about 0) instead of tracing
  one sinusoid per gene.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest
from scipy import stats
from scipy.optimize import least_squares

from metainformant.core.utils import errors
from metainformant.simulation.models.popgen import generate_linkage_disequilibrium_data
from metainformant.simulation.models.rna import simulate_differential_expression, simulate_time_series_expression

# ---------------------------------------------------------------------------
# simulate_differential_expression: signed fold changes
# ---------------------------------------------------------------------------


class TestDifferentialExpressionSignedFoldChanges:
    """A fold change is a signed ratio: positive multiplies, negative divides by its magnitude."""

    # (fold change, expected group-1 / group-0 ratio of the means)
    CASES = [
        (-4.0, 0.25),
        (-2.0, 0.5),
        (4.0, 4.0),
        (2.0, 2.0),
        (0.5, 0.5),
        (0.0, 1.0),
        (-1.0, 1.0),
        (8.0, 8.0),
    ]

    @pytest.fixture(scope="class")
    def simulated(self) -> tuple[np.ndarray, np.ndarray]:
        fold_changes = np.array([fc for fc, _ in self.CASES])
        return simulate_differential_expression(6000, 20, fold_changes, rng=random.Random(2024))

    def test_group_mean_ratio_follows_the_signed_ratio(self, simulated: tuple[np.ndarray, np.ndarray]) -> None:
        """The ratio of group means matches the fold change, judged by the delta-method standard error."""
        expr, labels = simulated
        group0, group1 = expr[labels == 0], expr[labels == 1]
        for gene, (fold_change, expected_ratio) in enumerate(self.CASES):
            mean0, mean1 = group0[:, gene].mean(), group1[:, gene].mean()
            var0, var1 = group0[:, gene].var(ddof=1), group1[:, gene].var(ddof=1)
            assert mean0 > 0 and mean1 > 0, f"gene {gene} is silent, the test would learn nothing"
            # Standard error of log(mean1 / mean0) for independent samples (delta method)
            se = math.sqrt(var0 / (len(group0) * mean0**2) + var1 / (len(group1) * mean1**2))
            z = (math.log(mean1 / mean0) - math.log(expected_ratio)) / se
            assert abs(z) < 4.5, f"fold change {fold_change:+.1f}: ratio {mean1 / mean0:.3f}, z={z:.2f}"

    def test_negative_fold_change_downregulates_group_one(self, simulated: tuple[np.ndarray, np.ndarray]) -> None:
        """Welch tests agree on the direction: negative fold changes lower group 1, positive ones raise it."""
        expr, labels = simulated
        group0, group1 = expr[labels == 0], expr[labels == 1]
        for gene, (fold_change, _) in enumerate(self.CASES):
            if fold_change == 0.0 or fold_change == -1.0:
                continue
            decreases = fold_change < 1.0  # -4, -2 and the ratio 0.5 all lower the group-1 mean
            alternative = "less" if decreases else "greater"
            p_value = stats.ttest_ind(group1[:, gene], group0[:, gene], equal_var=False, alternative=alternative).pvalue
            assert p_value < 1e-6, f"fold change {fold_change}: group 1 should be {'lower' if decreases else 'higher'}"

    def test_single_negative_fold_change_lowers_the_gene(self) -> None:
        """The reported case: ``fold_changes=[-4.0]`` used to raise the gene fourfold."""
        expr, labels = simulate_differential_expression(4000, 10, np.array([-4.0]), rng=random.Random(11))
        group0, group1 = expr[labels == 0, 0], expr[labels == 1, 0]
        assert group0.mean() > 0
        assert group1.mean() < group0.mean()
        assert group1.mean() / group0.mean() == pytest.approx(0.25, rel=0.2)

    def test_genes_without_a_fold_change_are_untouched(self, simulated: tuple[np.ndarray, np.ndarray]) -> None:
        """Genes beyond the supplied fold changes keep the same mean in both groups."""
        expr, labels = simulated
        group0, group1 = expr[labels == 0], expr[labels == 1]
        n_diff = len(self.CASES)
        pvalues = [
            stats.ttest_ind(group1[:, g], group0[:, g], equal_var=False).pvalue for g in range(n_diff, expr.shape[1])
        ]
        # Twelve untouched genes: none may differ at a stringent threshold
        assert min(pvalues) > 1e-4

    @pytest.mark.parametrize("bad", [-0.5, -0.99])
    def test_negative_fold_change_below_one_fold_is_rejected(self, bad: float) -> None:
        """A fold decrease below one fold would raise the mean, contradicting its sign."""
        with pytest.raises(errors.ValidationError, match="fold decrease"):
            simulate_differential_expression(10, 5, np.array([bad]), rng=random.Random(1))

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
    def test_non_finite_fold_change_is_rejected(self, bad: float) -> None:
        """NaN used to fall through both branches and silently mean 'no change'."""
        with pytest.raises(errors.ValidationError, match="finite"):
            simulate_differential_expression(10, 5, np.array([2.0, bad]), rng=random.Random(1))


# ---------------------------------------------------------------------------
# generate_linkage_disequilibrium_data
# ---------------------------------------------------------------------------


def _r_squared(genotypes: list[list[int]], left: int, right: int) -> float:
    """Squared Pearson correlation between two columns (scipy as the oracle)."""
    a = np.array([row[left] for row in genotypes], dtype=float)
    b = np.array([row[right] for row in genotypes], dtype=float)
    return float(stats.pearsonr(a, b)[0] ** 2)


class TestLinkageDisequilibriumGenerator:
    @pytest.mark.parametrize("target", [0.1, 0.5, 0.8])
    def test_adjacent_r_squared_hits_the_target(self, target: float) -> None:
        """The realised r-squared between adjacent sites equals ``r_squared_target`` (it used to be its square)."""
        genotypes = generate_linkage_disequilibrium_data(30000, 2, r_squared_target=target, rng=random.Random(101))
        # Sampling error of r-squared is about 2 r (1 - r^2) / sqrt(n) <= 0.01 here; 0.02 is a wide margin
        assert _r_squared(genotypes, 0, 1) == pytest.approx(target, abs=0.02)

    def test_every_adjacent_pair_reaches_the_target(self) -> None:
        genotypes = generate_linkage_disequilibrium_data(8000, 25, r_squared_target=0.5, rng=random.Random(7))
        values = [_r_squared(genotypes, k, k + 1) for k in range(24)]
        assert np.mean(values) == pytest.approx(0.5, abs=0.02)
        assert min(values) > 0.4 and max(values) < 0.6

    def test_correlation_decays_geometrically_with_distance(self) -> None:
        """Equal-frequency sites follow a Markov chain: r-squared at lag k is target ** k."""
        genotypes = generate_linkage_disequilibrium_data(40000, 6, r_squared_target=0.5, rng=random.Random(5))
        for lag in (1, 2, 3):
            assert _r_squared(genotypes, 0, lag) == pytest.approx(0.5**lag, abs=0.02)

    def test_default_target_keeps_the_historic_r_squared(self) -> None:
        """With no target the generator has always produced r-squared near 0.25 (copy probability 0.5)."""
        genotypes = generate_linkage_disequilibrium_data(30000, 2, rng=random.Random(3))
        assert _r_squared(genotypes, 0, 1) == pytest.approx(0.25, abs=0.02)

    def test_requested_allele_frequencies_are_preserved(self) -> None:
        """Each site keeps its own frequency; the copy chain used to blend in the previous site's."""
        freqs = [0.3, 0.4, 0.2, 0.5, 0.1]
        n = 40000
        genotypes = generate_linkage_disequilibrium_data(n, allele_frequencies=freqs, rng=random.Random(2))
        observed = np.mean(genotypes, axis=0)
        for site, (p, seen) in enumerate(zip(freqs, observed)):
            assert abs(seen - p) < 5 * math.sqrt(p * (1 - p) / n), f"site {site}: {seen:.4f} vs {p}"

    def test_unequal_frequencies_reach_an_attainable_target(self) -> None:
        """Different frequencies can still carry the target when they allow it (D = r sqrt(p q p' q'))."""
        freqs = [0.3, 0.35, 0.4, 0.3]
        n = 60000
        genotypes = generate_linkage_disequilibrium_data(
            n, allele_frequencies=freqs, r_squared_target=0.09, rng=random.Random(8)
        )
        for k in range(3):
            assert _r_squared(genotypes, k, k + 1) == pytest.approx(0.09, abs=0.012)
        for p, seen in zip(freqs, np.mean(genotypes, axis=0)):
            assert abs(seen - p) < 5 * math.sqrt(p * (1 - p) / n)

    def test_unattainable_target_is_capped_at_the_maximum_r_squared(self) -> None:
        """A rare allele next to a common one limits r-squared; the pair is capped, not distorted."""
        p1, p2, n = 0.5, 0.1, 60000
        # Closed form: D <= min(p1 (1 - p2), (1 - p1) p2), so r-squared <= D_max^2 / (p1 q1 p2 q2)
        d_max = min(p1 * (1 - p2), (1 - p1) * p2)
        r2_max = d_max**2 / (p1 * (1 - p1) * p2 * (1 - p2))
        assert r2_max == pytest.approx(1 / 9)

        genotypes = generate_linkage_disequilibrium_data(
            n, allele_frequencies=[p1, p2], r_squared_target=0.5, rng=random.Random(4)
        )
        assert _r_squared(genotypes, 0, 1) == pytest.approx(r2_max, abs=0.012)
        for p, seen in zip((p1, p2), np.mean(genotypes, axis=0)):
            assert abs(seen - p) < 5 * math.sqrt(p * (1 - p) / n)

    def test_monomorphic_sites_are_allowed(self) -> None:
        genotypes = generate_linkage_disequilibrium_data(
            50, allele_frequencies=[0.0, 0.5, 1.0], r_squared_target=0.5, rng=random.Random(9)
        )
        assert all(row[0] == 0 and row[2] == 1 for row in genotypes)

    def test_selection_coefficient_sets_the_adjacent_correlation(self) -> None:
        """A non-zero selection coefficient replaces sqrt(target): |s| is the correlation r."""
        genotypes = generate_linkage_disequilibrium_data(
            40000, 2, selection_coefficient=0.4, r_squared_target=0.9, rng=random.Random(6)
        )
        assert math.sqrt(_r_squared(genotypes, 0, 1)) == pytest.approx(0.4, abs=0.02)

    def test_recombination_rate_is_validated_but_does_not_change_the_output(self) -> None:
        """The documented contract: the rate is accepted and checked, but the model has no recombination step."""
        kwargs = dict(n_sites=12, r_squared_target=0.5)
        low = generate_linkage_disequilibrium_data(200, recombination_rate=0.0, rng=random.Random(12), **kwargs)
        high = generate_linkage_disequilibrium_data(200, recombination_rate=0.4, rng=random.Random(12), **kwargs)
        assert low == high
        with pytest.raises(errors.ValidationError, match="recombination_rate"):
            generate_linkage_disequilibrium_data(10, 3, recombination_rate=-0.1)

    def test_output_is_reproducible_zero_one_alleles(self) -> None:
        first = generate_linkage_disequilibrium_data(30, 8, r_squared_target=0.4, rng=random.Random(21))
        second = generate_linkage_disequilibrium_data(30, 8, r_squared_target=0.4, rng=random.Random(21))
        assert first == second
        assert len(first) == 30 and all(len(row) == 8 for row in first)
        assert all(type(allele) is int and allele in (0, 1) for row in first for allele in row)

    @pytest.mark.parametrize("target", [-0.1, 1.5, float("nan")])
    def test_r_squared_target_must_be_a_proportion(self, target: float) -> None:
        with pytest.raises(errors.ValidationError, match="r_squared_target"):
            generate_linkage_disequilibrium_data(10, 3, r_squared_target=target)

    @pytest.mark.parametrize("freqs", [[0.5, 1.2], [-0.1, 0.5], [0.5, float("nan")]])
    def test_allele_frequencies_must_be_proportions(self, freqs: list[float]) -> None:
        with pytest.raises(errors.ValidationError, match="allele_frequencies"):
            generate_linkage_disequilibrium_data(10, allele_frequencies=freqs)


# ---------------------------------------------------------------------------
# simulate_time_series_expression: one coherent curve per gene
# ---------------------------------------------------------------------------

N_TIME = 400
TIMES = np.linspace(0, 4 * np.pi, N_TIME)  # the sampling times of the documented model


def _fit_rectified_sinusoid(counts: np.ndarray, frequency: float) -> tuple[float, float, float, np.ndarray]:
    """Fit ``max(0, baseline + amplitude * sin(frequency * t + phase))`` to a gene's counts.

    Returns the baseline, amplitude, phase and the Pearson residuals. When the counts are Poisson draws
    around one such curve, the Pearson dispersion ``sum(residual ** 2) / (n - 3)`` is 1; any structure the
    curve does not hold (a parameter that moves in time, a wrong shape) pushes it far above 1.
    """

    def residuals(params: np.ndarray) -> np.ndarray:
        return np.asarray(np.maximum(0.0, params[0] + params[1] * np.sin(frequency * TIMES + params[2])) - counts)

    best = None
    for phase0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        start = [max(counts.mean(), 1.0), max((counts.max() - counts.min()) / 2, 1.0), phase0]
        fit = least_squares(residuals, start, bounds=([0.0, 0.0, -np.inf], [np.inf, np.inf, np.inf]))
        if best is None or fit.cost < best.cost:
            best = fit
    baseline, amplitude, phase = best.x
    mean = np.maximum(0.0, baseline + amplitude * np.sin(frequency * TIMES + phase))
    positive = mean > 0  # where the curve is floored at zero the counts are exactly zero
    pearson = (counts[positive] - mean[positive]) / np.sqrt(mean[positive])
    return float(baseline), float(amplitude), float(phase % (2 * np.pi)), pearson


def _lag_one_autocorrelation(series: np.ndarray) -> float:
    return float(np.corrcoef(series[:-1], series[1:])[0, 1])


class TestTimeSeriesExpression:
    """Each gene is one smooth sinusoid in time; only the counting noise is drawn per observation."""

    FREQUENCIES = np.linspace(0.5, 2.0, 16)

    @pytest.fixture(scope="class")
    def series(self) -> np.ndarray:
        counts = simulate_time_series_expression(
            N_TIME, len(self.FREQUENCIES), oscillation_freq=self.FREQUENCIES, rng=random.Random(7)
        )
        return counts.astype(float)

    @pytest.fixture(scope="class")
    def fits(self, series: np.ndarray) -> list[tuple[float, float, float, np.ndarray]]:
        return [_fit_rectified_sinusoid(series[:, gene], f) for gene, f in enumerate(self.FREQUENCIES)]

    def test_successive_time_points_are_strongly_correlated(self, series: np.ndarray) -> None:
        """Phase, amplitude and baseline were redrawn at every time point, so lag-1 autocorrelation was about 0."""
        lag_one = np.array([_lag_one_autocorrelation(series[:, g]) for g in range(series.shape[1])])
        assert np.median(lag_one) > 0.9
        assert lag_one.min() > 0.8

    def test_the_correlation_is_temporal_not_a_property_of_the_values(self, series: np.ndarray) -> None:
        """Control: the same values in shuffled time order are uncorrelated at lag 1."""
        shuffled = np.random.default_rng(0).permutation(series, axis=0)
        lag_one = np.array([_lag_one_autocorrelation(shuffled[:, g]) for g in range(series.shape[1])])
        assert abs(np.median(lag_one)) < 0.1

    def test_each_gene_is_one_sinusoid_with_constant_parameters(
        self, fits: list[tuple[float, float, float, np.ndarray]]
    ) -> None:
        """A constant-parameter sinusoid leaves exactly Poisson scatter (dispersion 1); the old data left about 35."""
        dispersions = np.array([float(np.sum(f[3] ** 2) / (f[3].size - 3)) for f in fits])
        assert np.median(dispersions) == pytest.approx(1.0, abs=0.15)
        assert dispersions.max() < 2.0

    def test_fitted_parameters_lie_in_the_documented_ranges(
        self, fits: list[tuple[float, float, float, np.ndarray]]
    ) -> None:
        baselines = np.array([f[0] for f in fits])
        amplitudes = np.array([f[1] for f in fits])
        assert baselines.min() > 3 and baselines.max() < 22  # documented: uniform 5 to 20
        assert amplitudes.min() > 8 and amplitudes.max() < 105  # documented: uniform 10 to 100

    def test_parameters_belong_to_the_gene_and_differ_between_genes(
        self, fits: list[tuple[float, float, float, np.ndarray]]
    ) -> None:
        baselines = np.array([f[0] for f in fits])
        amplitudes = np.array([f[1] for f in fits])
        phases = np.array([f[2] for f in fits])
        assert np.ptp(baselines) > 5 and np.ptp(amplitudes) > 30
        # 16 phases drawn uniformly have a mean resultant length near 0.2; one shared phase would give 1
        assert abs(np.mean(np.exp(1j * phases))) < 0.85

    def test_the_noise_left_over_is_independent_between_observations(
        self, fits: list[tuple[float, float, float, np.ndarray]]
    ) -> None:
        lag_one = np.array([_lag_one_autocorrelation(f[3]) for f in fits])
        assert abs(np.median(lag_one)) < 0.1
        assert np.abs(lag_one).max() < 0.25

    def test_the_dominant_frequency_is_the_supplied_oscillation_frequency(self) -> None:
        """Frequencies 0.5, 1, 1.5 and 2 are 1, 2, 3 and 4 whole cycles over the record: the spectrum peaks there."""
        frequencies = np.array([0.5, 1.0, 1.5, 2.0])
        counts = simulate_time_series_expression(N_TIME, 4, oscillation_freq=frequencies, rng=random.Random(11))
        for gene, frequency in enumerate(frequencies):
            values = counts[:, gene].astype(float)
            spectrum = np.abs(np.fft.rfft(values - values.mean()))
            assert int(np.argmax(spectrum[1:])) + 1 == round(2 * frequency)

    def test_counts_are_non_negative_integers(self, series: np.ndarray) -> None:
        assert np.all(series >= 0) and np.all(series == np.round(series))
