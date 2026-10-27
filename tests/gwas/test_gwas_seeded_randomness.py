"""Regression tests: GWAS analysis code that samples at random uses its own seeded generator.

Defect (audit of global randomness): four places drew from the module-level ``random`` state, which made their
output different on every call and made them consume (and depend on) the state of unrelated code:

* ``summary_stats.calculate_ld_decay`` chose its variant pairs with ``random.randint``;
* ``ld_decay.ld_decay_plot`` and ``ldsr.ldsr_plot`` scattered a ``random.sample`` of their points;
* ``power.subsample_convergence`` drew its variant subsamples with ``random.sample``.

Each now takes ``seed`` (default 42) and draws from ``random.Random(seed)``. The tests check the three
properties that matter, without any test double: the same seed gives the same output whatever the global state
is, a different seed gives a different draw, and the global state is left exactly as it was found.
"""

from __future__ import annotations

import random
from typing import Any

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from metainformant.gwas.analysis.association import association_test_linear  # noqa: E402
from metainformant.gwas.analysis.ld_decay import ld_decay_plot  # noqa: E402
from metainformant.gwas.analysis.ldsr import ldsr_plot  # noqa: E402
from metainformant.gwas.analysis.power import subsample_convergence  # noqa: E402
from metainformant.gwas.analysis.summary_stats import calculate_ld_decay  # noqa: E402


def _scatter_offsets(fig: Any) -> np.ndarray:
    """The (x, y) points of the first scatter collection of a figure's first axes."""
    offsets = np.array(fig.axes[0].collections[0].get_offsets())
    plt.close(fig)
    return offsets


def _under_global_seed(global_seed: int, function: Any, /, *args: Any, **kwargs: Any) -> Any:
    """Call ``function`` after seeding the global ``random`` module (``seed`` stays free for the callee)."""
    random.seed(global_seed)
    return function(*args, **kwargs)


# ---------------------------------------------------------------------------
# summary_stats.calculate_ld_decay
# ---------------------------------------------------------------------------


def _ld_decay_inputs() -> tuple[list[list[int]], list[dict[str, Any]]]:
    """60 variants on one chromosome, 1 kb apart, in blocks of five that share a latent genotype pattern."""
    rng = np.random.default_rng(7)
    n_variants, n_samples = 60, 30
    genotypes: list[list[int]] = []
    for block in range(n_variants // 5):
        latent = rng.integers(0, 3, size=n_samples)
        for _ in range(5):
            noise = rng.random(n_samples) < 0.25
            row = np.where(noise, rng.integers(0, 3, size=n_samples), latent)
            genotypes.append([int(v) for v in row])
    variants = [{"chrom": "chr1", "pos": 1000 * (k + 1)} for k in range(n_variants)]
    return genotypes, variants


class TestCalculateLdDecay:
    def test_same_seed_gives_the_same_curve_whatever_the_global_state(self) -> None:
        genotypes, variants = _ld_decay_inputs()

        first = _under_global_seed(1, calculate_ld_decay, genotypes, variants, max_distance=60_000, bins=10, seed=5)
        second = _under_global_seed(2, calculate_ld_decay, genotypes, variants, max_distance=60_000, bins=10, seed=5)

        assert first == second
        assert first["n_pairs"] > 0

    def test_default_seed_is_reproducible(self) -> None:
        genotypes, variants = _ld_decay_inputs()

        first = _under_global_seed(11, calculate_ld_decay, genotypes, variants, max_distance=60_000, bins=10)
        second = _under_global_seed(12, calculate_ld_decay, genotypes, variants, max_distance=60_000, bins=10)

        assert first == second

    def test_different_seeds_draw_different_pairs(self) -> None:
        genotypes, variants = _ld_decay_inputs()

        one = calculate_ld_decay(genotypes, variants, max_distance=60_000, bins=10, seed=1)
        two = calculate_ld_decay(genotypes, variants, max_distance=60_000, bins=10, seed=2)

        assert one["r2_means"] != two["r2_means"]

    def test_global_random_state_is_left_untouched(self) -> None:
        genotypes, variants = _ld_decay_inputs()
        random.seed(2024)
        before = random.getstate()

        calculate_ld_decay(genotypes, variants, max_distance=60_000, bins=10)

        assert random.getstate() == before

    def test_the_curve_is_still_a_valid_r_squared_profile(self) -> None:
        """Nearby variants of one block are strongly linked: the first bin must beat the last non-empty one."""
        genotypes, variants = _ld_decay_inputs()

        result = calculate_ld_decay(genotypes, variants, max_distance=60_000, bins=6)

        means = result["r2_means"]
        assert all(0.0 <= m <= 1.0 for m in means)
        occupied = [m for m in means if m > 0.0]
        assert occupied[0] > occupied[-1]


# ---------------------------------------------------------------------------
# ld_decay.ld_decay_plot
# ---------------------------------------------------------------------------


def _ld_decay_result(n_pairs: int = 6000) -> dict[str, Any]:
    """A ``compute_ld_decay``-shaped result with more pairs than the plot displays."""
    rng = np.random.default_rng(3)
    pairs = [(int(d), float(r)) for d, r in zip(rng.integers(1, 500_000, n_pairs), rng.random(n_pairs))]
    return {
        "pairs": pairs,
        "binned": [(25_000.0, 0.5, 100), (75_000.0, 0.3, 80)],
        "r2_half_distance_bp": None,
        "n_pairs_computed": n_pairs,
    }


class TestLdDecayPlot:
    def test_same_seed_scatters_the_same_points_whatever_the_global_state(self) -> None:
        result = _ld_decay_result()

        first = _scatter_offsets(_under_global_seed(1, ld_decay_plot, result, seed=9))
        second = _scatter_offsets(_under_global_seed(2, ld_decay_plot, result, seed=9))

        assert first.shape == (5000, 2)
        np.testing.assert_array_equal(first, second)

    def test_different_seeds_scatter_different_points(self) -> None:
        result = _ld_decay_result()

        one = _scatter_offsets(ld_decay_plot(result, seed=1))
        two = _scatter_offsets(ld_decay_plot(result, seed=2))

        assert not np.array_equal(one, two)

    def test_global_random_state_is_left_untouched(self) -> None:
        result = _ld_decay_result()
        random.seed(2025)
        before = random.getstate()

        plt.close(ld_decay_plot(result))

        assert random.getstate() == before

    def test_the_displayed_points_are_a_subset_of_the_pairs(self) -> None:
        result = _ld_decay_result()
        shown = _scatter_offsets(ld_decay_plot(result))

        allowed = {(d / 1000, r) for d, r in result["pairs"]}
        assert all((float(x), float(y)) in allowed for x, y in shown)
        assert len({(float(x), float(y)) for x, y in shown}) == 5000  # sampled without replacement

    def test_few_pairs_are_all_shown(self) -> None:
        result = _ld_decay_result(n_pairs=40)

        assert _scatter_offsets(ld_decay_plot(result)).shape == (40, 2)


# ---------------------------------------------------------------------------
# ldsr.ldsr_plot
# ---------------------------------------------------------------------------


def _ldsr_inputs(n: int = 3500) -> tuple[list[float], list[float], dict[str, float]]:
    rng = np.random.default_rng(5)
    ld_scores = [float(v) for v in rng.random(n) * 10]
    chi2 = [float(1 + 0.2 * ld + rng.random()) for ld in ld_scores]
    return chi2, ld_scores, {"slope": 0.2, "intercept": 1.0, "h2_snp": 0.3, "attenuation_ratio": 0.05}


class TestLdsrPlot:
    def test_same_seed_scatters_the_same_points_whatever_the_global_state(self) -> None:
        chi2, ld, fit = _ldsr_inputs()

        first = _scatter_offsets(_under_global_seed(1, ldsr_plot, chi2, ld, fit, seed=4))
        second = _scatter_offsets(_under_global_seed(2, ldsr_plot, chi2, ld, fit, seed=4))

        assert first.shape == (3000, 2)
        np.testing.assert_array_equal(first, second)

    def test_different_seeds_scatter_different_points(self) -> None:
        chi2, ld, fit = _ldsr_inputs()

        one = _scatter_offsets(ldsr_plot(chi2, ld, fit, seed=1))
        two = _scatter_offsets(ldsr_plot(chi2, ld, fit, seed=2))

        assert not np.array_equal(one, two)

    def test_global_random_state_is_left_untouched(self) -> None:
        chi2, ld, fit = _ldsr_inputs()
        random.seed(2026)
        before = random.getstate()

        plt.close(ldsr_plot(chi2, ld, fit))

        assert random.getstate() == before


# ---------------------------------------------------------------------------
# power.subsample_convergence
# ---------------------------------------------------------------------------


def _convergence_inputs() -> tuple[list[list[int]], list[float]]:
    rng = np.random.default_rng(11)
    n_variants, n_samples = 40, 24
    genotypes = [[int(v) for v in rng.integers(0, 3, n_samples)] for _ in range(n_variants)]
    traits = [float(v) for v in rng.normal(0.0, 1.0, n_samples)]
    return genotypes, traits


class TestSubsampleConvergence:
    KWARGS = {"fractions": [0.25, 0.5], "n_replicates": 4, "alpha": 0.05}

    def test_same_seed_gives_the_same_replicates_whatever_the_global_state(self) -> None:
        genotypes, traits = _convergence_inputs()

        first = _under_global_seed(1, subsample_convergence, genotypes, traits, seed=3, **self.KWARGS)
        second = _under_global_seed(2, subsample_convergence, genotypes, traits, seed=3, **self.KWARGS)

        assert first["status"] == "success"
        assert first == second

    def test_default_seed_is_reproducible(self) -> None:
        genotypes, traits = _convergence_inputs()

        first = _under_global_seed(21, subsample_convergence, genotypes, traits, **self.KWARGS)
        second = _under_global_seed(22, subsample_convergence, genotypes, traits, **self.KWARGS)

        assert first == second

    def test_different_seeds_draw_different_subsamples(self) -> None:
        genotypes, traits = _convergence_inputs()

        one = subsample_convergence(genotypes, traits, seed=1, **self.KWARGS)
        two = subsample_convergence(genotypes, traits, seed=2, **self.KWARGS)

        assert one["metrics"]["mean_neg_log10_p"] != two["metrics"]["mean_neg_log10_p"]

    def test_global_random_state_is_left_untouched(self) -> None:
        genotypes, traits = _convergence_inputs()
        random.seed(2027)
        before = random.getstate()

        subsample_convergence(genotypes, traits, **self.KWARGS)

        assert random.getstate() == before

    def test_the_full_data_fraction_matches_a_direct_recomputation(self) -> None:
        """Fraction 1.0 involves no sampling: every replicate is the full scan, which can be redone directly."""
        genotypes, traits = _convergence_inputs()
        alpha = 0.2

        result = subsample_convergence(genotypes, traits, fractions=[1.0], n_replicates=3, alpha=alpha, seed=8)

        p_values = [association_test_linear(row, traits)["p_value"] for row in genotypes]
        expected_hits = sum(1 for p in p_values if p < alpha)
        expected_mean_nlp = sum(-np.log10(max(p, 1e-300)) for p in p_values) / len(p_values)
        assert 0 < expected_hits < len(p_values)  # a non-trivial count, so the comparison means something
        hits = result["metrics"]["n_significant"][0]
        assert hits["values"] == [expected_hits] * 3
        assert hits["std"] == 0.0
        assert result["metrics"]["mean_neg_log10_p"][0]["mean"] == pytest.approx(expected_mean_nlp, rel=1e-12)

    def test_empty_matrix_is_still_an_error_result(self) -> None:
        assert subsample_convergence([], [])["status"] == "error"
