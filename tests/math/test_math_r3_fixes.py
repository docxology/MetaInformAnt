"""Regression tests for math defects fixed in the round-3 audit.

REAL IMPLEMENTATION POLICY: every test runs the real functions and checks them against independent
oracles (closed forms from the literature, scipy, and direct simulation); nothing is mocked or patched.

Covered here:
- ``effective_size_from_family_size_variance`` died on an integer census size with a bare AssertionError.
- ``effective_size_sex_ratio`` jumped by a factor of two at a sex ratio of exactly 0.5 and exceeded the
  census size away from it.
- ``fisher_exact_test`` held a made-up p-value (``chi_square / 3.84``) in a fallback branch.
- ``ld_decay_r2`` estimated the decay rate from the first and last bin of a coarse binning.
- ``ddm_mean_decision_time`` / ``ddm_analytic_accuracy`` used rough stand-ins for the closed forms.
- ``metropolis_hastings``, ``simulate_coalescent`` and ``calculate_confidence_intervals`` drew from, or
  reseeded, a global random generator with no way to seed them locally.
- ``compare_statistics`` ran Welch's t-test whatever ``test_type`` said (``"kruskal"`` included), echoed
  the requested name, and turned any SciPy error into a p-value of 1.0.
- ``calculate_confidence_intervals`` answered a dictionary of point values with an invented +-10% band.
"""

from __future__ import annotations

import itertools
import math
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats
from scipy.optimize import curve_fit

import metainformant
from metainformant.math.bayesian.inference import metropolis_hastings
from metainformant.math.core.utilities import fisher_exact_test
from metainformant.math.decision_theory.ddm import ddm_analytic_accuracy, ddm_mean_decision_time
from metainformant.math.population_genetics.coalescent import simulate_coalescent
from metainformant.math.population_genetics.effective_size import effective_size_sex_ratio
from metainformant.math.population_genetics.ld import ld_decay_r2
from metainformant.math.population_genetics.statistics import (
    calculate_confidence_intervals,
    compare_statistics,
    effective_size_from_family_size_variance,
)

# ---------------------------------------------------------------------------
# effective_size_from_family_size_variance
# ---------------------------------------------------------------------------


class TestFamilySizeVariance:
    def test_integer_census_size_is_accepted(self) -> None:
        """``effective_size_from_family_size_variance(1000, Vk=5)`` used to raise a bare AssertionError."""
        assert effective_size_from_family_size_variance(1000, Vk=5) == pytest.approx((4 * 1000 - 2) / (5 + 2))

    @pytest.mark.parametrize("census", [1000, 1000.0, np.int64(1000), np.float32(1000)])
    def test_int_float_and_numpy_sizes_agree(self, census: float) -> None:
        assert effective_size_from_family_size_variance(census, Vk=2.0) == pytest.approx((4 * 1000 - 2) / (2.0 + 2))

    @pytest.mark.parametrize("census,vk", [(50, 0.0), (1000, 2.0), (1000, 5.0), (12, 7.5)])
    def test_matches_crow_kimura_general_formula(self, census: int, vk: float) -> None:
        """For a stable population (mean offspring number 2) Crow and Kimura's general
        (N k - 1) / (k - 1 + Vk / k) reduces to (4N - 2) / (Vk + 2)."""
        k = 2.0
        general = (census * k - 1.0) / (k - 1.0 + vk / k)
        assert effective_size_from_family_size_variance(census, Vk=vk) == pytest.approx(general)

    def test_list_form_is_unchanged(self) -> None:
        assert effective_size_from_family_size_variance([1, 5]) == pytest.approx(1.0)
        assert effective_size_from_family_size_variance([2, 2, 2, 2]) == float("inf")
        assert effective_size_from_family_size_variance([]) == 0.0

    def test_tuple_and_array_offspring_counts_are_accepted(self) -> None:
        assert effective_size_from_family_size_variance((1, 5)) == pytest.approx(1.0)
        assert effective_size_from_family_size_variance(np.array([1, 5])) == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "args",
        [([2, 3, 2], 5.0), ("1000", 5.0), (True, 5.0), (100.0, None), (None, None)],
    )
    def test_wrong_types_raise_type_error_not_assertion_error(self, args: tuple) -> None:
        with pytest.raises(TypeError):
            effective_size_from_family_size_variance(*args)

    @pytest.mark.parametrize("census,vk", [(0, 1.0), (-5, 1.0), (float("nan"), 1.0), (100, -1.0), (100, float("nan"))])
    def test_impossible_values_raise_value_error(self, census: float, vk: float) -> None:
        with pytest.raises(ValueError):
            effective_size_from_family_size_variance(census, Vk=vk)


# ---------------------------------------------------------------------------
# effective_size_sex_ratio
# ---------------------------------------------------------------------------


class TestSexRatioEffectiveSize:
    @pytest.mark.parametrize("males,females", [(100, 100), (40, 60), (10, 90), (1, 500), (300.5, 7.25)])
    def test_wright_formula(self, males: float, females: float) -> None:
        assert effective_size_sex_ratio(males, females) == pytest.approx(4 * males * females / (males + females))

    def test_continuous_through_a_sex_ratio_of_one_half(self) -> None:
        """Ne used to jump from 400.00001 at a sex ratio of 0.4999 to 200 at exactly 0.5."""
        at_half = effective_size_sex_ratio(100, 100, sex_ratio=0.5)
        assert at_half == pytest.approx(200.0)
        for ratio in (0.4999, 0.49999999, 0.5001, 0.45, 0.3, 0.9):
            assert effective_size_sex_ratio(100, 100, sex_ratio=ratio) == pytest.approx(at_half, rel=1e-9)

    def test_never_exceeds_the_census_size(self) -> None:
        """Unequal sexes can only lower Ne; the old adjustment gave 800 for 100 males and 300 females."""
        for males in (1, 5, 40, 100, 1000):
            for females in (1, 7, 60, 300, 5000):
                for ratio in (0.01, 0.1, 0.25, 0.5, 0.75, 0.99):
                    ne = effective_size_sex_ratio(males, females, sex_ratio=ratio)
                    assert ne <= males + females + 1e-9, (males, females, ratio, ne)

    def test_symmetric_in_the_two_sexes(self) -> None:
        assert effective_size_sex_ratio(10, 90) == pytest.approx(effective_size_sex_ratio(90, 10))

    @pytest.mark.parametrize("ratio", [None, 0.1])
    def test_matches_allele_frequency_drift_simulation(self, ratio: float | None) -> None:
        """One generation of drift: p' = (p_m + p_f) / 2 with binomial samples from 2 Nm male and 2 Nf female genes.

        Var(p') = p q / (2 Ne) defines the variance effective size. The sex ratio of the breeders is 0.1
        here, so passing it (the consistent value) must not change the answer.
        """
        males, females, p, replicates = 10, 90, 0.5, 400_000
        rng = np.random.default_rng(3)
        male_freq = rng.binomial(2 * males, p, replicates) / (2 * males)
        female_freq = rng.binomial(2 * females, p, replicates) / (2 * females)
        offspring = (male_freq + female_freq) / 2
        ne_from_drift = p * (1 - p) / (2 * offspring.var(ddof=1))
        kwargs = {} if ratio is None else {"sex_ratio": ratio}
        assert effective_size_sex_ratio(males, females, **kwargs) == pytest.approx(ne_from_drift, rel=0.02)

    def test_validation_is_unchanged(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            effective_size_sex_ratio(0.0, 100.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            effective_size_sex_ratio(100.0, 100.0, sex_ratio=1.0)


# ---------------------------------------------------------------------------
# fisher_exact_test
# ---------------------------------------------------------------------------


class TestFisherExactTest:
    @pytest.mark.parametrize("table", [(8, 2, 1, 9), (10, 2, 3, 15), (5, 0, 0, 5), (1, 1, 1, 1), (0, 1, 1, 1)])
    def test_matches_scipy_exact_test(self, table: tuple[int, int, int, int]) -> None:
        a, b, c, d = table
        expected = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        odds_ratio, p_value = fisher_exact_test(a, b, c, d)
        assert p_value == pytest.approx(float(expected.pvalue), rel=1e-12)
        assert odds_ratio == pytest.approx(float(expected.statistic), rel=1e-12)

    def test_empty_margin_gives_p_value_one(self) -> None:
        odds_ratio, p_value = fisher_exact_test(0, 0, 3, 4)
        assert math.isnan(odds_ratio)
        assert p_value == 1.0

    def test_without_scipy_the_function_raises_instead_of_inventing_a_p_value(self, tmp_path: Path) -> None:
        """A fallback branch used to return ``min(1.0, chi_square / 3.84)``, which grows as the evidence strengthens.

        scipy is blocked after the package is imported (as an uninstalled scipy would be at call time) in a
        fresh interpreter; the function must fail loudly, never hand back a number.
        """
        code = (
            "import sys\n"
            "from metainformant.math.core.utilities import fisher_exact_test\n"
            "sys.modules['scipy'] = None\n"
            "sys.modules['scipy.stats'] = None\n"
            "try:\n"
            "    result = fisher_exact_test(20, 0, 0, 20)\n"
            "except ImportError as exc:\n"
            "    print('IMPORT_ERROR', exc)\n"
            "else:\n"
            "    print('RETURNED', result)\n"
        )
        # The test session prepends the repository's src to PYTHONPATH for children, so pin the child to the
        # very package this process imported.
        package_root = str(Path(metainformant.__file__).resolve().parent.parent)
        env = {
            **os.environ,
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONPATH": package_root + os.pathsep + os.environ.get("PYTHONPATH", ""),
        }
        proc = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, cwd=tmp_path, env=env, timeout=120
        )
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.startswith("IMPORT_ERROR"), proc.stdout
        assert "scipy is required" in proc.stdout


# ---------------------------------------------------------------------------
# ld_decay_r2: exponential fit
# ---------------------------------------------------------------------------


class TestLdDecayFit:
    def test_recovers_exponential_parameters_exactly(self) -> None:
        distances = np.linspace(1, 5000, 400)
        fit = ld_decay_r2(list(distances), list(0.8 * np.exp(-0.002 * distances)))
        assert fit["decay_rate"] == pytest.approx(0.002, rel=1e-6)
        assert fit["initial_r2"] == pytest.approx(0.8, rel=1e-6)
        assert fit["half_decay_distance"] == pytest.approx(math.log(2) / 0.002, rel=1e-6)

    def test_agrees_with_scipy_curve_fit_on_noisy_pairs(self) -> None:
        rng = np.random.default_rng(7)
        distances = rng.uniform(1, 3000, 800)
        r2 = np.clip(0.7 * np.exp(-0.0015 * distances) + rng.normal(0, 0.03, distances.size), 0, 1)
        fit = ld_decay_r2(list(distances), list(r2))
        (a, k), _ = curve_fit(
            lambda x, a, k: a * np.exp(-k * x), distances, r2, p0=[0.5, 1e-3], bounds=([0, 0], [np.inf, np.inf])
        )
        assert fit["decay_rate"] == pytest.approx(k, rel=1e-4)
        assert fit["initial_r2"] == pytest.approx(a, rel=1e-4)
        # ... and lands near the generating values
        assert fit["decay_rate"] == pytest.approx(0.0015, rel=0.1)

    def test_pairs_closer_than_one_unit_are_used(self) -> None:
        """Distances in Mb: the old binning started at 1 and discarded most of these pairs."""
        distances = np.linspace(0.01, 3.0, 200)
        fit = ld_decay_r2(list(distances), list(0.6 * np.exp(-3.0 * distances)))
        assert fit["decay_rate"] == pytest.approx(3.0, rel=1e-6)
        assert fit["initial_r2"] == pytest.approx(0.6, rel=1e-6)

    def test_rate_scales_inversely_with_the_distance_unit(self) -> None:
        distances = np.linspace(1, 100, 80)
        r2 = 0.5 * np.exp(-0.05 * distances)
        per_unit = ld_decay_r2(list(distances), list(r2))["decay_rate"]
        per_thousand = ld_decay_r2(list(distances * 1000.0), list(r2))["decay_rate"]
        assert per_thousand == pytest.approx(per_unit / 1000.0, rel=1e-6)

    def test_geometric_decay_has_rate_minus_log_ratio(self) -> None:
        """r-squared = 0.9 * 0.8 ** d is exponential with k = -ln(0.8)."""
        distances = np.arange(1, 40, dtype=float)
        fit = ld_decay_r2(list(distances), list(0.9 * 0.8**distances))
        assert fit["decay_rate"] == pytest.approx(-math.log(0.8), rel=1e-6)

    def test_max_distance_excludes_distant_pairs(self) -> None:
        distances = np.linspace(1, 6000, 600)
        r2 = 0.8 * np.exp(-0.002 * distances)
        r2[distances > 3000] = 0.5  # contamination beyond the cut-off
        fit = ld_decay_r2(list(distances), list(r2), max_distance=3000)
        assert fit["decay_rate"] == pytest.approx(0.002, rel=1e-6)

    @pytest.mark.parametrize("r2", [[0.2, 0.2, 0.2, 0.2], [0.1, 0.2, 0.3, 0.4]])
    def test_flat_or_rising_data_report_no_decay(self, r2: list[float]) -> None:
        fit = ld_decay_r2([1.0, 10.0, 100.0, 1000.0], r2)
        assert fit["decay_rate"] == 0.0
        assert fit["half_decay_distance"] == float("inf")
        assert fit["initial_r2"] == pytest.approx(float(np.mean(r2)))

    def test_all_zero_r_squared_report_no_decay(self) -> None:
        fit = ld_decay_r2([1.0, 10.0, 100.0], [0.0, 0.0, 0.0])
        assert fit == {"decay_rate": 0.0, "half_decay_distance": float("inf"), "initial_r2": 0.0}

    def test_too_few_distinct_distances_give_the_empty_result(self) -> None:
        assert ld_decay_r2([1.0, 2.0], [0.5, 0.2]) == {"decay_rate": 0.0, "half_decay_distance": float("inf")}
        assert ld_decay_r2([5.0, 5.0, 5.0, 5.0], [0.5, 0.4, 0.3, 0.2]) == {
            "decay_rate": 0.0,
            "half_decay_distance": float("inf"),
        }
        assert ld_decay_r2([1.0, 2.0, 3.0], [0.5, 0.4, 0.3], max_distance=0.5) == {
            "decay_rate": 0.0,
            "half_decay_distance": float("inf"),
        }

    def test_invalid_inputs_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            ld_decay_r2([1.0, 2.0, 3.0], [0.5, 0.4])
        with pytest.raises(ValueError, match="non-negative"):
            ld_decay_r2([-1.0, 2.0, 3.0], [0.5, 0.4, 0.3])
        with pytest.raises(ValueError, match="finite"):
            ld_decay_r2([1.0, 2.0, 3.0], [0.5, float("nan"), 0.3])

    def test_scalar_form_is_unchanged(self) -> None:
        assert ld_decay_r2(0.8, recombination_rate=0.1, generations=5) == pytest.approx(0.8 * 0.9**10)


# ---------------------------------------------------------------------------
# DDM closed forms against simulation
# ---------------------------------------------------------------------------


def _simulate_ddm(drift: float, boundary: float, noise: float, n_trials: int, dt: float, seed: int):
    """Euler-Maruyama paths from the midpoint until they leave (-a/2, a/2); returns (times, hit_upper).

    Stepping misses excursions between steps, so each barrier is pulled in by 0.5826 s sqrt(dt)
    (Broadie, Glasserman and Kou 1997), which makes the stepped walk match the continuous process to O(dt).
    """
    rng = np.random.default_rng(seed)
    inner = boundary / 2 - 0.5826 * noise * math.sqrt(dt)
    step_sd = noise * math.sqrt(dt)
    position = np.zeros(n_trials)
    elapsed = np.zeros(n_trials)
    hit_upper = np.zeros(n_trials, dtype=bool)
    active = np.arange(n_trials)
    while active.size:
        position[active] += drift * dt + step_sd * rng.standard_normal(active.size)
        elapsed[active] += dt
        upper = position[active] >= inner
        lower = position[active] <= -inner
        hit_upper[active[upper]] = True
        active = active[~(upper | lower)]
    return elapsed, hit_upper


class TestDdmClosedForms:
    CASES = [(1.0, 2.0, 1.0), (0.5, 1.5, 2.0), (0.0, 1.0, 1.0), (-0.8, 1.2, 0.7)]

    @pytest.mark.parametrize("drift,boundary,noise", CASES)
    def test_closed_forms_match_simulation(self, drift: float, boundary: float, noise: float) -> None:
        n_trials = 40000
        times, upper = _simulate_ddm(drift, boundary, noise, n_trials, dt=1e-3, seed=11)

        accuracy = ddm_analytic_accuracy(drift, boundary, noise)
        assert stats.binomtest(int(upper.sum()), n_trials, accuracy).pvalue > 1e-3

        expected_time = ddm_mean_decision_time(drift, boundary, noise)
        standard_error = times.std(ddof=1) / math.sqrt(n_trials)
        assert abs(times.mean() - expected_time) < 4.5 * standard_error

    def test_mean_decision_time_known_values(self) -> None:
        # (a / 2v) tanh(v a / 2 s^2) at v = 1, a = 2, s = 1 is tanh(1)
        assert ddm_mean_decision_time(1.0, 2.0) == pytest.approx(math.tanh(1.0))
        # Bogacz et al. (2006) with thresholds +-z: (z / A) tanh(A z / c^2), z = a / 2
        z, a_drift, c = 0.9, 0.4, 1.3
        assert ddm_mean_decision_time(a_drift, 2 * z, c) == pytest.approx((z / a_drift) * math.tanh(a_drift * z / c**2))

    def test_accuracy_known_values(self) -> None:
        # Bogacz et al. (2006): error rate 1 / (1 + exp(2 A z / c^2)) with thresholds +-z
        z, a_drift, c = 0.9, 0.4, 1.3
        error_rate = 1.0 / (1.0 + math.exp(2 * a_drift * z / c**2))
        assert ddm_analytic_accuracy(a_drift, 2 * z, c) == pytest.approx(1.0 - error_rate)

    def test_no_drift_limit_is_the_brownian_exit_time(self) -> None:
        assert ddm_mean_decision_time(0.0, 3.0, 1.5) == pytest.approx(3.0**2 / (4 * 1.5**2))
        assert ddm_mean_decision_time(1e-12, 3.0, 1.5) == pytest.approx(ddm_mean_decision_time(0.0, 3.0, 1.5))

    def test_symmetry_in_the_drift_sign(self) -> None:
        for drift in (0.3, 1.0, 2.5):
            assert ddm_mean_decision_time(-drift, 1.7, 0.8) == pytest.approx(ddm_mean_decision_time(drift, 1.7, 0.8))
            assert ddm_analytic_accuracy(-drift, 1.7, 0.8) == pytest.approx(
                1.0 - ddm_analytic_accuracy(drift, 1.7, 0.8)
            )
        assert ddm_analytic_accuracy(0.0, 1.7, 0.8) == 0.5

    def test_noise_enters_through_the_scale_invariance_of_the_process(self) -> None:
        """Scaling v, a and s together rescales the state space only, so time and accuracy are unchanged."""
        for factor in (0.1, 3.0, 40.0):
            assert ddm_mean_decision_time(0.7 * factor, 1.4 * factor, 0.9 * factor) == pytest.approx(
                ddm_mean_decision_time(0.7, 1.4, 0.9)
            )
            assert ddm_analytic_accuracy(0.7 * factor, 1.4 * factor, 0.9 * factor) == pytest.approx(
                ddm_analytic_accuracy(0.7, 1.4, 0.9)
            )

    def test_extreme_drift_saturates_without_overflow(self) -> None:
        assert ddm_analytic_accuracy(-1000.0, 2.0) == 0.0
        assert ddm_analytic_accuracy(1000.0, 2.0) == 1.0
        assert ddm_mean_decision_time(1000.0, 2.0) == pytest.approx(2.0 / (2 * 1000.0), rel=1e-9)

    @pytest.mark.parametrize("kwargs", [{"noise_sd": 0.0}, {"noise_sd": -1.0}, {"boundary": 0.0}])
    def test_invalid_parameters_raise(self, kwargs: dict) -> None:
        params = {"drift_rate": 0.5, "boundary": 1.0, "noise_sd": 1.0} | kwargs
        with pytest.raises(ValueError):
            ddm_mean_decision_time(**params)
        with pytest.raises(ValueError):
            ddm_analytic_accuracy(**params)

    @pytest.mark.parametrize("drift", [float("nan"), float("inf")])
    def test_non_finite_drift_raises(self, drift: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            ddm_mean_decision_time(drift, 1.0)
        with pytest.raises(ValueError, match="finite"):
            ddm_analytic_accuracy(drift, 1.0)


# ---------------------------------------------------------------------------
# Seeding: no unseedable global randomness
# ---------------------------------------------------------------------------


def _standard_normal_log_posterior(params: list[float]) -> float:
    return -0.5 * params[0] ** 2


def _same_numpy_state(first: tuple, second: tuple) -> bool:
    return first[0] == second[0] and np.array_equal(first[1], second[1]) and first[2:] == second[2:]


class TestMetropolisHastingsSeeding:
    def run(self, **kwargs):
        return metropolis_hastings(_standard_normal_log_posterior, [0.0], n_samples=400, burn_in=100, **kwargs)

    def test_same_seed_gives_the_same_chain_whatever_the_global_state(self) -> None:
        first = self.run(random_state=11)
        random.seed(98765)  # a different global state must not matter
        second = self.run(random_state=11)
        assert first["samples"] == second["samples"]
        assert first["acceptance_rate"] == second["acceptance_rate"]

    def test_different_seeds_give_different_chains(self) -> None:
        assert self.run(random_state=1)["samples"] != self.run(random_state=2)["samples"]

    def test_seeded_run_leaves_the_global_generator_untouched(self) -> None:
        random.seed(5)
        before = random.getstate()
        self.run(random_state=3)
        assert random.getstate() == before

    def test_random_instance_is_used_as_given(self) -> None:
        by_seed = self.run(random_state=21)
        by_instance = self.run(random_state=random.Random(21))
        assert by_instance["samples"] == by_seed["samples"]
        rng = random.Random(21)
        self.run(random_state=rng)
        assert rng.getstate() != random.Random(21).getstate()  # the caller's generator advanced

    def test_chain_targets_the_posterior(self) -> None:
        """A seeded chain is a valid sample: it passes a Kolmogorov-Smirnov test against N(0, 1) after thinning."""
        chain = metropolis_hastings(
            _standard_normal_log_posterior, [0.0], n_samples=60000, burn_in=2000, proposal_scale=2.0, random_state=17
        )
        draws = np.array([s[0] for s in chain["samples"]])[::40]  # thin to near-independent draws
        assert stats.kstest(draws, "norm").pvalue > 1e-3

    def test_default_still_follows_the_global_generator(self) -> None:
        random.seed(20260905)
        first = self.run()
        random.seed(20260905)
        second = self.run()
        assert first["samples"] == second["samples"]


class TestCoalescentSeeding:
    def test_same_seed_gives_the_same_genealogy_whatever_the_global_state(self) -> None:
        np.random.seed(1)
        first = simulate_coalescent(8, effective_size=1000.0, mutation_rate=1e-3, random_state=5)
        np.random.seed(999)
        second = simulate_coalescent(8, effective_size=1000.0, mutation_rate=1e-3, random_state=5)
        assert first == second

    def test_generator_and_random_state_instances_are_accepted(self) -> None:
        by_seed = simulate_coalescent(8, effective_size=1000.0, mutation_rate=1e-3, random_state=5)
        by_generator = simulate_coalescent(
            8, effective_size=1000.0, mutation_rate=1e-3, random_state=np.random.default_rng(5)
        )
        assert by_generator == by_seed
        legacy = simulate_coalescent(
            8, effective_size=1000.0, mutation_rate=1e-3, random_state=np.random.RandomState(5)
        )
        assert legacy["n_samples"] == 8 and len(legacy["coalescence_times"]) == 7

    def test_seeded_run_leaves_the_global_generator_untouched(self) -> None:
        np.random.seed(7)
        before = np.random.get_state()
        simulate_coalescent(10, random_state=3)
        assert _same_numpy_state(before, np.random.get_state())

    def test_default_still_follows_the_global_generator(self) -> None:
        np.random.seed(7)
        first = simulate_coalescent(6)
        np.random.seed(7)
        second = simulate_coalescent(6)
        assert first == second

    def test_expected_time_to_mrca_matches_theory(self) -> None:
        """E[T_MRCA] = 4 Ne (1 - 1/n) generations for the diploid coalescent; the seeded mean must agree."""
        n, ne, replicates = 5, 100.0, 4000
        times = np.array(
            [simulate_coalescent(n, effective_size=ne, random_state=seed)["total_time"] for seed in range(replicates)]
        )
        expected = 4 * ne * (1 - 1 / n)
        assert abs(times.mean() - expected) < 4.5 * times.std(ddof=1) / math.sqrt(replicates)


class TestBootstrapConfidenceIntervalSeeding:
    VALUES = [1.0, 2.0, 3.5, 4.0, 5.5, 7.0]

    def test_does_not_reseed_the_global_numpy_generator(self) -> None:
        """The bootstrap used ``np.random.seed(42)``, silently resetting every caller's random stream."""
        np.random.seed(2026)
        before = np.random.get_state()
        calculate_confidence_intervals(self.VALUES, method="bootstrap")
        assert _same_numpy_state(before, np.random.get_state())

    def test_interval_is_unchanged_and_reproducible(self) -> None:
        # Values recorded from the implementation that seeded the global generator with 42
        expected = (2.3333333333333335, 5.504166666666666)
        np.random.seed(1)
        assert calculate_confidence_intervals(self.VALUES, method="bootstrap") == pytest.approx(expected, rel=1e-12)
        np.random.seed(2)
        assert calculate_confidence_intervals(self.VALUES, method="bootstrap") == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# compare_statistics: the requested test is the test that runs
# ---------------------------------------------------------------------------

GROUP_A = [4.1, 5.3, 6.2, 5.8, 4.9]
GROUP_B = [6.8, 7.4, 5.95, 8.1, 7.7, 6.5]
COMPARISON_TYPES = ("ttest", "mannwhitney", "kruskal")


def _welch_by_hand(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Welch's t statistic and two-sided p-value from the closed form (Welch-Satterthwaite degrees of freedom)."""
    n1, n2 = a.size, b.size
    v1, v2 = a.var(ddof=1), b.var(ddof=1)
    se2 = v1 / n1 + v2 / n2
    t = (a.mean() - b.mean()) / math.sqrt(se2)
    df = se2**2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    return float(t), float(2 * stats.t.sf(abs(t), df))


def _kruskal_by_hand(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Kruskal-Wallis H from pooled mid-ranks with the tie correction, and its chi-square(1) p-value."""
    pooled = np.concatenate([a, b])
    n = pooled.size
    ranks = stats.rankdata(pooled)
    rank_sum_a, rank_sum_b = ranks[: a.size].sum(), ranks[a.size :].sum()
    h = 12 / (n * (n + 1)) * (rank_sum_a**2 / a.size + rank_sum_b**2 / b.size) - 3 * (n + 1)
    ties = np.unique(pooled, return_counts=True)[1]
    h /= 1 - np.sum(ties**3 - ties) / (n**3 - n)
    return float(h), float(stats.chi2.sf(h, 1))


def _mann_whitney_by_enumeration(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """U of the first sample and the exact two-sided p-value from all splits of the pooled ranks (no ties)."""
    n1, n2 = a.size, b.size
    ranks = stats.rankdata(np.concatenate([a, b]))
    u_observed = ranks[:n1].sum() - n1 * (n1 + 1) / 2
    every_u = np.array([sum(chosen) - n1 * (n1 + 1) / 2 for chosen in itertools.combinations(ranks, n1)])
    tail = np.mean(every_u <= min(u_observed, n1 * n2 - u_observed))
    return float(u_observed), float(min(1.0, 2 * tail))


class TestCompareStatisticsRunsTheRequestedTest:
    a = np.array(GROUP_A)
    b = np.array(GROUP_B)

    @pytest.mark.parametrize(
        "test_type,reference",
        [
            ("ttest", lambda a, b: stats.ttest_ind(a, b, equal_var=False)),
            ("mannwhitney", lambda a, b: stats.mannwhitneyu(a, b, alternative="two-sided")),
            ("kruskal", lambda a, b: stats.kruskal(a, b)),
        ],
    )
    def test_each_advertised_test_matches_scipy(self, test_type: str, reference) -> None:
        expected = reference(self.a, self.b)
        result = compare_statistics(GROUP_A, GROUP_B, test_type=test_type)
        assert result["test_type"] == test_type
        assert result["test_statistic"] == pytest.approx(float(expected.statistic), rel=1e-12)
        assert result["p_value"] == pytest.approx(float(expected.pvalue), rel=1e-12)

    def test_kruskal_is_the_kruskal_wallis_statistic_not_the_t_statistic(self) -> None:
        """``test_type="kruskal"`` used to return Welch's t (about -3.71 here) under the name 'kruskal'."""
        result = compare_statistics(GROUP_A, GROUP_B, test_type="kruskal")
        h, p = _kruskal_by_hand(self.a, self.b)
        assert result["test_statistic"] == pytest.approx(h, rel=1e-12)
        assert result["p_value"] == pytest.approx(p, rel=1e-12)
        assert result["test_statistic"] > 0  # an H statistic is never negative, a t statistic can be
        assert result["test_statistic"] != pytest.approx(compare_statistics(GROUP_A, GROUP_B)["test_statistic"])

    def test_welch_ttest_matches_the_closed_form(self) -> None:
        t, p = _welch_by_hand(self.a, self.b)
        result = compare_statistics(GROUP_A, GROUP_B, test_type="ttest")
        assert result["test_statistic"] == pytest.approx(t, rel=1e-12)
        assert result["p_value"] == pytest.approx(p, rel=1e-12)

    def test_ttest_is_the_default(self) -> None:
        assert compare_statistics(GROUP_A, GROUP_B) == compare_statistics(GROUP_A, GROUP_B, test_type="ttest")

    def test_mann_whitney_matches_the_exact_enumeration(self) -> None:
        u, p = _mann_whitney_by_enumeration(self.a, self.b)
        result = compare_statistics(GROUP_A, GROUP_B, test_type="mannwhitney")
        assert result["test_statistic"] == pytest.approx(u)
        assert result["p_value"] == pytest.approx(p, rel=1e-12)

    def test_the_three_tests_give_three_different_answers(self) -> None:
        statistics = {t: compare_statistics(GROUP_A, GROUP_B, test_type=t)["test_statistic"] for t in COMPARISON_TYPES}
        assert len({round(value, 6) for value in statistics.values()}) == 3

    def test_inputs_may_be_arrays_or_lists(self) -> None:
        for test_type in COMPARISON_TYPES:
            assert compare_statistics(self.a, self.b, test_type=test_type) == compare_statistics(
                GROUP_A, GROUP_B, test_type=test_type
            )

    @pytest.mark.parametrize("bad", ["anova", "Kruskal", "wilcoxon", "welch", ""])
    def test_unknown_test_type_raises(self, bad: str) -> None:
        """An unknown name used to run Welch's t-test and be echoed back as if it had been run."""
        with pytest.raises(ValueError, match="Unknown test_type"):
            compare_statistics(GROUP_A, GROUP_B, test_type=bad)
        with pytest.raises(ValueError, match="Unknown test_type"):
            compare_statistics({"pi": 0.5}, {"pi": 0.1}, test_type=bad)

    @pytest.mark.parametrize("test_type", COMPARISON_TYPES)
    @pytest.mark.parametrize(
        "first,second",
        [
            ([], [1.0, 2.0, 3.0]),
            ([1.0, 2.0, 3.0], []),
            ([1.0, float("nan"), 3.0], [2.0, 3.0, 4.0]),
            ([1.0, float("inf"), 3.0], [2.0, 3.0, 4.0]),
            ([[1.0, 2.0], [3.0, 4.0]], [1.0, 2.0, 3.0]),
        ],
    )
    def test_untestable_input_raises_instead_of_returning_a_p_value(
        self, test_type: str, first: list, second: list
    ) -> None:
        """SciPy errors used to be swallowed into ``p_value = 1.0`` (and NaN data into a NaN p-value)."""
        with pytest.raises(ValueError):
            compare_statistics(first, second, test_type=test_type)

    def test_single_observation_cannot_be_compared_with_welch(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            compare_statistics([1.0], [2.0, 3.0, 4.0], test_type="ttest")

    def test_constant_groups_make_welch_undefined(self) -> None:
        with pytest.raises(ValueError, match="constant"):
            compare_statistics([2.0, 2.0, 2.0], [2.0, 2.0, 2.0], test_type="ttest")

    @pytest.mark.parametrize("test_type", ["mannwhitney", "kruskal"])
    def test_identical_observations_make_a_rank_test_undefined(self, test_type: str) -> None:
        with pytest.raises(ValueError, match="identical"):
            compare_statistics([2.0, 2.0, 2.0], [2.0, 2.0], test_type=test_type)

    def test_one_observation_per_group_is_enough_for_the_rank_tests(self) -> None:
        for test_type in ("mannwhitney", "kruskal"):
            result = compare_statistics([1.0], [2.0], test_type=test_type)
            assert 0.0 <= result["p_value"] <= 1.0

    @pytest.mark.parametrize("test_type", ["mannwhitney", "kruskal"])
    def test_dictionaries_hold_no_observations_to_rank(self, test_type: str) -> None:
        """A dictionary comparison is the approximate normal one; it must not be passed off as a rank test."""
        with pytest.raises(ValueError, match="ranks observations"):
            compare_statistics({"pi": 0.5}, {"pi": 0.1}, test_type=test_type)

    def test_dictionary_comparison_still_works_for_the_default_test(self) -> None:
        result = compare_statistics({"pi": 0.5, "variance": 0.01}, {"pi": 0.1, "variance": 0.01})
        assert result["pi"]["significant"] is True
        assert result["pi"]["difference"] == pytest.approx(0.4)


# ---------------------------------------------------------------------------
# calculate_confidence_intervals: intervals come from data, never from a made-up band
# ---------------------------------------------------------------------------


class TestConfidenceIntervalsComeFromData:
    OBSERVATIONS = {
        "pi": [0.0090, 0.0112, 0.0101, 0.0124, 0.0083, 0.0107, 0.0096, 0.0118],
        "theta": [0.0201, 0.0243, 0.0178, 0.0225, 0.0212, 0.0190, 0.0236, 0.0207],
    }

    @pytest.mark.parametrize("points", [{"pi": 0.01, "theta": 0.01}, {"d": -0.5}, {"x": 0.0}, {"y": np.float64(0.25)}])
    def test_dictionary_of_point_values_is_rejected(self, points: dict) -> None:
        """A lone estimate holds no sampling information; the +-10% band it used to get was invented."""
        with pytest.raises(ValueError, match="no sampling information"):
            calculate_confidence_intervals(points)

    def test_the_error_names_the_offending_statistic(self) -> None:
        with pytest.raises(ValueError, match="'theta'"):
            calculate_confidence_intervals({"pi": [0.1, 0.2, 0.3], "theta": 0.01})

    @pytest.mark.parametrize("method", ["bootstrap", "normal", "percentile"])
    def test_dictionary_of_observations_gives_the_interval_of_each_list(self, method: str) -> None:
        intervals = calculate_confidence_intervals(self.OBSERVATIONS, method=method)
        assert set(intervals) == {"pi", "theta"}
        for name, observations in self.OBSERVATIONS.items():
            lower, upper = calculate_confidence_intervals(observations, method=method)
            assert intervals[name] == {"ci_lower": lower, "ci_upper": upper}
            assert lower < upper

    def test_normal_interval_matches_mean_plus_minus_z_standard_errors(self) -> None:
        for name, observations in self.OBSERVATIONS.items():
            interval = calculate_confidence_intervals({name: observations}, confidence_level=0.9, method="normal")[name]
            centre, half_width = np.mean(observations), stats.norm.ppf(0.95) * stats.sem(observations)
            assert interval["ci_lower"] == pytest.approx(centre - half_width, rel=1e-12)
            assert interval["ci_upper"] == pytest.approx(centre + half_width, rel=1e-12)

    def test_interval_is_estimated_not_proportional_to_the_value(self) -> None:
        """Same mean, different spread: the width follows the data. A fixed +-10% band could not."""
        tight = {"x": [9.9, 10.0, 10.1, 10.0, 9.95, 10.05]}
        loose = {"x": [5.0, 15.0, 8.0, 12.0, 6.0, 14.0]}
        assert np.mean(tight["x"]) == pytest.approx(np.mean(loose["x"]))
        loose_ci = calculate_confidence_intervals(loose, method="normal")["x"]
        tight_ci = calculate_confidence_intervals(tight, method="normal")["x"]
        assert loose_ci["ci_upper"] - loose_ci["ci_lower"] > 10 * (tight_ci["ci_upper"] - tight_ci["ci_lower"])

    def test_interval_narrows_as_observations_accumulate(self) -> None:
        rng = np.random.default_rng(5)
        small = {"x": rng.normal(5.0, 2.0, 10)}
        large = {"x": rng.normal(5.0, 2.0, 1000)}
        for method in ("bootstrap", "normal"):
            narrow = calculate_confidence_intervals(large, method=method)["x"]
            wide = calculate_confidence_intervals(small, method=method)["x"]
            assert narrow["ci_upper"] - narrow["ci_lower"] < 0.3 * (wide["ci_upper"] - wide["ci_lower"])

    @pytest.mark.parametrize("method", ["normal", "bootstrap"])
    def test_nominal_coverage_of_the_true_mean(self, method: str) -> None:
        """A 95% interval from the dictionary form must cover the true mean about 95% of the time."""
        rng = np.random.default_rng(2026)
        true_mean, replicates = 5.0, 300
        covered = 0
        for _ in range(replicates):
            interval = calculate_confidence_intervals({"x": rng.normal(true_mean, 2.0, 40)}, method=method)["x"]
            covered += interval["ci_lower"] <= true_mean <= interval["ci_upper"]
        assert stats.binomtest(int(covered), replicates, 0.95).pvalue > 1e-3
        assert covered / replicates > 0.90

    @pytest.mark.parametrize("method", ["bootstrap", "normal", "percentile"])
    def test_array_input_gives_the_list_interval(self, method: str) -> None:
        """A numpy array used to fail on the truth value of the array itself."""
        values = [1.0, 2.0, 3.5, 4.0, 5.5, 7.0]
        assert calculate_confidence_intervals(np.array(values), method=method) == calculate_confidence_intervals(
            values, method=method
        )

    def test_integer_array_and_list_agree_with_floats(self) -> None:
        assert calculate_confidence_intervals(np.array([1, 2, 3, 4, 5]), method="normal") == pytest.approx(
            calculate_confidence_intervals([1.0, 2.0, 3.0, 4.0, 5.0], method="normal")
        )

    def test_dictionary_intervals_are_reproducible_and_leave_numpy_state_alone(self) -> None:
        np.random.seed(77)
        before = np.random.get_state()
        first = calculate_confidence_intervals(self.OBSERVATIONS)
        second = calculate_confidence_intervals(self.OBSERVATIONS)
        assert first == second
        assert _same_numpy_state(before, np.random.get_state())

    def test_dictionary_entries_may_be_arrays(self) -> None:
        as_arrays = {name: np.array(values) for name, values in self.OBSERVATIONS.items()}
        assert calculate_confidence_intervals(as_arrays) == calculate_confidence_intervals(self.OBSERVATIONS)

    @pytest.mark.parametrize(
        "bad,message",
        [
            ({}, "empty"),
            ({"pi": []}, "empty"),
            ({"pi": [[1.0, 2.0], [3.0, 4.0]]}, "one-dimensional"),
        ],
    )
    def test_empty_or_malformed_dictionaries_raise(self, bad: dict, message: str) -> None:
        with pytest.raises(ValueError, match=message):
            calculate_confidence_intervals(bad)

    def test_method_and_level_are_validated_for_dictionaries_too(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            calculate_confidence_intervals(self.OBSERVATIONS, method="jackknife")
        with pytest.raises(ValueError, match="between 0 and 1"):
            calculate_confidence_intervals(self.OBSERVATIONS, confidence_level=1.5)

    def test_two_dimensional_array_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="one-dimensional"):
            calculate_confidence_intervals(np.ones((3, 3)))
