"""Simulated null distributions for the neutrality tests (Tajima's D, Fu and Li's D*/F*, Fay and Wu's H).

Zero-mocks. The p-values of ``p_value_method="simulation"`` are checked on samples drawn from the same neutral
coalescent (uniform when the sample is neutral, valid when it holds few mutations) and on samples from strong
expansions and bottlenecks (small), and the Monte Carlo arithmetic and the conditioning on the number of mutations
are checked against hand-computed and exactly integrated values.
"""

from __future__ import annotations

import math
import random
import subprocess
import sys
from typing import Any, Callable, Tuple

import numpy as np
import pytest
from scipy import integrate, stats

from metainformant.dna.integration import neutrality_null
from metainformant.dna.population import analysis, core
from metainformant.simulation.models.popgen import (
    generate_population_sequences,
    simulate_bottleneck_population,
    simulate_population_expansion,
)

Z = 4.0  # tolerance in standard errors
TESTS: dict[str, Callable[..., Tuple[float, float]]] = {
    "tajima_d": analysis.calculate_tajima_d,
    "fu_li_d_star": analysis.calculate_fu_li_d,
    "fu_li_f_star": analysis.calculate_fu_li_f,
}


def _neutral_p_values(test: str, theta_per_site: float, samples: int, simulations: int) -> np.ndarray:
    """Simulated p-values of ``samples`` neutral 20-sequence samples drawn from the coalescent."""
    p_values = []
    for index in range(samples):
        sequences = generate_population_sequences(20, 600, theta=theta_per_site, rng=random.Random(1000 + index))
        _, p_value = TESTS[test](sequences, p_value_method="simulation", n_simulations=simulations, random_state=index)
        p_values.append(p_value)
    return np.array(p_values)


# --- Monte Carlo arithmetic --------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "observed, expected",
    [
        (3.0, 0.8),  # lower tail (1 + 4) / 5 = 1, upper tail (1 + 1) / 5 = 0.4: twice the smaller
        (10.0, 0.4),  # beyond every simulated value: the floor 2 / (R + 1)
        (-10.0, 0.4),
        (1.5, 1.0),  # both tails 0.6: twice 0.6 is capped at 1
        (0.0, 0.8),
    ],
)
def test_monte_carlo_p_value_uses_the_plus_one_correction_and_the_smaller_tail(
    observed: float, expected: float
) -> None:
    simulated = np.array([0.0, 1.0, 2.0, 3.0])
    assert neutrality_null.monte_carlo_two_sided_p(observed, simulated) == pytest.approx(expected)


def test_monte_carlo_ties_count_in_both_tails_up_to_rounding() -> None:
    simulated = np.array([1.0, 1.0 + 1e-13, 1.0 - 1e-13, 5.0, -5.0])
    # Three simulated values equal the observed one: lower (1 + 4) / 6, upper (1 + 4) / 6.
    assert neutrality_null.monte_carlo_two_sided_p(1.0, simulated) == 1.0
    assert neutrality_null.monte_carlo_two_sided_p(4.0, simulated) == pytest.approx(2 * 2 / 6)


# --- calibration -------------------------------------------------------------------------------------------


@pytest.mark.parametrize("test", sorted(TESTS))
def test_simulated_p_values_are_uniform_for_neutral_samples(test: str) -> None:
    """About 64 mutations per sample: the p-values of neutral samples are uniform on (0, 1).

    Measured over 2000 samples (not run here): P(p <= 0.05) is 0.047, 0.043 and 0.050 for D, D* and F*, against
    0.024, 0.057 and 0.050 for the normal approximation, which also gives no D below 0.01.
    """
    samples = 250
    p_values = _neutral_p_values(test, 0.03, samples, 400)

    assert stats.kstest(p_values, "uniform").statistic < 1.95 / math.sqrt(samples)  # 0.1 % critical value
    assert abs(p_values.mean() - 0.5) <= Z * math.sqrt(1 / 12 / samples)
    assert np.mean(p_values <= 0.05) <= 0.05 + Z * math.sqrt(0.05 * 0.95 / samples)


@pytest.mark.parametrize("test", sorted(TESTS))
def test_simulated_p_values_stay_valid_with_few_mutations(test: str) -> None:
    """About 10 mutations: the null is discrete, so the test is conservative but never rejects too often."""
    samples = 250
    p_values = _neutral_p_values(test, 0.005, samples, 400)
    for level in (0.05, 0.1):
        assert np.mean(p_values <= level) <= level + Z * math.sqrt(level * (1 - level) / samples)


def test_simulated_fay_wu_h_p_values_are_uniform_for_perfectly_polarized_neutral_samples() -> None:
    samples, p_values = 250, []
    for index in range(samples):
        ancestor = "".join(random.Random(900 + index).choices("ACGT", k=600))
        sequences = generate_population_sequences(
            20, reference_sequence=ancestor, theta=0.03, rng=random.Random(2000 + index)
        )
        _, p_value = analysis.calculate_fay_wu_h(
            sequences, ancestor, p_value_method="simulation", n_simulations=400, random_state=index
        )
        assert p_value is not None
        p_values.append(p_value)

    assert stats.kstest(p_values, "uniform").statistic < 1.95 / math.sqrt(samples)
    assert np.mean(np.array(p_values) <= 0.05) <= 0.05 + Z * math.sqrt(0.05 * 0.95 / samples)


# --- power -------------------------------------------------------------------------------------------------


def _expansion(rng: random.Random) -> list[str]:
    return simulate_population_expansion(
        n_sequences=40, sequence_length=2000, initial_size=100, expansion_factor=100.0, expansion_time=400, rng=rng
    )


def _bottleneck(rng: random.Random) -> list[str]:
    return simulate_bottleneck_population(
        n_sequences=40, sequence_length=3000, bottleneck_size=10, bottleneck_duration=10, rng=rng
    )


@pytest.mark.parametrize(
    "make, tests, sign",
    [
        (_expansion, sorted(TESTS), -1),
        (_bottleneck, ["fu_li_d_star", "fu_li_f_star"], +1),
    ],
)
def test_strongly_non_neutral_samples_have_small_simulated_p_values(
    make: Callable[[random.Random], list[str]], tests: list[str], sign: int
) -> None:
    """Median simulated p over 16 samples: 0.005-0.018 after expansions, 0.002-0.006 (D*, F*) after bottlenecks."""
    samples = [make(random.Random(500 + index)) for index in range(16)]
    for test in tests:
        results = [
            TESTS[test](sequences, p_value_method="simulation", n_simulations=1000, random_state=index)
            for index, sequences in enumerate(samples)
        ]
        statistics_, p_values = np.array([r[0] for r in results]), np.array([r[1] for r in results])
        assert np.mean(np.sign(statistics_) == sign) >= 0.8, test  # one locus is one random genealogy: not every sample
        assert np.median(p_values) < 0.02, test
        assert np.mean(p_values < 0.05) >= 0.7, test


# --- the options -------------------------------------------------------------------------------------------


def test_the_seed_makes_simulated_p_values_reproducible() -> None:
    sequences = generate_population_sequences(15, 400, theta=0.03, rng=random.Random(4))
    first = analysis.calculate_tajima_d(sequences, p_value_method="simulation", n_simulations=500, random_state=7)
    assert first == analysis.calculate_tajima_d(
        sequences, p_value_method="simulation", n_simulations=500, random_state=7
    )
    for make in (lambda: random.Random(7), lambda: np.random.default_rng(7)):
        assert analysis.calculate_tajima_d(
            sequences, p_value_method="simulation", n_simulations=500, random_state=make()
        ) == analysis.calculate_tajima_d(sequences, p_value_method="simulation", n_simulations=500, random_state=make())
    other = analysis.calculate_tajima_d(sequences, p_value_method="simulation", n_simulations=500, random_state=8)
    assert other[0] == first[0]  # the statistic does not depend on the simulation
    assert other[1] != first[1] and abs(other[1] - first[1]) < 0.15  # the p-value moves only by Monte Carlo error


def test_default_method_is_the_normal_approximation() -> None:
    sequences = generate_population_sequences(15, 400, theta=0.03, rng=random.Random(4))
    for name, function in TESTS.items():
        statistic, p_value = function(sequences)
        assert p_value == pytest.approx(math.erfc(abs(statistic) / math.sqrt(2)), rel=1e-12), name
        assert function(sequences, p_value_method="normal") == (statistic, p_value)
    suite = analysis.neutrality_test_suite(sequences)
    assert "p_value_method" not in suite and "n_simulations" not in suite
    assert suite["tajima_d_p_value"] == pytest.approx(math.erfc(abs(suite["tajima_d"]) / math.sqrt(2)), rel=1e-12)


def test_the_suite_passes_the_options_to_every_test() -> None:
    ancestor = "".join(random.Random(1).choices("ACGT", k=400))
    sequences = generate_population_sequences(12, reference_sequence=ancestor, theta=0.03, rng=random.Random(5))
    options = {"p_value_method": "simulation", "n_simulations": 500, "random_state": 7}
    suite = analysis.neutrality_test_suite(sequences, ancestor, **options)

    assert suite["p_value_method"] == "simulation" and suite["n_simulations"] == 500
    assert suite["tajima_d_p_value"] == analysis.calculate_tajima_d(sequences, **options)[1]
    assert suite["fu_li_d_star_p_value"] == analysis.calculate_fu_li_d(sequences, **options)[1]
    assert suite["fu_li_f_star_p_value"] == analysis.calculate_fu_li_f(sequences, **options)[1]
    h, h_p_value = analysis.calculate_fay_wu_h(sequences, ancestor, **options)
    assert (suite["fay_wu_h"], suite["fay_wu_h_p_value"]) == (h, h_p_value)
    normal = analysis.neutrality_test_suite(sequences, ancestor)
    assert suite["tajima_d"] == normal["tajima_d"] and suite["fay_wu_h"] == normal["fay_wu_h"]
    assert suite["tajima_d_p_value"] != normal["tajima_d_p_value"]
    assert set(suite["neutrality_summary"]) == set(normal["neutrality_summary"])
    # Without an outgroup the sites cannot be polarized, so there is no Fay and Wu p-value by either method.
    assert analysis.neutrality_test_suite(sequences, **options)["fay_wu_h_p_value"] is None


def test_the_options_are_validated_before_any_test_runs() -> None:
    sequences = generate_population_sequences(8, 100, theta=0.03, rng=random.Random(2))
    for function in (*TESTS.values(), analysis.calculate_fay_wu_h, analysis.neutrality_test_suite):
        with pytest.raises(ValueError, match="p_value_method"):
            function(sequences, p_value_method="bootstrap")
        for bad in (0, -5, 2.5, True):
            with pytest.raises(ValueError, match="n_simulations"):
                function(sequences, p_value_method="simulation", n_simulations=bad)
    with pytest.raises(ValueError, match="p_value_method"):
        analysis.neutrality_test_suite([], p_value_method="bootstrap")  # not swallowed by the per-test handlers
    with pytest.raises(ValueError, match="aligned"):
        analysis.calculate_tajima_d(["ACGT", "ACG", "ACGT", "ACGT"], p_value_method="simulation")


def test_trivial_samples_give_statistic_zero_and_p_value_one_by_either_method() -> None:
    same = ["ACGTACGT"] * 6
    for function in TESTS.values():
        for method in ("normal", "simulation"):
            assert function(same, p_value_method=method) == (0.0, 1.0)
            assert function(same[:3], p_value_method=method) == (0.0, 1.0)  # fewer than 4 sequences
    assert analysis.calculate_fay_wu_h(same, same[0], p_value_method="simulation") == (0.0, 1.0)


def test_p_values_cannot_beat_the_monte_carlo_resolution() -> None:
    sequences = _expansion(random.Random(3))
    for n_simulations in (9, 99):
        _, p_value = analysis.calculate_tajima_d(
            sequences, p_value_method="simulation", n_simulations=n_simulations, random_state=1
        )
        assert p_value >= 2 / (n_simulations + 1) - 1e-12


# --- conditioning on the number of mutations: an exact small case ------------------------------------------


def _doubleton_probability_for_four_sequences() -> float:
    """P(a single mutation is carried by 2 of 4 sequences) on an unweighted four-sequence coalescent tree.

    See ``tests/simulation/test_simulation_coalescent.py``: with X_j independent exponentials of rates
    (1.5, 1, 0.5), the probability is E[X_2 / sum X] / 3 + E[X_3 / sum X] / 3 (the second topology) and
    E[X_2 / sum X] / 3 * 2/3 (the first); each E[X_j / sum X] is a one-dimensional integral.
    """
    rates = [1.5, 1.0, 0.5]

    def share(j: int) -> float:
        def integrand(s: float) -> float:
            value = rates[j] / (rates[j] + s) ** 2
            for i, rate in enumerate(rates):
                if i != j:
                    value *= rate / (rate + s)
            return value

        return float(integrate.quad(integrand, 0.0, math.inf, epsabs=1e-13, epsrel=1e-13)[0])

    return (2 / 3) * (share(1) / 3) + (1 / 3) * (share(1) / 3 + share(2))


def test_four_sequence_p_values_match_the_exactly_integrated_conditional_null() -> None:
    """One mutation: a doubleton has p = 2 P(doubleton) = 0.5232 under the conditioned null (0.5455 for 1/i)."""
    doubleton = ["ACGT", "ACGT", "ACGA", "ACGA"]
    singleton = ["ACGT", "ACGT", "ACGT", "ACGA"]
    simulations = 200_000
    expected = 2 * _doubleton_probability_for_four_sequences()
    assert expected == pytest.approx(0.52325, abs=2e-5)

    _, p_value = analysis.calculate_tajima_d(
        doubleton, p_value_method="simulation", n_simulations=simulations, random_state=3
    )
    se = 2 * math.sqrt(expected / 2 * (1 - expected / 2) / simulations)
    assert abs(p_value - expected) <= Z * se
    assert abs(2 * 3 / 11 - expected) > 8 * se  # the 1/i law would be told apart
    # A singleton is the least extreme class on the upper side: every simulated value is at least as large.
    assert (
        analysis.calculate_tajima_d(singleton, p_value_method="simulation", n_simulations=1000, random_state=3)[1]
        == 1.0
    )
    # With n = 4 the standardized statistics coincide, so the three tests share their p-values.
    assert (
        analysis.calculate_fu_li_d(doubleton, p_value_method="simulation", n_simulations=1000, random_state=3)[1]
        == analysis.calculate_tajima_d(doubleton, p_value_method="simulation", n_simulations=1000, random_state=3)[1]
    )


def test_the_observed_and_simulated_statistics_share_one_computation() -> None:
    """The Monte Carlo comparison uses the statistic that core returns, to rounding."""
    sequences = generate_population_sequences(12, 300, theta=0.04, rng=random.Random(6))
    spectrum = core._site_frequency_spectrum(sequences)
    for test, function in (("tajima_d", core.tajimas_d), ("fu_li_d_star", core.fu_and_li_d_star_from_sequences)):
        from_spectrum = core._neutrality_statistics_from_spectra(spectrum[1:], 12, test)
        assert from_spectrum.shape == (1,)
        assert from_spectrum[0] == pytest.approx(function(sequences), rel=1e-12)
    several = core._neutrality_statistics_from_spectra(np.vstack([spectrum[1:], np.zeros(11)]), 12, "tajima_d")
    assert several[1] == 0.0 and several[0] == pytest.approx(core.tajimas_d(sequences), rel=1e-12)


# --- the integration adapter -------------------------------------------------------------------------------


def test_the_adapter_scores_each_group_of_tests_against_its_own_null() -> None:
    ancestor = "".join(random.Random(1).choices("ACGT", k=400))
    sequences = generate_population_sequences(12, reference_sequence=ancestor, theta=0.03, rng=random.Random(5))
    options: dict[str, Any] = {"n_simulations": 300, "random_state": 7}

    together = neutrality_null.simulated_neutrality_tests(
        sequences, ["tajima_d", "fay_wu_h", "fu_li_f_star"], outgroup=ancestor, **options
    )
    folded = neutrality_null.simulated_neutrality_tests(sequences, ["tajima_d", "fu_li_f_star"], **options)
    polarized = neutrality_null.simulated_neutrality_tests(sequences, ["fay_wu_h"], outgroup=ancestor, **options)

    assert set(together) == {"tajima_d", "fay_wu_h", "fu_li_f_star"}
    assert together == {**folded, **polarized}
    # The observed statistics are the ones core computes; the folded tests ignore the outgroup.
    assert together["tajima_d"][0] == pytest.approx(core.tajimas_d(sequences), rel=1e-12)
    assert together["fay_wu_h"][0] == pytest.approx(core.normalized_fay_wu_h(sequences, ancestor), rel=1e-12)


def test_the_adapter_rejects_unknown_tests_and_polarized_tests_without_an_outgroup() -> None:
    sequences = generate_population_sequences(8, 100, theta=0.03, rng=random.Random(2))
    with pytest.raises(ValueError, match="unknown neutrality test"):
        neutrality_null.simulated_neutrality_tests(sequences, ["tajima"], n_simulations=10)
    with pytest.raises(ValueError, match="outgroup"):
        neutrality_null.simulated_neutrality_tests(sequences, ["fay_wu_h"], n_simulations=10)


def test_the_simulation_domain_is_loaded_only_when_a_simulated_p_value_is_requested() -> None:
    """dna reaches simulation only through the integration adapter, which analysis imports lazily."""
    code = (
        "import sys\n"
        "from metainformant.dna.population import analysis\n"
        "seqs = ['ACGTACGT', 'ACGTACGA', 'ACGAACGT', 'ACGTTCGT', 'ACGTACCT']\n"
        "analysis.calculate_tajima_d(seqs)\n"
        "analysis.neutrality_test_suite(seqs)\n"
        "loaded = lambda: sum(name.startswith('metainformant.simulation') for name in sys.modules)\n"
        "before = loaded()\n"
        "analysis.calculate_tajima_d(seqs, p_value_method='simulation', n_simulations=20, random_state=1)\n"
        "print(before, loaded() > 0)\n"
    )
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert run.stdout.split() == ["0", "True"], run.stdout
