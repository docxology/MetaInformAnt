"""Validation of the seeded coalescent simulator against established neutral-model results.

Zero-mocks. Every check compares seeded simulations with a closed form (Watterson 1975, Tajima 1989,
Hudson 1990) or with numerical integration, and every tolerance is four standard errors measured from
the replicates themselves. The simulations are seeded, so the checks are deterministic.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest
from scipy import integrate

from metainformant.core.utils.errors import ValidationError
from metainformant.dna.population import core
from metainformant.simulation.models.popgen import (
    PopulationEpoch,
    _SizeHistory,
    haplotypes_to_sequences,
    simulate_coalescent_haplotypes,
    simulate_coalescent_spectra,
)

Z = 4.0  # tolerance in measured standard errors


def _harmonic(n: int) -> tuple[float, float]:
    """a_n and b_n: sums of 1/i and 1/i^2 for i < n."""
    return sum(1 / i for i in range(1, n)), sum(1 / i**2 for i in range(1, n))


def _pairwise(spectra: np.ndarray, n: int) -> np.ndarray:
    """Mean number of pairwise differences of every replicate."""
    i = np.arange(1, n)
    return spectra @ (i * (n - i) / (n * (n - 1) / 2))


def _tajima_variance_terms(n: int) -> tuple[float, float]:
    """c1, c2 of Tajima (1989): Var(pi - S/a_n) = c1 theta + c2 theta^2."""
    a1, a2 = _harmonic(n)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    return b1 - 1 / a1, b2 - (n + 2) / (a1 * n) + a2 / a1**2


def _tajima_d(spectra: np.ndarray, n: int) -> np.ndarray:
    """Tajima's (1989) D of every replicate from the published constants; 0 without mutations."""
    a1, a2 = _harmonic(n)
    c1, c2 = _tajima_variance_terms(n)
    segregating = spectra.sum(axis=1)
    variance = (c1 / a1) * segregating + (c2 / (a1**2 + a2)) * segregating * (segregating - 1)
    numerator = _pairwise(spectra, n) - segregating / a1
    return np.divide(numerator, np.sqrt(variance), out=np.zeros_like(numerator), where=variance > 0)


def _assert_mean(values: np.ndarray, expected: float, label: str) -> None:
    mean, se = float(np.mean(values)), float(np.std(values, ddof=1) / math.sqrt(len(values)))
    assert abs(mean - expected) <= Z * se, f"{label}: mean {mean:.5f} vs {expected:.5f} (standard error {se:.5f})"


def _assert_proportion(successes: int, trials: int, expected: float, label: str) -> None:
    se = math.sqrt(expected * (1 - expected) / trials)
    assert abs(successes / trials - expected) <= Z * se, f"{label}: {successes / trials:.5f} vs {expected:.5f}"


# --- constant size: Watterson, Tajima and Fu results -------------------------------------------------------


def test_segregating_sites_match_wattersons_mean_and_variance() -> None:
    n, theta, replicates = 10, 5.0, 40_000
    sites = simulate_coalescent_spectra(n, replicates, theta=theta, rng=101).sum(axis=1)
    a, b = _harmonic(n)

    _assert_mean(sites, theta * a, "E[S] = theta a_n")
    # Var(S) = theta a_n + theta^2 b_n; the standard error of a sample variance follows from the fourth moment.
    squares = (sites - sites.mean()) ** 2
    se = squares.std(ddof=1) / math.sqrt(replicates)
    assert abs(squares.mean() - (theta * a + theta**2 * b)) <= Z * se


def test_mean_pairwise_difference_is_theta() -> None:
    n, theta = 12, 7.0
    spectra = simulate_coalescent_spectra(n, 40_000, theta=theta, rng=102)
    _assert_mean(_pairwise(spectra, n), theta, "E[pi] = theta")


@pytest.mark.parametrize("n, theta", [(6, 3.0), (10, 5.0), (25, 12.0)])
def test_tajima_numerator_has_its_exact_mean_and_variance(n: int, theta: float) -> None:
    """pi - S/a_n has mean 0 and variance c1 theta + c2 theta^2 (Tajima 1989), so its z-score is standard."""
    spectra = simulate_coalescent_spectra(n, 40_000, theta=theta, rng=103 + n)
    c1, c2 = _tajima_variance_terms(n)
    z = (_pairwise(spectra, n) - spectra.sum(axis=1) / _harmonic(n)[0]) / math.sqrt(c1 * theta + c2 * theta**2)

    _assert_mean(z, 0.0, "mean z")
    squares = z**2
    se = squares.std(ddof=1) / math.sqrt(z.size)
    assert abs(squares.mean() - 1.0) <= Z * se, f"variance of z is {squares.mean():.4f} (standard error {se:.4f})"


def test_site_frequency_spectrum_is_theta_over_i() -> None:
    """E[xi_i] = theta / i (Fu 1995): the number of mutations carried by i of n sequences."""
    n, theta, replicates = 10, 6.0, 40_000
    spectra = simulate_coalescent_spectra(n, replicates, theta=theta, rng=104)
    assert spectra.shape == (replicates, n - 1)
    for i in range(1, n):
        _assert_mean(spectra[:, i - 1], theta / i, f"E[xi_{i}] = theta / {i}")


def test_probability_of_no_segregating_site_is_an_exact_product() -> None:
    """P(S = 0) = prod_{j<n} j / (j + theta): each of the n - 1 coalescent intervals must carry no mutation."""
    n, theta, replicates = 5, 2.0, 60_000
    sites = simulate_coalescent_spectra(n, replicates, theta=theta, rng=105).sum(axis=1)
    expected = math.prod(j / (j + theta) for j in range(1, n))
    _assert_proportion(int(np.sum(sites == 0)), replicates, expected, "P(S = 0)")


def test_mean_tajimas_d_is_near_zero_under_constant_size() -> None:
    """The estimated-variance ratio is slightly biased (about -0.1 for 20 sequences), so 'near 0' is |mean| < 0.15."""
    n, theta, replicates = 20, 10.0, 20_000
    d = _tajima_d(simulate_coalescent_spectra(n, replicates, theta=theta, rng=106), n)
    mean, se = d.mean(), d.std(ddof=1) / math.sqrt(replicates)

    assert abs(mean) < 0.15, f"mean D {mean:.3f} (standard error {se:.3f})"
    assert 0.7 < d.std() < 1.1  # the null spread of D is a little below 1


def test_the_old_star_genealogy_signature_is_not_reproduced() -> None:
    """Independent mutations from one ancestor give D near -2.6 for 20 sequences; the coalescent does not."""
    n, theta, replicates = 20, 10.0, 5_000
    d = _tajima_d(simulate_coalescent_spectra(n, replicates, theta=theta, rng=107), n)
    assert d.mean() - Z * d.std(ddof=1) / math.sqrt(replicates) > -0.5


def test_segregating_sites_conditioning_places_exactly_that_many_mutations() -> None:
    spectra = simulate_coalescent_spectra(8, 500, segregating_sites=13, rng=108)
    assert spectra.shape == (500, 7)
    assert (spectra.sum(axis=1) == 13).all()
    assert (spectra >= 0).all()
    assert simulate_coalescent_spectra(8, 5, segregating_sites=0, rng=1).sum() == 0
    assert simulate_coalescent_spectra(8, 5, theta=0.0, rng=1).sum() == 0


def _exact_class_probabilities_n4() -> list[float]:
    """P(a single mutation is carried by i = 1, 2, 3 of 4 sequences) when placed on a branch of an unweighted tree.

    With T_k ~ Exp(C(k, 2)) the total branch length is T = 4 T_4 + 3 T_3 + 2 T_2, and the branch lengths
    below i sequences are L_1 = 4 T_4 + 2 T_3 + T_2, L_2 = T_3, L_3 = T_2 for the topology in which the third
    merger joins a pair to a single lineage (probability 2/3), and L_1 = 4 T_4 + 2 T_3, L_2 = T_3 + 2 T_2,
    L_3 = 0 for two cherries (probability 1/3). The probability is E[L_i / T]; for independent exponentials
    X_j with rates r_j, E[X_j / sum X] = integral over s of r_j / (r_j + s)^2 times prod_{i != j} r_i / (r_i + s).
    """
    rates = [1.5, 1.0, 0.5]  # rates of 4 T_4, 3 T_3 and 2 T_2

    def share(j: int) -> float:
        def integrand(s: float) -> float:
            value = rates[j] / (rates[j] + s) ** 2
            for i, rate in enumerate(rates):
                if i != j:
                    value *= rate / (rate + s)
            return value

        return float(integrate.quad(integrand, 0.0, math.inf, epsabs=1e-13, epsrel=1e-13)[0])

    t3_share, t2_share = share(1), share(2)
    p3 = (2 / 3) * (t2_share / 2)
    p2 = (2 / 3) * (t3_share / 3) + (1 / 3) * (t3_share / 3 + t2_share)
    return [1 - p2 - p3, p2, p3]


def test_conditioned_mutation_is_placed_in_proportion_to_branch_length() -> None:
    """With S = 1 and n = 4 the class probabilities follow from integrating over the tree: (0.575, 0.262, 0.163).

    They are not the 1/i law (0.545, 0.273, 0.182 for i = 1, 2, 3 with the mutation weighted by the tree
    length), which tells apart conditioning on S from weighting every tree by its total length.
    """
    replicates = 200_000
    spectra = simulate_coalescent_spectra(4, replicates, segregating_sites=1, rng=109)
    expected = _exact_class_probabilities_n4()
    assert expected[1] == pytest.approx(0.26162, abs=2e-5)
    for i, probability in enumerate(expected):
        _assert_proportion(int(spectra[:, i].sum()), replicates, probability, f"P(carried by {i + 1} of 4)")
    assert abs(expected[1] - 3 / 11) > 8 * math.sqrt(0.27 * 0.73 / replicates)


# --- demography: rescaled coalescent rates -----------------------------------------------------------------


def test_recent_expansion_gives_negative_tajimas_d() -> None:
    n, theta, replicates = 20, 10.0, 20_000
    histories = {
        "instantaneous 100-fold expansion 0.5 units ago": [PopulationEpoch(0.0, 100.0), PopulationEpoch(0.5, 1.0)],
        "exponential 100-fold expansion over 0.5 units": [
            PopulationEpoch(0.0, 100.0, math.log(100.0) / 0.5),
            PopulationEpoch(0.5, 1.0),
        ],
    }
    for label, epochs in histories.items():
        d = _tajima_d(simulate_coalescent_spectra(n, replicates, theta=theta, epochs=epochs, rng=110), n)
        mean, se = d.mean(), d.std(ddof=1) / math.sqrt(replicates)
        assert mean + Z * se < -1.0, f"{label}: mean D {mean:.3f} (standard error {se:.3f})"


def test_bottleneck_observed_shortly_afterwards_gives_positive_d_and_fades() -> None:
    """A moderate bottleneck (size 0.1 for 0.05 units, about half a coalescent unit) leaves D > 0 at first.

    The excess of intermediate-frequency variants is then diluted by new mutations on the long external
    branches, so D falls through 0 and turns negative as the population recovers.
    """
    n, theta, replicates = 20, 20.0, 20_000

    def mean_d(recovery: float) -> tuple[float, float]:
        epochs = [PopulationEpoch(0.0, 1.0), PopulationEpoch(recovery, 0.1), PopulationEpoch(recovery + 0.05, 1.0)]
        d = _tajima_d(simulate_coalescent_spectra(n, replicates, theta=theta, epochs=epochs, rng=111), n)
        return float(d.mean()), float(d.std(ddof=1) / math.sqrt(replicates))

    soon, soon_se = mean_d(0.005)
    later, later_se = mean_d(0.2)
    assert soon - Z * soon_se > 0.3, f"mean D {soon:.3f} (standard error {soon_se:.3f})"
    assert later + Z * later_se < 0.0, f"mean D {later:.3f} (standard error {later_se:.3f})"


def _hazard(epochs: list[PopulationEpoch], t: float) -> float:
    """Integrated coalescent hazard int_0^t ds / size(s), by quadrature (independent of the simulator)."""
    total = 0.0
    ends = [epoch.start for epoch in epochs[1:]] + [math.inf]
    for epoch, end in zip(epochs, ends):
        if t <= epoch.start:
            break
        total += integrate.quad(
            lambda s, e=epoch: math.exp(e.growth_rate * (s - e.start)) / e.size, epoch.start, min(t, end)
        )[0]
    return total


def _pair_coalescence_moments(epochs: list[PopulationEpoch]) -> tuple[float, float]:
    """E[T] and E[T^2] of the coalescence time of two lineages, from P(T > t) = exp(-hazard(t))."""
    edges = [epoch.start for epoch in epochs] + [math.inf]
    first = second = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        first += integrate.quad(lambda t: math.exp(-_hazard(epochs, t)), low, high, limit=200)[0]
        second += integrate.quad(lambda t: 2 * t * math.exp(-_hazard(epochs, t)), low, high, limit=200)[0]
    return first, second


DEMOGRAPHIES = {
    "bottleneck": [PopulationEpoch(0.0, 1.0), PopulationEpoch(0.1, 0.05), PopulationEpoch(0.4, 1.0)],
    "exponential expansion": [PopulationEpoch(0.0, 50.0, math.log(50.0) / 0.5), PopulationEpoch(0.5, 1.0)],
    "growth then decline then growth": [
        PopulationEpoch(0.0, 4.0, 3.0),
        PopulationEpoch(0.3, 0.2, -2.0),
        PopulationEpoch(0.8, 1.5, 0.5),
        PopulationEpoch(1.5, 1.0),
    ],
}


@pytest.mark.parametrize("name", sorted(DEMOGRAPHIES))
def test_pairwise_coalescence_time_follows_the_integrated_hazard(name: str) -> None:
    """With n = 2 and mutation parameter theta, S given the coalescence time T is Poisson(theta T).

    Hence E[S] = theta E[T] and E[S(S - 1)] = theta^2 E[T^2], where E[T] and E[T^2] follow from
    P(T > t) = exp(-integral of 1 / size), evaluated here by numerical quadrature.
    """
    epochs, theta, replicates = DEMOGRAPHIES[name], 3.0, 200_000
    mean_time, second_moment = _pair_coalescence_moments(epochs)
    sites = simulate_coalescent_spectra(2, replicates, theta=theta, epochs=epochs, rng=112)[:, 0].astype(float)

    _assert_mean(sites, theta * mean_time, f"{name}: E[S] = theta E[T]")
    _assert_mean(sites * (sites - 1), theta**2 * second_moment, f"{name}: E[S(S-1)] = theta^2 E[T^2]")


@pytest.mark.parametrize("name", sorted(DEMOGRAPHIES))
def test_the_hazard_inverse_reaches_every_epoch_exactly(name: str) -> None:
    """Mapping an integrated hazard back to time lands where quadrature says the hazard equals it."""
    epochs = DEMOGRAPHIES[name]
    history = _SizeHistory(epochs)
    boundaries = [_hazard(epochs, epoch.start) for epoch in epochs[1:]]
    targets = np.array(
        sorted([0.001, 0.3, 1.0, 4.0, 25.0, 100.0] + [b * f for b in boundaries for f in (0.5, 1.0, 1.5)])
    )
    times = history.time_at_hazard(targets)

    assert (np.diff(times) > 0).all()
    for target, time in zip(targets, times):
        assert _hazard(epochs, float(time)) == pytest.approx(float(target), rel=1e-8, abs=1e-10)


def test_haplotype_and_spectrum_simulators_agree_under_a_bottleneck() -> None:
    """The topology-tracking haplotype simulator and the vectorized spectrum simulator sample one model."""
    n, theta, replicates = 12, 6.0, 5_000
    epochs = [PopulationEpoch(0.0, 1.0), PopulationEpoch(0.02, 0.1), PopulationEpoch(0.1, 1.0)]
    generator = np.random.default_rng(113)
    from_haplotypes = np.zeros((replicates, n - 1))
    for r in range(replicates):
        carriers = simulate_coalescent_haplotypes(n, theta, epochs=epochs, rng=generator).sum(axis=0, dtype=np.int64)
        from_haplotypes[r] = np.bincount(carriers, minlength=n)[1:n]
    from_spectra = simulate_coalescent_spectra(n, replicates, theta=theta, epochs=epochs, rng=114)

    for i in range(n - 1):
        difference = from_haplotypes[:, i].mean() - from_spectra[:, i].mean()
        se = math.sqrt(from_haplotypes[:, i].var(ddof=1) / replicates + from_spectra[:, i].var(ddof=1) / replicates)
        assert abs(difference) <= Z * se, f"class {i + 1}: difference {difference:.4f} (standard error {se:.4f})"
    sites_difference = from_haplotypes.sum(axis=1).mean() - from_spectra.sum(axis=1).mean()
    sites_se = math.sqrt(
        from_haplotypes.sum(axis=1).var(ddof=1) / replicates + from_spectra.sum(axis=1).var(ddof=1) / replicates
    )
    assert abs(sites_difference) <= Z * sites_se


def test_large_requests_are_simulated_in_blocks_without_changing_the_model() -> None:
    n, theta, replicates = 20, 4.0, 110_000  # more than 2 million replicate-samples: two blocks
    spectra = simulate_coalescent_spectra(n, replicates, theta=theta, rng=115)
    assert spectra.shape == (replicates, n - 1)
    _assert_mean(spectra.sum(axis=1), theta * _harmonic(n)[0], "E[S] over blocks")


# --- haplotype matrices and sequences ----------------------------------------------------------------------


def test_haplotype_matrix_is_a_valid_infinite_sites_sample() -> None:
    n = 15
    for seed in range(20):
        haplotypes = simulate_coalescent_haplotypes(n, 8.0, rng=random.Random(seed))
        sites = haplotypes.shape[1]
        assert haplotypes.shape[0] == n and haplotypes.dtype == np.uint8
        assert set(np.unique(haplotypes)) <= {0, 1}
        carriers = haplotypes.sum(axis=0)
        assert ((carriers >= 1) & (carriers <= n - 1)).all()  # every column is polymorphic
        # One genealogy and no recurrent mutation: any two sites show at most three of the four gametes.
        for a in range(sites):
            for b in range(a + 1, sites):
                assert len(set(2 * haplotypes[:, a] + haplotypes[:, b])) <= 3


def test_haplotypes_are_reproducible_for_every_kind_of_seed() -> None:
    for make in (lambda: 7, lambda: random.Random(7), lambda: np.random.default_rng(7)):
        assert np.array_equal(
            simulate_coalescent_haplotypes(9, 5.0, rng=make()), simulate_coalescent_haplotypes(9, 5.0, rng=make())
        )
        assert np.array_equal(
            simulate_coalescent_spectra(9, 4, theta=5.0, rng=make()),
            simulate_coalescent_spectra(9, 4, theta=5.0, rng=make()),
        )
    assert not np.array_equal(
        simulate_coalescent_spectra(9, 50, theta=5.0, rng=1), simulate_coalescent_spectra(9, 50, theta=5.0, rng=2)
    )
    assert simulate_coalescent_haplotypes(9, 0.0, rng=1).shape == (9, 0)


def test_sequences_are_the_haplotype_matrix_written_at_distinct_positions() -> None:
    reference = "ACGT" * 100
    haplotypes = simulate_coalescent_haplotypes(12, 6.0, rng=random.Random(3))
    sequences = haplotypes_to_sequences(haplotypes, reference_sequence=reference, rng=random.Random(4))

    assert len(sequences) == 12 and {len(s) for s in sequences} == {len(reference)}
    assert core.segregating_sites(sequences) == haplotypes.shape[1]  # a distinct position for every column
    changed = [i for i in range(len(reference)) if any(seq[i] != reference[i] for seq in sequences)]
    patterns = sorted(tuple(int(seq[i] != reference[i]) for seq in sequences) for i in changed)
    assert patterns == sorted(tuple(column) for column in haplotypes.T.tolist())  # columns map one to one onto sites
    for i in changed:  # each site has one derived base, and it is not the ancestral one
        assert len({seq[i] for seq in sequences} - {reference[i]}) == 1


def test_random_ancestral_sequences_are_seeded_and_use_all_four_bases() -> None:
    haplotypes = simulate_coalescent_haplotypes(8, 3.0, rng=random.Random(11))
    first = haplotypes_to_sequences(haplotypes, 300, rng=random.Random(12))
    assert first == haplotypes_to_sequences(haplotypes, 300, rng=random.Random(12))
    assert first != haplotypes_to_sequences(haplotypes, 300, rng=random.Random(13))
    assert set("".join(first)) == set("ACGT") and {len(s) for s in first} == {300}


def test_sequences_keep_the_reference_and_never_mutate_ambiguous_positions() -> None:
    reference = ("ACGTN" * 80) + "acgtn" * 20  # 500 positions, 100 of them N or n
    haplotypes = simulate_coalescent_haplotypes(10, 12.0, rng=random.Random(5))
    sequences = haplotypes_to_sequences(haplotypes, reference_sequence=reference, rng=random.Random(6))

    columns = list(zip(*sequences))
    variable = [i for i, column in enumerate(columns) if len(set(column)) > 1]
    assert len(variable) == haplotypes.shape[1]
    for i, column in enumerate(columns):
        if i not in variable:
            assert set(column) == {reference[i]}
        else:
            assert reference[i] not in "Nn"
            assert reference[i] in column  # the ancestral allele is still carried by some sequence
            derived = (set(column) - {reference[i]}).pop()
            assert derived != reference[i] and derived.isupper() == reference[i].isupper()


def test_too_many_mutations_for_the_sequence_are_rejected() -> None:
    haplotypes = simulate_coalescent_haplotypes(10, 40.0, rng=random.Random(8))
    assert haplotypes.shape[1] > 20
    with pytest.raises(ValueError, match="infinite-sites"):
        haplotypes_to_sequences(haplotypes, 20, rng=random.Random(9))


# --- input validation --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "epochs, message",
    [
        ([PopulationEpoch(0.1, 1.0)], "start at time 0"),
        ([PopulationEpoch(0.0, 1.0), PopulationEpoch(0.0, 2.0)], "strictly increasing"),
        ([PopulationEpoch(0.0, 1.0), PopulationEpoch(0.5, 0.0)], "positive"),
        ([PopulationEpoch(0.0, 1.0), PopulationEpoch(0.5, 1.0, -1.0)], "growth_rate >= 0"),
        ([PopulationEpoch(0.0, 1.0), PopulationEpoch(float("inf"), 1.0)], "finite"),
        ([PopulationEpoch(0.0, 1.0, 1e5), PopulationEpoch(0.5, 1.0)], "overflows"),
    ],
)
def test_invalid_demographic_histories_are_rejected(epochs: list[PopulationEpoch], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        simulate_coalescent_spectra(5, 3, theta=1.0, epochs=epochs, rng=1)


def test_invalid_simulation_requests_are_rejected() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        simulate_coalescent_spectra(5, 3, rng=1)
    with pytest.raises(ValueError, match="exactly one"):
        simulate_coalescent_spectra(5, 3, theta=1.0, segregating_sites=2, rng=1)
    with pytest.raises(ValueError, match="non-negative integer"):
        simulate_coalescent_spectra(5, 3, segregating_sites=-1, rng=1)
    with pytest.raises(ValueError, match="non-negative integer"):
        simulate_coalescent_spectra(5, 3, segregating_sites=2.5, rng=1)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="finite"):
        simulate_coalescent_haplotypes(5, float("nan"), rng=1)
    with pytest.raises(ValidationError, match="n_samples"):
        simulate_coalescent_haplotypes(1, 1.0, rng=1)
    with pytest.raises(ValueError, match="length"):
        haplotypes_to_sequences(np.zeros((3, 2), dtype=np.uint8), rng=1)
    with pytest.raises(ValueError, match="2-D"):
        haplotypes_to_sequences(np.zeros(4, dtype=np.uint8), 10, rng=1)
    with pytest.raises(ValueError, match="only 0"):
        haplotypes_to_sequences(np.full((3, 2), 2, dtype=np.uint8), 10, rng=1)
    with pytest.raises(ValueError, match="does not match"):
        haplotypes_to_sequences(np.zeros((3, 2), dtype=np.uint8), 5, reference_sequence="ACGTACGT", rng=1)
