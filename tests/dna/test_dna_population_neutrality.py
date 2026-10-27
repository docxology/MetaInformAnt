"""Neutrality tests (Tajima's D, Fu and Li's D*/F*, normalized Fay and Wu's H) against independent oracles.

Zero-mocks: closed forms are evaluated in exact rational arithmetic from the
published formulas (Tajima 1989; Fu and Li 1993 for D*; Zeng et al. 2006), and
F*, which has no closed form here, is checked against a seeded simulation of
the standard neutral coalescent.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from metainformant.dna.population import analysis, core

HAND = ["AAAA", "AAAA", "AACT", "AAAT"]
# Six sequences: singletons, doubletons, a 3/3 split, a tri-allelic column and an N column.
SIX = ["AAAAAAAA", "AAAAACAA", "CAAAAAAA", "AGAAAAGA", "AGTAAGGA", "AATTNAGC"]


def _harmonic(n: int) -> tuple[Fraction, Fraction]:
    return sum(Fraction(1, i) for i in range(1, n)), sum(Fraction(1, i * i) for i in range(1, n))


def _published_coefficients(n: int) -> dict[str, tuple[Fraction, Fraction]]:
    """(alpha, beta) of Var(N) = alpha*theta + beta*theta^2 from the published closed forms."""
    a, b = _harmonic(n)
    b1 = Fraction(n + 1, 3 * (n - 1))
    b2 = Fraction(2 * (n * n + n + 3), 9 * n * (n - 1))
    tajima = (b1 - 1 / a, b2 - Fraction(n + 2) / (a * n) + b / a**2)
    c = Fraction(2) * (n * a - 2 * (n - 1)) / ((n - 1) * (n - 2))
    a_next = a + Fraction(1, n)
    d = (
        c
        + Fraction(n - 2, (n - 1) ** 2)
        + Fraction(2, n - 1) * (Fraction(3, 2) - (2 * a_next - 3) / (n - 2) - Fraction(1, n))
    )
    v = (Fraction(n, n - 1) ** 2 * b + a * a * d - 2 * n * a * (a + 1) / Fraction((n - 1) ** 2)) / (a * a + b)
    u = Fraction(n, n - 1) * (a - Fraction(n, n - 1)) - v
    b_next = b + Fraction(1, n * n)
    zeng_beta = (18 * n * n * (3 * n + 2) * b_next - (88 * n**3 + 9 * n * n - 13 * n + 6)) / Fraction(
        9 * n * (n - 1) ** 2
    )
    return {
        "tajima_d": tajima,
        "fu_li_d_star": (a * (u + v), v * (a * a + b)),
        "fay_wu_h": (Fraction(n - 2, 6 * (n - 1)), zeng_beta),
    }


def _standardized(numerator: Fraction, alpha: Fraction, beta: Fraction, n: int, mutations: int) -> float:
    a, b = _harmonic(n)
    variance = alpha * mutations / a + beta * mutations * (mutations - 1) / (a * a + b)
    return float(numerator) / math.sqrt(variance)


def test_tajima_d_matches_tajima_1989_on_a_hand_example() -> None:
    # k = 7/6 pairwise differences, S = 2, n = 4.
    alpha, beta = _published_coefficients(4)["tajima_d"]
    expected = _standardized(Fraction(7, 6) - 2 / _harmonic(4)[0], alpha, beta, 4, 2)

    assert expected == pytest.approx(0.5915801398995606, rel=1e-12)
    assert core.tajimas_d(HAND) == pytest.approx(expected, rel=1e-12)
    tajima_d, p_value = analysis.calculate_tajima_d(HAND)
    assert tajima_d == pytest.approx(expected, rel=1e-12)
    assert p_value == pytest.approx(math.erfc(expected / math.sqrt(2)), rel=1e-12)


@pytest.mark.parametrize("n", [4, 5, 6, 9, 10, 17, 40])
def test_variance_coefficients_reproduce_published_closed_forms(n: int) -> None:
    computed = core._neutrality_variance_coefficients(n)
    for test, (alpha, beta) in _published_coefficients(n).items():
        assert computed[test] == pytest.approx((float(alpha), float(beta)), rel=1e-10), test
    f_star_alpha = Fraction(n + 1, 3 * (n - 1)) + Fraction(n - 5, n)  # sum_i w_i^2 / i, derived by hand
    assert computed["fu_li_f_star"][0] == pytest.approx(float(f_star_alpha), rel=1e-12)


def _coalescent_branch_lengths(n: int, replicates: int, rng: np.random.Generator) -> np.ndarray:
    """Total branch length subtending i = 1..n-1 leaves; time in units of 2N generations."""
    sizes = np.ones((replicates, n), dtype=np.int64)
    lengths = np.zeros((replicates, n + 1))
    rows = np.arange(replicates)
    for k in range(n, 1, -1):
        epoch = rng.exponential(2.0 / (k * (k - 1)), size=replicates)
        for slot in range(k):
            np.add.at(lengths, (rows, sizes[:, slot]), epoch)
        first = rng.integers(0, k, size=replicates)
        second = rng.integers(0, k - 1, size=replicates)
        second += second >= first
        sizes[rows, first] += sizes[rows, second]
        sizes[rows, second] = sizes[:, k - 1]
    return lengths[:, 1:n]


def test_genealogical_variance_matches_coalescent_simulation() -> None:
    """beta = Var(sum_i w_i L_i) / 4, since mutations fall at rate theta/2 per unit branch length."""
    n = 8
    lengths = _coalescent_branch_lengths(n, 60_000, np.random.default_rng(20261026))
    assert lengths.mean(axis=0) * np.arange(1, n) == pytest.approx(np.full(n - 1, 2.0), rel=0.02)

    weights = core._neutrality_weights(n)
    for test, (_, beta) in core._neutrality_variance_coefficients(n).items():
        simulated = float(np.var(lengths @ weights[test]) / 4.0)
        assert simulated == pytest.approx(beta, rel=0.06), test


def test_statistics_on_six_sequences_match_exact_values() -> None:
    # Folded spectrum: five singleton mutations (two at the tri-allelic column), two doubletons
    # and one 3/3 split; the column holding an N is excluded.
    assert core._site_frequency_spectrum(SIX).tolist() == [0, 5, 2, 1, 0, 0]
    coefficients = _published_coefficients(6)
    expected_tajima = _standardized(_six_pairwise() - 8 / _harmonic(6)[0], *coefficients["tajima_d"], 6, 8)
    expected_d_star = _standardized(Fraction(6, 5) * 8 - _harmonic(6)[0] * 5, *coefficients["fu_li_d_star"], 6, 8)

    assert core.tajimas_d(SIX) == pytest.approx(expected_tajima, rel=1e-12)
    assert core.fu_and_li_d_star_from_sequences(SIX) == pytest.approx(expected_d_star, rel=1e-12)
    assert core.fu_and_li_f_star_from_sequences(SIX) == pytest.approx(-0.300133743312, rel=1e-10)
    assert expected_tajima == pytest.approx(-0.286578219237, rel=1e-10)


def _six_pairwise() -> Fraction:
    """Mean pairwise differences of SIX from its spectrum: sum_j j (n - j) eta_j / C(n, 2)."""
    return Fraction(1 * 5 * 5 + 2 * 4 * 2 + 3 * 3 * 1, 15)


def test_four_sequence_statistics_coincide() -> None:
    """With n = 4 the folded spectrum has two classes, so every zero-mean statistic agrees once standardized."""
    values = [
        core.tajimas_d(HAND),
        core.fu_and_li_d_star_from_sequences(HAND),
        core.fu_and_li_f_star_from_sequences(HAND),
    ]

    assert values == pytest.approx([values[0]] * 3, rel=1e-12)


def test_gapped_and_ambiguous_columns_are_excluded() -> None:
    padded = [seq + extra for seq, extra in zip(SIX, "AC-AGA")]
    ambiguous = [seq + extra for seq, extra in zip(SIX, "ACNAGA")]

    for statistic in (core.tajimas_d, core.fu_and_li_d_star_from_sequences, core.fu_and_li_f_star_from_sequences):
        assert statistic(padded) == pytest.approx(statistic(SIX), rel=1e-12)
        assert statistic(ambiguous) == pytest.approx(statistic(SIX), rel=1e-12)


def test_normalized_fay_wu_h_uses_the_outgroup_polarized_spectrum() -> None:
    outgroup = "AGAAAAAA"  # column 2's ancestral G makes its derived A a four-copy mutation
    assert core._site_frequency_spectrum(SIX, outgroup).tolist() == [0, 5, 1, 1, 1, 0]
    alpha, beta = _published_coefficients(6)["fay_wu_h"]
    theta_pi, theta_l = Fraction(10, 3), Fraction(14, 5)
    expected = _standardized(theta_pi - theta_l, alpha, beta, 6, 8)

    assert core.normalized_fay_wu_h(SIX, outgroup) == pytest.approx(expected, rel=1e-12)
    h, p_value = analysis.calculate_fay_wu_h(SIX, outgroup)
    assert h == pytest.approx(core.fay_wu_h_from_sequences(SIX, outgroup))
    assert p_value == pytest.approx(math.erfc(expected / math.sqrt(2)), rel=1e-12)
    assert analysis.calculate_fay_wu_h(SIX)[1] is None
    assert analysis.neutrality_test_suite(SIX, outgroup)["fay_wu_h_p_value"] == pytest.approx(p_value)
    with pytest.raises(ValueError, match="aligned"):
        core.normalized_fay_wu_h(SIX, outgroup[:-1])


def test_monomorphic_samples_give_zero_statistics_and_unit_p_values() -> None:
    same = ["ACGT"] * 5
    for calculate in (analysis.calculate_tajima_d, analysis.calculate_fu_li_d, analysis.calculate_fu_li_f):
        assert calculate(same) == (0.0, 1.0)
    for statistic in (core.tajimas_d, core.fu_and_li_d_star_from_sequences, core.fu_and_li_f_star_from_sequences):
        with pytest.raises(ValueError, match="at least 4"):
            statistic(same[:3])
        with pytest.raises(ValueError, match="aligned"):
            statistic(same[:4] + ["ACG"])


def test_genotype_heterozygosity_counts_only_dosage_one() -> None:
    """0/1/2 alternate-allele counts: 2 is homozygous, and values outside 0..2 are missing."""
    homozygous = analysis.calculate_summary_statistics(genotype_matrix=[[2, 0], [2, 0], [0, 2]])
    assert homozygous["mean_heterozygosity"] == 0.0

    mixed = analysis.calculate_summary_statistics(genotype_matrix=[[1, 2], [1, -1], [0, 1], [2, 1]])
    # Locus 1: two heterozygotes of four calls; locus 2: two of three calls.
    assert mixed["mean_heterozygosity"] == pytest.approx((2 / 4 + 2 / 3) / 2)
    assert mixed["heterozygosity_variance"] == pytest.approx(((2 / 4 - 7 / 12) ** 2 + (2 / 3 - 7 / 12) ** 2) / 2)
