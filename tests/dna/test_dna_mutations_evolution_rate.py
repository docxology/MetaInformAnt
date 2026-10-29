"""Regression tests for ``simulate_sequence_evolution`` and the divergence/classification fixes in
``metainformant.dna.variation.mutations``.

REAL IMPLEMENTATION POLICY: the real functions run against independent statistical oracles (scipy's exact
binomial test and chi-square test, and hand-counted Hamming distances between consecutive generations);
nothing is mocked or patched.

``simulate_sequence_evolution`` applied exactly one point mutation per generation whatever ``mutation_rate``
was and had no way to seed its draws locally. It now draws the number of mutations each generation from
Binomial(length, rate) through its ``rng``.
"""

from __future__ import annotations

import math
import random

import pytest

from metainformant.dna.variation import mutations
from metainformant.dna.variation.mutations import (
    analyze_mutation_spectrum,
    classify_mutations,
    detect_selection_signatures,
    simulate_sequence_evolution,
)

stats = pytest.importorskip("scipy.stats")


def _hamming(first: str, second: str) -> int:
    return sum(1 for a, b in zip(first, second) if a != b)


def _substitutions_per_generation(lineage: list[str]) -> list[int]:
    """A generation changes ``k`` distinct positions to other bases, so consecutive sequences differ at ``k`` sites."""
    return [_hamming(lineage[i], lineage[i + 1]) for i in range(len(lineage) - 1)]


class TestMutationCountFollowsTheRate:
    def test_total_substitutions_match_rate_times_length_times_generations(self) -> None:
        """G independent Binomial(L, r) counts sum to Binomial(L G, r): mean r L G, s.e. sqrt(L G r (1 - r))."""
        length, generations, rate = 2000, 500, 0.002
        lineage = simulate_sequence_evolution("ACGT" * (length // 4), generations, rate, rng=random.Random(11))
        total = sum(_substitutions_per_generation(lineage))

        expected = rate * length * generations  # 2000
        standard_error = math.sqrt(length * generations * rate * (1 - rate))  # about 44.7
        assert abs(total - expected) < 4 * standard_error
        assert stats.binomtest(total, length * generations, rate).pvalue > 1e-3

    @pytest.mark.parametrize("rate", [0.0005, 0.01, 0.2])
    def test_mean_count_scales_with_the_rate(self, rate: float) -> None:
        length, generations = 1000, 400
        lineage = simulate_sequence_evolution("A" * length, generations, rate, rng=random.Random(5))
        counts = _substitutions_per_generation(lineage)
        standard_error = math.sqrt(length * rate * (1 - rate) / generations)  # of the mean count per generation
        assert abs(sum(counts) / generations - rate * length) < 4 * standard_error

    def test_per_generation_counts_are_binomially_distributed(self) -> None:
        """A fixed one-per-generation count (the old behaviour) fails this chi-square test of the binomial counts."""
        length, generations, rate = 50, 4000, 0.1
        lineage = simulate_sequence_evolution("ACGTA" * (length // 5), generations, rate, rng=random.Random(21))
        counts = _substitutions_per_generation(lineage)

        # Pool the tails so every expected cell holds at least five generations
        observed = [0] * 12
        for count in counts:
            observed[min(count, 11)] += 1
        probabilities = [stats.binom.pmf(k, length, rate) for k in range(11)]
        probabilities.append(1.0 - sum(probabilities))
        expected = [p * generations for p in probabilities]
        keep = [i for i, e in enumerate(expected) if e >= 5]
        observed_kept = [observed[i] for i in keep[:-1]] + [sum(observed[i] for i in range(keep[-1], 12))]
        expected_kept = [expected[i] for i in keep[:-1]] + [sum(expected[i] for i in range(keep[-1], 12))]
        assert stats.chisquare(observed_kept, expected_kept).pvalue > 1e-3

        # Binomial dispersion: variance L r (1 - r) = 4.5 (a constant count would have variance 0)
        mean = sum(counts) / generations
        variance = sum((c - mean) ** 2 for c in counts) / (generations - 1)
        assert variance == pytest.approx(length * rate * (1 - rate), rel=0.15)

    def test_the_count_varies_from_generation_to_generation(self) -> None:
        lineage = simulate_sequence_evolution("ACGT" * 25, 200, 0.05, rng=random.Random(3))
        assert len(set(_substitutions_per_generation(lineage))) > 5

    def test_rate_zero_changes_nothing(self) -> None:
        """Exactly one mutation per generation, whatever the rate, used to change the sequence even at rate 0."""
        lineage = simulate_sequence_evolution("ACGTACGT", 25, 0.0, rng=random.Random(1))
        assert lineage == ["ACGTACGT"] * 26

    def test_rate_one_changes_every_site_every_generation(self) -> None:
        sequence = "ACGTACGTAC"
        lineage = simulate_sequence_evolution(sequence, 10, 1.0, rng=random.Random(2))
        assert _substitutions_per_generation(lineage) == [len(sequence)] * 10

    def test_replacement_bases_always_differ_from_the_base_they_replace(self) -> None:
        lineage = simulate_sequence_evolution("A" * 300, 1, 0.5, rng=random.Random(8))
        replacements = [b for a, b in zip(lineage[0], lineage[1]) if a != b]
        assert replacements and set(replacements) <= {"C", "G", "T"}

    def test_long_sequences_with_small_rates_are_cheap_to_simulate(self) -> None:
        """The count is drawn by jumping between mutations: a megabase at rate 1e-5 costs about ten draws."""
        lineage = simulate_sequence_evolution("A" * 1_000_000, 5, 1e-5, rng=random.Random(4))
        counts = _substitutions_per_generation(lineage)
        assert all(2 <= count <= 25 for count in counts)  # Binomial(1e6, 1e-5): mean 10, standard deviation 3.16


class TestRandomSource:
    SEQUENCE = "ACGT" * 25

    def test_same_rng_seed_gives_the_same_lineage_whatever_the_global_state(self) -> None:
        random.seed(111)
        first = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05, rng=random.Random(42))
        random.seed(222)  # a different global state must not matter
        second = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05, rng=random.Random(42))
        assert first == second

    def test_different_seeds_give_different_lineages(self) -> None:
        first = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05, rng=random.Random(1))
        second = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05, rng=random.Random(2))
        assert first != second

    def test_local_generator_leaves_the_global_state_untouched(self) -> None:
        random.seed(5)
        before = random.getstate()
        simulate_sequence_evolution(self.SEQUENCE, 20, 0.05, rng=random.Random(1))
        assert random.getstate() == before

    def test_without_rng_the_global_state_still_controls_the_result(self) -> None:
        random.seed(7)
        first = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05)
        random.seed(7)
        second = simulate_sequence_evolution(self.SEQUENCE, 20, 0.05)
        assert first == second

    def test_every_draw_comes_from_the_given_generator(self) -> None:
        """Counts, positions and replacement bases all use ``rng``: the global state is neither read nor advanced."""
        random.seed(99)
        a = simulate_sequence_evolution(self.SEQUENCE, 10, 0.1, rng=random.Random(6))
        random.seed(100)
        b = simulate_sequence_evolution(self.SEQUENCE, 10, 0.1, rng=random.Random(6))
        assert a == b


class TestContract:
    @pytest.mark.parametrize("rate", [-0.1, 1.5, float("nan"), float("inf")])
    def test_rate_must_be_a_probability(self, rate: float) -> None:
        with pytest.raises(ValueError, match="mutation_rate"):
            simulate_sequence_evolution("ACGT", 3, rate)

    @pytest.mark.parametrize("sequence,generations", [("", 5), ("ACGT", 0), ("ACGT", -1)])
    def test_rate_is_checked_even_when_there_is_nothing_to_simulate(self, sequence: str, generations: int) -> None:
        with pytest.raises(ValueError, match="mutation_rate"):
            simulate_sequence_evolution(sequence, generations, 2.0)

    def test_trivial_lineages_are_unchanged(self) -> None:
        assert simulate_sequence_evolution("", 5) == [""]
        assert simulate_sequence_evolution("ACGT", 0) == ["ACGT"]
        assert simulate_sequence_evolution("ACGT", -3) == ["ACGT"]

    def test_first_entry_is_the_input_and_the_rest_are_upper_case(self) -> None:
        lineage = simulate_sequence_evolution("acgtacgt", 3, 0.2, rng=random.Random(9))
        assert lineage[0] == "acgtacgt"
        assert all(entry == entry.upper() and len(entry) == 8 for entry in lineage[1:])
        assert len(lineage) == 4


class TestBinomialDraw:
    def test_single_trial_is_a_bernoulli_draw(self) -> None:
        rng = random.Random(13)
        successes = sum(mutations._binomial_count(1, 0.3, rng) for _ in range(20000))
        assert stats.binomtest(successes, 20000, 0.3).pvalue > 1e-3

    def test_sparse_large_count_has_the_binomial_mean_and_variance(self) -> None:
        rng = random.Random(14)
        draws = [mutations._binomial_count(1_000_000, 1e-4, rng) for _ in range(600)]
        mean = sum(draws) / len(draws)
        variance = sum((d - mean) ** 2 for d in draws) / (len(draws) - 1)
        assert abs(mean - 100.0) < 4 * math.sqrt(100.0 * (1 - 1e-4) / len(draws))
        assert variance == pytest.approx(100.0 * (1 - 1e-4), rel=0.25)

    def test_counts_never_exceed_the_number_of_trials(self) -> None:
        rng = random.Random(15)
        assert max(mutations._binomial_count(7, 0.99, rng) for _ in range(500)) <= 7

    def test_degenerate_probabilities(self) -> None:
        rng = random.Random(16)
        assert mutations._binomial_count(50, 0.0, rng) == 0
        assert mutations._binomial_count(50, 1.0, rng) == 50
        assert mutations._binomial_count(0, 0.5, rng) == 0


class TestDivergenceAndClassification:
    """Placeholder and fabricated values in the other functions of the module."""

    def test_divergence_from_references_is_a_hand_counted_mean(self) -> None:
        # ATCGATCG vs ATCGATCA: 1 of 8 sites differs; vs TCGAATCG: 4 of 8 (A/T, T/C, C/G, G/A)
        result = detect_selection_signatures("ATCGATCG", ["ATCGATCA", "TCGAATCG"])
        assert result["diversity"] == pytest.approx((1 / 8 + 4 / 8) / 2)
        assert result["conservation"] == pytest.approx(1.0 - (1 / 8 + 4 / 8) / 2)

    def test_only_measured_quantities_are_reported(self) -> None:
        """``synonymous_sites`` and ``nonsynonymous_sites`` were always 0, whatever the sequences."""
        assert set(detect_selection_signatures("ATCG", ["ATCG"])) == {"diversity", "conservation"}

    def test_references_of_another_length_are_ignored_when_others_remain(self) -> None:
        result = detect_selection_signatures("ATCG", ["ATCGAAA", "TCG", "ATCA"])
        assert result["diversity"] == pytest.approx(0.25)

    @pytest.mark.parametrize("references", [[], ["ATCGAAA", "TCG"]])
    def test_no_comparable_reference_is_an_error_not_perfect_conservation(self, references: list[str]) -> None:
        """No reference used to give ``conservation = 1.0`` (empty list) or ``0.0`` (other lengths) out of thin air."""
        with pytest.raises(ValueError, match="same length"):
            detect_selection_signatures("ATCG", references)

    def test_ambiguous_bases_have_no_substitution_type(self) -> None:
        """An N used to count as a transversion, and a codon with an N as 'synonymous' (X == X)."""
        result = classify_mutations("ATGNAA", "ATGCAC")
        assert result["total"] == 2
        assert result["ambiguous"] == 1  # N -> C
        assert result["transversions"] == 1  # A -> C in the codon NAA ...
        assert result["transitions"] == 0
        assert result["synonymous"] == 0 and result["nonsynonymous"] == 0  # ... whose amino acid is unknown

    def test_types_account_for_every_difference(self) -> None:
        ancestral, derived = "ATGCCCAANTTA-GAC", "ATACCGAAGTTAAGAT"
        result = classify_mutations(ancestral, derived)
        assert result["total"] == _hamming(ancestral.upper(), derived.upper())
        assert result["transitions"] + result["transversions"] + result["ambiguous"] == result["total"]

    def test_unambiguous_input_is_classified_as_before(self) -> None:
        result = classify_mutations("TTTAAA", "TTCAAA")
        assert result == {
            "transitions": 1,
            "transversions": 0,
            "synonymous": 1,
            "nonsynonymous": 0,
            "ambiguous": 0,
            "total": 1,
        }

    def test_ti_tv_ratio_is_undefined_without_mutations(self) -> None:
        """0 transitions over 0 transversions used to be reported as an infinite ratio."""
        assert math.isnan(analyze_mutation_spectrum("ATCGATCG", "ATCGATCG")["ti_tv_ratio"])

    def test_ti_tv_ratio_is_infinite_with_transitions_only(self) -> None:
        assert analyze_mutation_spectrum("ATCA", "ATCG")["ti_tv_ratio"] == float("inf")  # A -> G

    def test_ti_tv_ratio_is_the_hand_counted_quotient(self) -> None:
        # reference ACGT -> GCGA: A->G (transition) and T->A (transversion); -> GTGA adds C->T (transition)
        assert analyze_mutation_spectrum("GCGA", "ACGT")["ti_tv_ratio"] == pytest.approx(1.0)
        assert analyze_mutation_spectrum("GTGA", "ACGT")["ti_tv_ratio"] == pytest.approx(2.0)
