"""Regression tests for ``generate_point_mutations``: a local random generator and a working ``mutation_rate``.

REAL IMPLEMENTATION POLICY: the real function runs against independent statistical oracles (scipy's exact
binomial test and chi-square test); nothing is mocked or patched.

The function drew from Python's global ``random`` state with no way to seed it locally, and its
``mutation_rate`` argument had no effect on the output.
"""

from __future__ import annotations

import random

import pytest
from scipy import stats

from metainformant.dna.variation.mutations import generate_point_mutations

SEQUENCE = "A" * 20000


def _changed_positions(original: str, mutated: str) -> list[int]:
    return [i for i, (a, b) in enumerate(zip(original, mutated)) if a != b]


class TestMutationRate:
    @pytest.mark.parametrize("rate", [0.01, 0.05, 0.3])
    def test_mutated_fraction_follows_the_rate(self, rate: float) -> None:
        """Every site mutates independently with probability ``mutation_rate`` (exact binomial test)."""
        mutated = generate_point_mutations(SEQUENCE, mutation_rate=rate, rng=random.Random(1))
        changed = _changed_positions(SEQUENCE, mutated)
        assert len(mutated) == len(SEQUENCE)
        assert stats.binomtest(len(changed), len(SEQUENCE), rate).pvalue > 1e-3

    def test_replacement_bases_are_uniform_over_the_other_three(self) -> None:
        mutated = generate_point_mutations(SEQUENCE, mutation_rate=0.3, rng=random.Random(2))
        replacements = [mutated[i] for i in _changed_positions(SEQUENCE, mutated)]
        assert set(replacements) == {"C", "G", "T"}  # never the original base
        counts = [replacements.count(base) for base in "CGT"]
        assert stats.chisquare(counts).pvalue > 1e-3

    def test_mutations_are_spread_along_the_whole_sequence(self) -> None:
        """Mutations fall uniformly along the sequence, not clustered at one end (equal halves, chi-square test)."""
        mutated = generate_point_mutations(SEQUENCE, mutation_rate=0.05, rng=random.Random(3))
        changed = _changed_positions(SEQUENCE, mutated)
        halves = [sum(1 for i in changed if i < 10000), sum(1 for i in changed if i >= 10000)]
        assert stats.chisquare(halves).pvalue > 1e-3

    def test_rate_zero_changes_nothing_and_rate_one_changes_everything(self) -> None:
        sequence = "ACGTACGTAC"
        assert generate_point_mutations(sequence, mutation_rate=0.0, rng=random.Random(4)) == sequence
        everything = generate_point_mutations(sequence, mutation_rate=1.0, rng=random.Random(4))
        assert all(a != b for a, b in zip(sequence, everything))

    @pytest.mark.parametrize("rate", [-0.1, 1.5, float("nan")])
    def test_rate_must_be_a_probability(self, rate: float) -> None:
        with pytest.raises(ValueError, match="mutation_rate"):
            generate_point_mutations("ACGT", mutation_rate=rate)

    def test_rate_is_not_checked_when_an_exact_count_is_requested(self) -> None:
        """With ``num_mutations`` the rate does not apply (and was never validated), so existing callers still work."""
        result = generate_point_mutations("ACGTACGT", num_mutations=2, mutation_rate=7.0, rng=random.Random(5))
        assert len(_changed_positions("ACGTACGT", result)) == 2


class TestExactCount:
    def test_exact_number_of_distinct_positions_changes(self) -> None:
        mutated = generate_point_mutations("A" * 50, num_mutations=10, rng=random.Random(6))
        changed = _changed_positions("A" * 50, mutated)
        assert len(changed) == 10
        assert all(mutated[i] in "CGT" for i in changed)

    def test_exact_count_ignores_the_rate(self) -> None:
        for rate in (0.0, 0.001, 1.0):
            mutated = generate_point_mutations("ACGT" * 10, num_mutations=3, mutation_rate=rate, rng=random.Random(7))
            assert len(_changed_positions("ACGT" * 10, mutated)) == 3

    def test_zero_mutations_returns_the_uppercased_sequence(self) -> None:
        assert generate_point_mutations("acgtacgt", num_mutations=0, rng=random.Random(8)) == "ACGTACGT"

    @pytest.mark.parametrize("count", [-1, 5])
    def test_impossible_counts_raise(self, count: int) -> None:
        with pytest.raises(ValueError):
            generate_point_mutations("ACGT", num_mutations=count)

    def test_empty_sequence_is_returned_unchanged(self) -> None:
        assert generate_point_mutations("", num_mutations=3) == ""
        assert generate_point_mutations("", mutation_rate=2.0) == ""


class TestRandomSource:
    def test_same_rng_seed_gives_the_same_result_whatever_the_global_state(self) -> None:
        random.seed(111)
        first = generate_point_mutations("ACGT" * 25, num_mutations=12, rng=random.Random(42))
        random.seed(222)  # a different global state must not matter
        second = generate_point_mutations("ACGT" * 25, num_mutations=12, rng=random.Random(42))
        assert first == second

    def test_rate_mode_is_reproducible_with_a_local_generator(self) -> None:
        first = generate_point_mutations(SEQUENCE[:2000], mutation_rate=0.1, rng=random.Random(9))
        random.seed(3)
        second = generate_point_mutations(SEQUENCE[:2000], mutation_rate=0.1, rng=random.Random(9))
        assert first == second

    def test_different_seeds_give_different_results(self) -> None:
        first = generate_point_mutations("ACGT" * 25, num_mutations=12, rng=random.Random(1))
        second = generate_point_mutations("ACGT" * 25, num_mutations=12, rng=random.Random(2))
        assert first != second

    def test_local_generator_leaves_the_global_state_untouched(self) -> None:
        random.seed(5)
        before = random.getstate()
        generate_point_mutations("ACGT" * 25, num_mutations=12, rng=random.Random(1))
        generate_point_mutations("ACGT" * 25, mutation_rate=0.2, rng=random.Random(1))
        assert random.getstate() == before

    def test_without_rng_the_global_state_still_controls_the_result(self) -> None:
        """Back-compatible default: ``random.seed`` makes a call without ``rng`` reproducible."""
        random.seed(7)
        first = generate_point_mutations("ACGTACGTACGT", num_mutations=4)
        random.seed(7)
        second = generate_point_mutations("ACGTACGTACGT", num_mutations=4)
        assert first == second
