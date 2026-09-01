"""Tests for metainformant.core.utils.seeds (deterministic seed derivation)."""

from __future__ import annotations

import pytest

from metainformant.core.utils.seeds import deterministic_replicate_seeds


class TestDeterministicReplicateSeeds:
    def test_count_and_range(self) -> None:
        seeds = deterministic_replicate_seeds(1234, 8)
        assert len(seeds) == 8
        assert all(0 <= s < 2**63 for s in seeds)

    def test_deterministic_across_calls(self) -> None:
        assert deterministic_replicate_seeds(7, 4) == deterministic_replicate_seeds(7, 4)

    def test_order_independent_distinct_values(self) -> None:
        seeds = deterministic_replicate_seeds(99, 10)
        assert len(set(seeds)) == 10

    def test_matches_math_domain_reexport(self) -> None:
        """The math.population_genetics re-export stays behaviorally identical."""
        from metainformant.math.population_genetics.statistics import (
            deterministic_replicate_seeds as math_seeds,
        )

        assert math_seeds(5, 3) == deterministic_replicate_seeds(5, 3)

    def test_nonpositive_count_raises(self) -> None:
        with pytest.raises(ValueError, match="n_replicates"):
            deterministic_replicate_seeds(1, 0)
