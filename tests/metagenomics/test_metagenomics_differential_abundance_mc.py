"""Tests for the Dirichlet Monte Carlo implementation of aldex2_like DA."""

from __future__ import annotations

import pytest

from metainformant.metagenomics.comparative.differential_abundance import (
    differential_abundance,
)

COUNTS = [
    [500, 10, 30],
    [600, 12, 28],
    [550, 9, 31],
    [520, 11, 29],
    [100, 40, 25],
    [90, 45, 27],
    [110, 38, 26],
    [95, 42, 24],
]
GROUPS = [0, 0, 0, 0, 1, 1, 1, 1]
TAXA = ["taxA", "taxB", "taxC"]


def test_same_seed_reproduces_exactly() -> None:
    """Identical seeds must give bit-identical Monte Carlo results."""
    res1 = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=64, seed=11)
    res2 = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=64, seed=11)

    assert res1 == res2


def test_different_seed_gives_different_draws() -> None:
    """Monte Carlo draws are stochastic: another seed perturbs p-values."""
    res1 = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=2, seed=11)
    res2 = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=2, seed=12)

    p1 = [r["p_value"] for r in res1]
    p2 = [r["p_value"] for r in res2]
    assert p1 != p2


def test_strong_difference_is_detected() -> None:
    """The strongly shifted taxon is significant; the null taxon is not."""
    results = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=128, seed=7)
    by_taxon = {r["taxon"]: r for r in results}

    assert by_taxon["taxA"]["adjusted_p"] < 0.05
    # taxC counts are nearly identical between groups: no signal.
    assert by_taxon["taxC"]["p_value"] > 0.01


def test_result_contract_preserved() -> None:
    """Output keys stay backward compatible and values stay sane."""
    results = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=32, seed=3)

    assert len(results) == 3
    for row in results:
        assert set(row) >= {
            "taxon",
            "log2fc",
            "p_value",
            "adjusted_p",
            "effect_size",
            "mean_group1",
            "mean_group2",
        }
        assert 0.0 <= row["p_value"] <= 1.0
        assert 0.0 <= row["adjusted_p"] <= 1.0
        # taxA: 500-600 counts in group1 vs 90-110 in group2 -> positive log2(g1/g2).
        if row["taxon"] == "taxA":
            assert row["log2fc"] > 0


def test_n_monte_carlo_validation() -> None:
    """n_monte_carlo must be a positive number of draws."""
    with pytest.raises(ValueError, match="n_monte_carlo"):
        differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=0)

    # A single draw is still a valid (if noisy) Monte Carlo instance.
    results = differential_abundance(COUNTS, GROUPS, TAXA, n_monte_carlo=1, seed=5)
    assert len(results) == 3
    assert all(0.0 <= r["p_value"] <= 1.0 for r in results)
