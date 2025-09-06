from __future__ import annotations

from metainformant.math import fixation_probability, hardy_weinberg_genotype_freqs, mutation_update, selection_update


def test_hardy_weinberg_genotype_freqs_basic():
    p2, two_pq, q2 = hardy_weinberg_genotype_freqs(0.25)
    assert abs(p2 - 0.0625) < 1e-9
    assert abs(two_pq - 0.375) < 1e-9
    assert abs(q2 - 0.5625) < 1e-9


def test_selection_update_balancing_and_directional():
    # Directional selection favoring A (w_AA > w_aa)
    p_next = selection_update(0.2, fitness_AA=1.2, fitness_Aa=1.1, fitness_aa=1.0)
    assert p_next > 0.2

    # No selection: mean fitness cancels out
    p_same = selection_update(0.2, 1.0, 1.0, 1.0)
    assert abs(p_same - 0.2) < 1e-12


def test_mutation_update_forward_and_back():
    p_next = mutation_update(0.5, mu=0.01, nu=0.02)
    # Expected: 0.5*(0.99) + 0.5*0.02 = 0.495 + 0.01 = 0.505
    assert abs(p_next - 0.505) < 1e-12


def test_fixation_probability_limits():
    # Neutral case equals initial frequency
    assert abs(fixation_probability(0.1, 1000, 0.0) - 0.1) < 1e-12
    # Certain boundaries
    assert fixation_probability(0.0, 1000, 0.1) == 0.0
    assert fixation_probability(1.0, 1000, -0.1) == 1.0
