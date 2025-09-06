from __future__ import annotations

from metainformant.math import (
    equilibrium_heterozygosity_infinite_alleles,
    heterozygosity_decay,
    inbreeding_coefficient,
    island_model_update,
    mutation_selection_balance_dominant,
    mutation_selection_balance_recessive,
)


def test_heterozygosity_and_inbreeding_drift():
    H0 = 0.5
    Ne = 100.0
    t = 10
    expected_Ht = H0 * ((1.0 - 1.0 / (2.0 * Ne)) ** t)
    assert abs(heterozygosity_decay(H0, Ne, t) - expected_Ht) < 1e-12

    Ft = inbreeding_coefficient(Ne, t)
    expected_Ft = 1.0 - ((1.0 - 1.0 / (2.0 * Ne)) ** t)
    assert abs(Ft - expected_Ft) < 1e-12


def test_equilibrium_heterozygosity_infinite_alleles():
    Ne = 1000.0
    mu = 1e-5
    Heq = equilibrium_heterozygosity_infinite_alleles(Ne, mu)
    expected = (4.0 * Ne * mu) / (1.0 + 4.0 * Ne * mu)
    assert abs(Heq - expected) < 1e-15


def test_island_model_update_and_mutation_selection_balance():
    p = 0.2
    m = 0.1
    pm = 0.8
    assert abs(island_model_update(p, m, pm) - 0.26) < 1e-12

    # Mutation-selection balance approximations
    mu = 1e-6
    s = 1e-2
    q_rec = mutation_selection_balance_recessive(mu, s)
    q_dom = mutation_selection_balance_dominant(mu, s)
    assert abs(q_rec - (mu / s) ** 0.5) < 1e-18
    assert abs(q_dom - (mu / s)) < 1e-18
