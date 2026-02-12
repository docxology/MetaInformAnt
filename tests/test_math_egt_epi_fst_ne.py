from __future__ import annotations

from metainformant.math.epidemiology.models import basic_reproduction_number, sir_step
from metainformant.math.evolutionary_dynamics.core import replicator_derivative, replicator_step
from metainformant.math.population_genetics.effective_size import effective_size_sex_ratio, harmonic_mean_effective_size
from metainformant.math.population_genetics.fst import fst_from_allele_freqs, fst_from_heterozygosity


def test_replicator_step_and_derivative():
    x = [0.5, 0.5]
    A = [[1.0, 0.0], [0.0, 1.5]]
    x_next = replicator_step(x, A)
    # Strategy 2 should increase in frequency
    assert x_next[1] > x[1]
    dx = replicator_derivative(x, A)
    assert dx[1] > 0


def test_sir_and_r0():
    S, I, R = 0.99, 0.01, 0.0
    Sn, In, Rn = sir_step(S, I, R, beta=0.5, gamma=0.25, dt=0.1)
    assert Sn >= 0 and In >= 0 and Rn >= 0
    assert basic_reproduction_number(0.5, 0.25) == 2.0


def test_fst():
    fst = fst_from_allele_freqs([0.2, 0.8])
    # Hs = mean(2p(1-p)) = mean(0.32, 0.32) = 0.32; Ht = 2*(0.5)*(0.5)=0.5 => Fst=0.36
    assert abs(fst - 0.36) < 1e-12
    assert fst_from_heterozygosity(0.2, 0.5) == (0.5 - 0.2) / 0.5


def test_effective_size():
    Ne = harmonic_mean_effective_size([100, 50, 200])
    assert Ne > 0
    assert abs(effective_size_sex_ratio(40, 60) - (4 * 40 * 60) / (40 + 60)) < 1e-12
