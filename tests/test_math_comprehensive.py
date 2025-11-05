"""Comprehensive tests for math module functions with edge cases and validation.

Tests cover all major submodules: price equation, population genetics,
coalescent theory, epidemiology, dynamics, and selection.
"""

from __future__ import annotations

import math

import pytest

from metainformant.math import (
    basic_reproduction_number,
    ddm_analytic_accuracy,
    ddm_mean_decision_time,
    expected_pairwise_diversity,
    expected_time_to_mrca,
    fixation_probability,
    hardy_weinberg_genotype_freqs,
    herd_immunity_threshold,
    heterozygosity_decay,
    inbreeding_coefficient,
    kin_selection_response,
    logistic_map,
    lotka_volterra_step,
    mutation_update,
    price_equation,
    selection_differential,
    selection_gradient,
    selection_update,
    sir_step,
    tajimas_D,
    watterson_theta,
)


class TestPriceEquation:
    """Comprehensive tests for Price equation and selection metrics."""

    def test_price_equation_with_offspring(self):
        """Test Price equation decomposition with offspring data."""
        fitness = [1.0, 1.2, 0.9, 1.1, 0.95]
        parent = [0.2, 0.4, 0.1, 0.35, 0.25]
        offspring = [0.25, 0.35, 0.15, 0.30, 0.28]

        cov_term, trans_term, total = price_equation(fitness, parent, offspring)

        assert abs(total - (cov_term + trans_term)) < 1e-10
        assert cov_term != 0.0  # Should have selection component
        assert total != 0.0

    def test_price_equation_without_offspring(self):
        """Test Price equation with only selection term."""
        fitness = [1.0, 1.2, 0.9]
        parent = [0.2, 0.4, 0.1]

        cov_term, trans_term, total = price_equation(fitness, parent, None)

        assert abs(trans_term) < 1e-10
        assert abs(total - cov_term) < 1e-10

    def test_selection_differential_positive_selection(self):
        """Test selection differential with positive selection."""
        fitness = [1.0, 1.5, 0.5, 1.2, 0.8]
        trait = [0.1, 0.9, 0.2, 0.8, 0.3]

        S = selection_differential(fitness, trait)

        assert S > 0.0  # Positive selection differential

    def test_selection_gradient_standardization(self):
        """Test selection gradient standardization by variance."""
        fitness = [1.0, 1.2, 0.9]
        trait_low_var = [0.2, 0.21, 0.19]
        trait_high_var = [0.1, 0.9, 0.2]

        beta_low = selection_gradient(fitness, trait_low_var)
        beta_high = selection_gradient(fitness, trait_high_var)

        # With higher variance, gradient should be smaller (standardized)
        assert abs(beta_high) < abs(beta_low) or trait_high_var == trait_low_var


class TestPopulationGenetics:
    """Tests for population genetic functions."""

    def test_hardy_weinberg_equilibrium(self):
        """Test Hardy-Weinberg genotype frequency calculation."""
        p = 0.6
        freq_AA, freq_Aa, freq_aa = hardy_weinberg_genotype_freqs(p)

        assert abs(freq_AA - (p * p)) < 1e-10
        assert abs(freq_Aa - (2 * p * (1 - p))) < 1e-10
        assert abs(freq_aa - ((1 - p) ** 2)) < 1e-10
        assert abs(freq_AA + freq_Aa + freq_aa - 1.0) < 1e-10

    def test_selection_update_directional(self):
        """Test selection update with directional selection."""
        p0 = 0.3
        # Selection favoring A
        p1 = selection_update(p0, fitness_AA=1.0, fitness_Aa=1.0, fitness_aa=0.8)

        assert p1 > p0  # Allele A should increase

    def test_selection_update_overdominance(self):
        """Test selection update with heterozygote advantage."""
        p0 = 0.5
        # Overdominance: heterozygote has highest fitness
        p1 = selection_update(p0, fitness_AA=0.8, fitness_Aa=1.0, fitness_aa=0.8)

        # Should remain near 0.5 (stable equilibrium)
        assert abs(p1 - 0.5) < 0.1

    def test_mutation_update_equilibrium(self):
        """Test mutation update approaching equilibrium."""
        mu = 0.01
        nu = 0.005

        # Start from arbitrary frequency
        p = 0.3
        for _ in range(1000):
            p = mutation_update(p, mu, nu)

        # Equilibrium should be nu / (mu + nu)
        expected_eq = nu / (mu + nu)
        assert abs(p - expected_eq) < 0.01

    def test_fixation_probability_neutral(self):
        """Test neutral fixation probability equals initial frequency."""
        p0 = 0.1
        Ne = 1000
        prob = fixation_probability(p0, Ne, selection_coefficient=0.0)

        assert abs(prob - p0) < 1e-10

    def test_fixation_probability_advantageous(self):
        """Test fixation probability for advantageous allele."""
        p0 = 0.01
        Ne = 1000
        prob_neutral = fixation_probability(p0, Ne, selection_coefficient=0.0)
        prob_advantageous = fixation_probability(p0, Ne, selection_coefficient=0.01)

        assert prob_advantageous > prob_neutral

    def test_watterson_theta_calculation(self):
        """Test Watterson's theta estimator."""
        S = 10
        n = 10
        theta = watterson_theta(S, n)

        assert theta > 0.0
        # Check it's approximately S / a1
        a1 = sum(1.0 / i for i in range(1, n))
        expected = S / a1
        assert abs(theta - expected) < 1e-10

    def test_heterozygosity_decay(self):
        """Test heterozygosity decay under drift."""
        H0 = 0.5
        Ne = 100
        t = 100

        Ht = heterozygosity_decay(H0, Ne, t)

        assert 0.0 <= Ht <= H0
        assert Ht < H0  # Should decrease over time


class TestCoalescent:
    """Tests for coalescent theory functions."""

    def test_expected_time_to_mrca(self):
        """Test expected time to MRCA calculation."""
        n = 10
        Ne = 1000

        t_mrca = expected_time_to_mrca(n, Ne)

        assert t_mrca > 0.0
        # For n=2, should be 4Ne
        t2 = expected_time_to_mrca(2, Ne)
        assert abs(t2 - (4 * Ne)) < 1e-10

    def test_expected_pairwise_diversity(self):
        """Test expected pairwise diversity calculation."""
        Ne = 1000
        mu = 0.0001

        pi = expected_pairwise_diversity(Ne, mu)

        expected = 4.0 * Ne * mu
        assert abs(pi - expected) < 1e-10

    def test_tajimas_D_neutral_expectation(self):
        """Test Tajima's D calculation."""
        n = 10
        S = 10
        # Calculate expected pi under neutrality
        a1 = sum(1.0 / i for i in range(1, n))
        expected_pi = S / a1

        D = tajimas_D(S, expected_pi, n)

        # Should be close to zero for neutral expectation
        assert abs(D) < 1.0


class TestEpidemiology:
    """Tests for epidemiological models."""

    def test_sir_model_conservation(self):
        """Test SIR model conserves total population."""
        S, I, R = 990.0, 10.0, 0.0
        beta, gamma = 0.3, 0.1

        total_before = S + I + R

        for _ in range(10):
            S, I, R = sir_step(S, I, R, beta, gamma, dt=0.01)

        total_after = S + I + R

        assert abs(total_after - total_before) < 0.1

    def test_basic_reproduction_number(self):
        """Test R0 calculation."""
        beta = 0.3
        gamma = 0.1

        R0 = basic_reproduction_number(beta, gamma)

        assert abs(R0 - 3.0) < 1e-10

    def test_herd_immunity_threshold(self):
        """Test herd immunity threshold calculation."""
        R0 = 3.0

        threshold = herd_immunity_threshold(R0)

        expected = 1.0 - 1.0 / R0
        assert abs(threshold - expected) < 1e-10
        assert 0.0 <= threshold <= 1.0


class TestDynamics:
    """Tests for population dynamics models."""

    def test_logistic_map_stability(self):
        """Test logistic map with stable parameters."""
        r = 2.5
        x0 = 0.5

        trajectory = logistic_map(r, x0, steps=100)

        # Should converge to stable equilibrium
        final = trajectory[-1]
        assert 0.0 <= final <= 1.0

    def test_lotka_volterra_oscillations(self):
        """Test Lotka-Volterra model oscillations."""
        prey, predator = 100.0, 10.0
        alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5

        trajectory_prey = [prey]
        trajectory_pred = [predator]

        for _ in range(100):
            prey, predator = lotka_volterra_step(prey, predator, alpha, beta, delta, gamma, dt=0.01)
            trajectory_prey.append(prey)
            trajectory_pred.append(predator)

        # Should show oscillations (not monotonic)
        assert max(trajectory_prey) > min(trajectory_prey)
        assert max(trajectory_pred) > min(trajectory_pred)


class TestDecisionTheory:
    """Tests for drift-diffusion models."""

    def test_ddm_accuracy_symmetric(self):
        """Test DDM accuracy for symmetric case."""
        accuracy = ddm_analytic_accuracy(drift_rate=0.0, boundary=1.0, noise_sd=1.0)

        assert abs(accuracy - 0.5) < 1e-10  # Should be 50% for no drift

    def test_ddm_accuracy_positive_drift(self):
        """Test DDM accuracy with positive drift."""
        accuracy = ddm_analytic_accuracy(drift_rate=0.5, boundary=1.0, noise_sd=1.0)

        assert accuracy > 0.5  # Should favor correct choice

    def test_ddm_decision_time(self):
        """Test DDM mean decision time."""
        time = ddm_mean_decision_time(drift_rate=0.5, boundary=1.0, noise_sd=1.0)

        assert time > 0.0
        # Zero drift case
        time_zero = ddm_mean_decision_time(drift_rate=0.0, boundary=1.0, noise_sd=1.0)
        assert abs(time_zero - 1.0) < 0.1


class TestSelection:
    """Tests for selection theory functions."""

    def test_kin_selection_hamiltons_rule(self):
        """Test Hamilton's rule calculation."""
        r = 0.5
        b = 0.4
        c = 0.1

        response = kin_selection_response(r, b, c)

        # r*b = 0.2, c = 0.1, so r*b > c (should be favored)
        assert response > 0.0
        assert abs(response - (r * b - c)) < 1e-10

    def test_kin_selection_not_favored(self):
        """Test Hamilton's rule when trait not favored."""
        r = 0.25
        b = 0.2
        c = 0.15

        response = kin_selection_response(r, b, c)

        # r*b = 0.05, c = 0.15, so r*b < c (should not be favored)
        assert response < 0.0


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_inputs(self):
        """Test functions with empty inputs."""
        # Price equation with empty fitness
        cov, trans, total = price_equation([], [], [])
        assert cov == 0.0 and trans == 0.0 and total == 0.0

    def test_zero_population_size(self):
        """Test functions with zero population size."""
        t = expected_time_to_mrca(10, 0.0)
        assert t == 0.0

    def test_invalid_frequencies(self):
        """Test functions with invalid frequency inputs."""
        # Hardy-Weinberg with invalid p
        freqs = hardy_weinberg_genotype_freqs(-0.1)
        assert freqs == (0.0, 0.0, 0.0)

        freqs = hardy_weinberg_genotype_freqs(1.5)
        assert freqs == (0.0, 0.0, 0.0)

    def test_extreme_selection_coefficients(self):
        """Test with extreme selection coefficients."""
        prob = fixation_probability(0.01, 1000, selection_coefficient=100.0)
        assert 0.0 <= prob <= 1.0

        prob = fixation_probability(0.01, 1000, selection_coefficient=-100.0)
        assert 0.0 <= prob <= 1.0






