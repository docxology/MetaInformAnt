"""Tests for population genetics simulation methods."""

from __future__ import annotations

import random

import numpy as np
import pytest

from metainformant.simulation.models.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    simulate_admixture,
    simulate_bottleneck_population,
    generate_two_populations,
    simulate_population_expansion,
)


class TestGeneratePopulationSequences:
    """Test population sequence generation."""

    def test_basic_generation(self):
        """Test basic sequence generation."""
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=100,
            mutation_rate=0.01,
        )
        assert len(seqs) == 10
        assert all(len(s) == 100 for s in seqs)

    def test_with_nucleotide_diversity(self):
        """Test generation with target nucleotide diversity."""
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=1000,
            nucleotide_diversity=0.01,
        )
        assert len(seqs) == 10

        # Check that sequences are not all identical
        unique_seqs = len(set(seqs))
        assert unique_seqs > 1

    def test_with_wattersons_theta(self):
        """Test generation with target Watterson's theta."""
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=1000,
            wattersons_theta=0.01,
        )
        assert len(seqs) == 10

    def test_with_reference_sequence(self):
        """Test generation with provided reference sequence."""
        ref_seq = "A" * 100
        seqs = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            reference_sequence=ref_seq,
            mutation_rate=0.01,
        )
        assert len(seqs) == 5
        assert all(len(s) == 100 for s in seqs)

    def test_reproducibility(self):
        """Test that results are reproducible with seed."""
        rng1 = random.Random(42)
        rng2 = random.Random(42)

        seqs1 = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            mutation_rate=0.01,
            rng=rng1,
        )
        seqs2 = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            mutation_rate=0.01,
            rng=rng2,
        )

        assert seqs1 == seqs2


class TestGenerateTwoPopulations:
    """Test two-population generation."""

    def test_basic_generation(self):
        """Test basic two-population generation."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.1,
        )
        assert len(pop1) == 10
        assert len(pop2) == 10
        assert all(len(s) == 1000 for s in pop1 + pop2)

    def test_population_differentiation(self):
        """Test that populations are differentiated."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.5,  # High differentiation
        )

        # Check that populations are different
        # (not all sequences identical between populations)
        pop1_consensus = set(pop1)
        pop2_consensus = set(pop2)

        # Should have some differences
        assert len(pop1_consensus) > 1 or len(pop2_consensus) > 1

    def test_high_fst(self):
        """Test high Fst scenario."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.8,  # Very high differentiation
        )
        assert len(pop1) == 10
        assert len(pop2) == 10


class TestGenerateGenotypeMatrix:
    """Test genotype matrix generation."""

    def test_basic_generation(self):
        """Test basic genotype matrix generation."""
        genotypes = generate_genotype_matrix(
            n_individuals=10,
            n_sites=5,
        )
        assert len(genotypes) == 10
        assert all(len(row) == 5 for row in genotypes)
        assert all(0 <= g <= 2 for row in genotypes for g in row)

    def test_with_allele_frequencies(self):
        """Test generation with specified allele frequencies."""
        freqs = [0.2, 0.3, 0.4, 0.1, 0.5]
        genotypes = generate_genotype_matrix(
            n_individuals=100,
            n_sites=5,
            allele_frequencies=freqs,
        )
        assert len(genotypes) == 100

        # Check that frequencies are approximately correct
        for site_idx in range(5):
            site_genotypes = [row[site_idx] for row in genotypes]
            # Average genotype / 2 ≈ allele frequency
            avg_genotype = sum(site_genotypes) / len(site_genotypes)
            estimated_freq = avg_genotype / 2.0
            assert abs(estimated_freq - freqs[site_idx]) < 0.2  # Allow some variance

    def test_hwe_vs_non_hwe(self):
        """Test Hardy-Weinberg equilibrium vs non-HWE."""
        genotypes_hwe = generate_genotype_matrix(
            n_individuals=100,
            n_sites=1,
            allele_frequencies=[0.5],
            hwe=True,
        )
        genotypes_non_hwe = generate_genotype_matrix(
            n_individuals=100,
            n_sites=1,
            allele_frequencies=[0.5],
            hwe=False,
        )

        # Both should have similar allele frequencies
        # But HWE should have more heterozygotes (2pq = 0.5 for p=0.5)
        hwe_het = sum(1 for row in genotypes_hwe if row[0] == 1)
        sum(1 for row in genotypes_non_hwe if row[0] == 1)

        # HWE should have ~50% heterozygotes for p=0.5
        assert hwe_het > 30  # Should be around 50

    def test_haploid(self):
        """Test haploid genotype generation."""
        genotypes = generate_genotype_matrix(
            n_individuals=10,
            n_sites=5,
            ploidy=1,
        )
        assert len(genotypes) == 10
        assert all(0 <= g <= 1 for row in genotypes for g in row)


class TestSimulateBottleneckPopulation:
    """Test bottleneck population simulation."""

    def test_basic_bottleneck(self):
        """Test basic bottleneck simulation."""
        seqs = simulate_bottleneck_population(
            n_sequences=20,
            sequence_length=1000,
            bottleneck_size=5,
            bottleneck_duration=10,
        )
        assert len(seqs) == 20
        assert all(len(s) == 1000 for s in seqs)

    def test_bottleneck_reduces_diversity(self):
        """Test that bottleneck reduces diversity."""
        seqs = simulate_bottleneck_population(
            n_sequences=20,
            sequence_length=1000,
            pre_bottleneck_diversity=0.01,
            bottleneck_size=2,  # Severe bottleneck
            bottleneck_duration=10,
        )
        # Should have reduced diversity (many sequences similar)
        unique_seqs = len(set(seqs))
        # With severe bottleneck, many sequences should be similar
        assert unique_seqs < len(seqs)


class TestSimulatePopulationExpansion:
    """Test population expansion simulation."""

    def test_basic_expansion(self):
        """Test basic expansion simulation."""
        seqs = simulate_population_expansion(
            n_sequences=20,
            sequence_length=1000,
            expansion_factor=10.0,
        )
        assert len(seqs) == 20
        assert all(len(s) == 1000 for s in seqs)

    def test_expansion_increases_sample_size(self):
        """Test that expansion creates more sequences."""
        seqs = simulate_population_expansion(
            n_sequences=50,
            sequence_length=1000,
            expansion_factor=10.0,
        )
        assert len(seqs) == 50


class TestGenerateSiteFrequencySpectrum:
    """Test site frequency spectrum generation."""

    def test_folded_sfs(self):
        """Test folded SFS generation."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=True,
        )
        assert len(sfs) == 5  # n//2 for n=10
        assert sum(sfs) == 100

    def test_unfolded_sfs(self):
        """Test unfolded SFS generation."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=False,
        )
        assert len(sfs) == 9  # n-1 for n=10
        assert sum(sfs) == 100

    def test_sfs_rare_alleles(self):
        """Test that SFS has more rare alleles."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=True,
        )
        # Under neutral model, rare alleles should be more common
        # (first bin should have more sites)
        assert sfs[0] > 0  # Should have some rare alleles


class TestGenerateLinkageDisequilibriumData:
    """Test linkage disequilibrium data generation."""

    def test_basic_ld_generation(self):
        """Test basic LD data generation."""
        genotypes = generate_linkage_disequilibrium_data(
            n_individuals=100,
            n_sites=10,
            r_squared_target=0.5,
        )
        assert len(genotypes) == 100
        assert all(len(row) == 10 for row in genotypes)
        assert all(0 <= g <= 2 for row in genotypes for g in row)

    def test_with_allele_frequencies(self):
        """Test LD generation with specified frequencies."""
        freqs = [0.3, 0.4, 0.2, 0.5, 0.1]
        genotypes = generate_linkage_disequilibrium_data(
            n_individuals=100,
            n_sites=5,
            allele_frequencies=freqs,
        )
        assert len(genotypes) == 100


class TestBottleneckStatsInterface:
    """Test the statistics interface of simulate_bottleneck_population.

    Regression tests: this interface previously called ``rng.poisson`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_stats_interface_returns_trajectory(self):
        """Test that the default interface returns simulation statistics."""
        result = simulate_bottleneck_population(
            initial_size=100,
            bottleneck_size=10,
            final_size=100,
            generations=6,
            rng=random.Random(42),
        )
        assert isinstance(result, dict)
        assert len(result["population_sizes"]) == 6
        assert result["population_sizes"][:2] == [100, 100]  # pre-bottleneck
        assert result["bottleneck_start"] == 2
        assert result["bottleneck_end"] == 4
        assert result["population_sizes"][3] == 10  # during bottleneck
        assert result["total_mutations"] >= 0
        assert 0.0 < result["final_diversity"] <= 1.0

    def test_stats_interface_is_deterministic(self):
        """Test reproducibility with a seeded RNG."""
        kwargs = dict(initial_size=50, bottleneck_size=5, final_size=50, generations=6)
        result_a = simulate_bottleneck_population(**kwargs, rng=random.Random(7))
        result_b = simulate_bottleneck_population(**kwargs, rng=random.Random(7))
        assert result_a == result_b


class TestExpansionStatsInterface:
    """Test the statistics interface of simulate_population_expansion.

    Regression tests: this interface previously called ``rng.poisson`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_stats_interface_returns_trajectory(self):
        """Test that the default interface returns simulation statistics."""
        result = simulate_population_expansion(
            initial_size=100,
            final_size=1000,
            expansion_time=5,
            rng=random.Random(42),
        )
        assert isinstance(result, dict)
        assert len(result["population_sizes"]) == 6  # expansion_time + 1
        assert result["population_sizes"][0] == 100
        assert result["population_sizes"][-1] == pytest.approx(1000, abs=1)
        assert result["growth_rate"] > 0
        assert result["total_mutations"] >= 0
        assert result["final_diversity"] <= 1.0


class TestSFSExpansionModel:
    """Test the expansion demographic model of generate_site_frequency_spectrum.

    Regression tests: this model previously called ``rng.zipf`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_expansion_model_default_alpha(self):
        """Test that the default alpha (1.0) produces a valid SFS."""
        sfs = generate_site_frequency_spectrum(
            n_samples=10,
            n_sites=100,
            demographic_model="expansion",
            rng=random.Random(42),
        )
        assert sum(sfs) == 100
        assert all(count >= 0 for count in sfs)

    def test_expansion_model_favors_rare_alleles(self):
        """Test that zipf-distributed frequencies favor low-frequency variants."""
        sfs = generate_site_frequency_spectrum(
            n_samples=50,
            n_sites=500,
            demographic_model="expansion",
            parameters={"alpha": 2.0},
            rng=random.Random(42),
        )
        assert sfs[0] > sfs[-1]  # singletons more common than high-frequency variants

    def test_expansion_model_custom_alpha_changes_spectrum(self):
        """Test that a provided alpha parameter is honored."""
        kwargs = dict(n_samples=10, n_sites=100, demographic_model="expansion")
        sfs_default = generate_site_frequency_spectrum(**kwargs, rng=random.Random(42))
        sfs_alpha2 = generate_site_frequency_spectrum(**kwargs, parameters={"alpha": 2.0}, rng=random.Random(42))
        assert sfs_default != sfs_alpha2


class TestSimulateAdmixture:
    """Test population admixture simulation.

    Regression tests: the drift step previously called ``rng.normal`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_basic_admixture(self):
        """Test that admixture runs and produces a frequency trajectory."""
        proportions = np.array([[0.9, 0.1], [0.1, 0.9]])
        result = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(42))
        assert len(result["frequency_trajectory"]) == 4  # generations + 1
        assert len(result["ancestral_frequencies"]) == 2
        assert result["generations"] == 3
        assert result["migration_rate"] == 0.01
        for freqs in result["frequency_trajectory"]:
            assert all(0.0 <= f <= 1.0 for f in freqs)

    def test_admixture_is_deterministic(self):
        """Test reproducibility with a seeded RNG."""
        proportions = np.array([[0.9, 0.1], [0.1, 0.9]])
        result_a = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(7))
        result_b = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(7))
        assert all(
            np.array_equal(a, b) for a, b in zip(result_a["frequency_trajectory"], result_b["frequency_trajectory"])
        )
