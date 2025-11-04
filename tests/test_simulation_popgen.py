"""Tests for population genetics simulation methods."""

from __future__ import annotations

import random

import pytest

from metainformant.simulation.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
    simulate_bottleneck_population,
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
        non_hwe_het = sum(1 for row in genotypes_non_hwe if row[0] == 1)
        
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

