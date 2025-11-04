"""Comprehensive tests for DNA population genetics module."""

from __future__ import annotations

import pytest

from metainformant.dna import population


class TestAlleleFrequencies:
    """Tests for allele_frequencies function."""

    def test_basic_allele_frequencies(self):
        """Test basic allele frequency calculation."""
        genotype_matrix = [
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        freqs = population.allele_frequencies(genotype_matrix)
        assert freqs == [0.5, 0.75, 0.5]

    def test_empty_matrix(self):
        """Test with empty genotype matrix."""
        freqs = population.allele_frequencies([])
        assert freqs == []

    def test_single_individual(self):
        """Test with single individual."""
        genotype_matrix = [[0, 1, 0]]
        freqs = population.allele_frequencies(genotype_matrix)
        assert freqs == [0.0, 1.0, 0.0]

    def test_all_zeros(self):
        """Test with all zeros (no alternate alleles)."""
        genotype_matrix = [[0, 0], [0, 0]]
        freqs = population.allele_frequencies(genotype_matrix)
        assert freqs == [0.0, 0.0]

    def test_all_ones(self):
        """Test with all ones (all alternate alleles)."""
        genotype_matrix = [[1, 1], [1, 1]]
        freqs = population.allele_frequencies(genotype_matrix)
        assert freqs == [1.0, 1.0]


class TestObservedHeterozygosity:
    """Tests for observed_heterozygosity function."""

    def test_basic_heterozygosity(self):
        """Test basic heterozygosity calculation."""
        genotypes = [(0, 0), (0, 1), (1, 1), (1, 0)]
        h_obs = population.observed_heterozygosity(genotypes)
        assert h_obs == 0.5

    def test_empty_genotypes(self):
        """Test with empty genotype list."""
        h_obs = population.observed_heterozygosity([])
        assert h_obs == 0.0

    def test_all_homozygous(self):
        """Test with all homozygous genotypes."""
        genotypes = [(0, 0), (1, 1), (0, 0)]
        h_obs = population.observed_heterozygosity(genotypes)
        assert h_obs == 0.0

    def test_all_heterozygous(self):
        """Test with all heterozygous genotypes."""
        genotypes = [(0, 1), (1, 0), (0, 1)]
        h_obs = population.observed_heterozygosity(genotypes)
        assert h_obs == 1.0


class TestNucleotideDiversity:
    """Tests for nucleotide_diversity function."""

    def test_basic_diversity(self):
        """Test basic nucleotide diversity calculation."""
        seqs = ["AAAA", "AAAT"]
        pi = population.nucleotide_diversity(seqs)
        assert abs(pi - 0.25) < 1e-9

    def test_three_sequences(self):
        """Test with three sequences."""
        seqs = ["AAAA", "AAAT", "AATT"]
        pi = population.nucleotide_diversity(seqs)
        assert pi > 0
        assert pi < 1

    def test_single_sequence(self):
        """Test with single sequence."""
        seqs = ["AAAA"]
        pi = population.nucleotide_diversity(seqs)
        assert pi == 0.0

    def test_empty_sequences(self):
        """Test with empty sequence list."""
        seqs = []
        pi = population.nucleotide_diversity(seqs)
        assert pi == 0.0

    def test_identical_sequences(self):
        """Test with all identical sequences."""
        seqs = ["AAAA", "AAAA", "AAAA"]
        pi = population.nucleotide_diversity(seqs)
        assert pi == 0.0

    def test_different_length_sequences(self):
        """Test with sequences of different lengths (truncates to shortest)."""
        seqs = ["AAAA", "AAATCG"]
        pi = population.nucleotide_diversity(seqs)
        # Should truncate to length 4
        assert abs(pi - 0.25) < 1e-9

    def test_zero_length_sequences(self):
        """Test with zero-length sequences."""
        seqs = ["", ""]
        pi = population.nucleotide_diversity(seqs)
        assert pi == 0.0


class TestTajimasD:
    """Tests for tajimas_d function."""

    def test_basic_tajimas_d(self):
        """Test basic Tajima's D calculation."""
        seqs = ["AAAA", "AAAT", "AATT"]
        d = population.tajimas_d(seqs)
        assert isinstance(d, float)

    def test_no_segregating_sites(self):
        """Test with no segregating sites."""
        seqs = ["AAAA", "AAAA", "AAAA"]
        d = population.tajimas_d(seqs)
        assert d == 0.0

    def test_single_sequence(self):
        """Test with single sequence."""
        seqs = ["AAAA"]
        d = population.tajimas_d(seqs)
        assert d == 0.0

    def test_insufficient_sequences(self):
        """Test with less than 2 sequences."""
        seqs = ["AAAA"]
        d = population.tajimas_d(seqs)
        assert d == 0.0


class TestHudsonFst:
    """Tests for hudson_fst function."""

    def test_fixed_differences(self):
        """Test with fixed differences between populations."""
        pop1 = ["AAAA", "AAAA"]
        pop2 = ["TTTT", "TTTT"]
        fst = population.hudson_fst(pop1, pop2)
        assert abs(fst - 1.0) < 1e-9

    def test_identical_populations(self):
        """Test with identical populations."""
        pop1 = ["AAAA", "AAAT"]
        pop2 = ["AAAA", "AAAT"]
        fst = population.hudson_fst(pop1, pop2)
        assert abs(fst - 0.0) < 1e-6  # Should be very close to 0

    def test_empty_population(self):
        """Test with empty population."""
        pop1 = ["AAAA"]
        pop2 = []
        fst = population.hudson_fst(pop1, pop2)
        assert fst == 0.0

    def test_different_length_sequences(self):
        """Test with sequences of different lengths."""
        pop1 = ["AAAA", "AAAT"]
        pop2 = ["TTTT", "TTTTCG"]
        fst = population.hudson_fst(pop1, pop2)
        # Should truncate to shortest length
        assert 0.0 <= fst <= 1.0


class TestSegregatingSites:
    """Tests for segregating_sites function."""

    def test_basic_segregating_sites(self):
        """Test basic segregating sites count."""
        seqs = ["AAAA", "AAAT", "AATT"]
        S = population.segregating_sites(seqs)
        assert S >= 1

    def test_no_segregating_sites(self):
        """Test with no segregating sites."""
        seqs = ["AAAA", "AAAA"]
        S = population.segregating_sites(seqs)
        assert S == 0

    def test_single_sequence(self):
        """Test with single sequence."""
        seqs = ["AAAA"]
        S = population.segregating_sites(seqs)
        assert S == 0

    def test_empty_sequences(self):
        """Test with empty sequence list."""
        seqs = []
        S = population.segregating_sites(seqs)
        assert S == 0


class TestWattersonsTheta:
    """Tests for wattersons_theta function."""

    def test_basic_wattersons_theta(self):
        """Test basic Watterson's theta calculation."""
        seqs = ["AAAA", "AAAT", "AATT"]
        theta_w = population.wattersons_theta(seqs)
        assert theta_w > 0

    def test_no_segregating_sites(self):
        """Test with no segregating sites."""
        seqs = ["AAAA", "AAAA"]
        theta_w = population.wattersons_theta(seqs)
        assert theta_w == 0.0

    def test_single_sequence(self):
        """Test with single sequence."""
        seqs = ["AAAA"]
        theta_w = population.wattersons_theta(seqs)
        assert theta_w == 0.0

    def test_empty_sequences(self):
        """Test with empty sequence list."""
        seqs = []
        theta_w = population.wattersons_theta(seqs)
        assert theta_w == 0.0


class TestIntegration:
    """Integration tests combining multiple functions."""

    def test_diversity_and_segregating_sites(self):
        """Test that diversity and segregating sites are related."""
        seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        pi = population.nucleotide_diversity(seqs)
        S = population.segregating_sites(seqs)
        theta_w = population.wattersons_theta(seqs)
        
        assert pi >= 0
        assert S >= 0
        assert theta_w >= 0
        
        # If there are segregating sites, diversity should be positive
        if S > 0:
            assert pi > 0

    def test_fst_and_diversity_relationship(self):
        """Test that Fst relates to diversity within and between populations."""
        pop1 = ["AAAA", "AAAT"]
        pop2 = ["TTTT", "TTTA"]
        
        fst = population.hudson_fst(pop1, pop2)
        pi1 = population.nucleotide_diversity(pop1)
        pi2 = population.nucleotide_diversity(pop2)
        
        # Fst should be high when populations are very different
        assert 0.0 <= fst <= 1.0
        # Diversity within populations should be non-negative
        assert pi1 >= 0
        assert pi2 >= 0

