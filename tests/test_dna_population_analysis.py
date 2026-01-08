"""Tests for population genetics analysis orchestrators."""

from __future__ import annotations

import pytest

from metainformant.dna.population.analysis import (
    calculate_fst,
    calculate_nucleotide_diversity,
    calculate_tajima_d,
    calculate_wattersons_theta,
)


class TestCalculateSummaryStatistics:
    """Test comprehensive summary statistics calculator."""

    def test_sequence_based_statistics(self):
        """Test calculation from sequence data."""
        seqs = ["AAAA", "AAAT", "AATT"]
        stats = calculate_summary_statistics(sequences=seqs)
        
        assert "nucleotide_diversity" in stats
        assert "segregating_sites" in stats
        assert "wattersons_theta" in stats
        assert "tajimas_d" in stats
        assert stats["sample_size"] == 3
        assert stats["sequence_length"] == 4

    def test_genotype_matrix_statistics(self):
        """Test calculation from genotype matrix."""
        genotypes = [
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 1],
        ]
        stats = calculate_summary_statistics(genotype_matrix=genotypes)
        
        assert "allele_frequencies" in stats
        assert "observed_heterozygosity" in stats
        assert stats["sample_size"] == 3
        assert stats["num_sites"] == 3

    def test_both_sequence_and_genotype(self):
        """Test with both sequence and genotype data."""
        seqs = ["AAAA", "AAAT"]
        genotypes = [[0, 1], [0, 1]]
        stats = calculate_summary_statistics(sequences=seqs, genotype_matrix=genotypes)
        
        assert "nucleotide_diversity" in stats
        assert "allele_frequencies" in stats

    def test_empty_input(self):
        """Test with empty input."""
        stats = calculate_summary_statistics()
        assert stats == {}


class TestComparePopulations:
    """Test population comparison function."""

    def test_sequence_based_comparison(self):
        """Test comparison using sequences."""
        pop1 = ["AAAA", "AAAA", "AAAT"]
        pop2 = ["TTTT", "TTTT", "TTTA"]
        
        result = compare_populations(pop1_sequences=pop1, pop2_sequences=pop2)
        
        assert "pop1_stats" in result
        assert "pop2_stats" in result
        assert "fst" in result
        assert "differentiation" in result
        assert result["fst"] > 0.5  # High differentiation
        assert result["differentiation"] == "high"

    def test_identical_populations(self):
        """Test comparison of identical populations."""
        pop1 = ["AAAA", "AAAA"]
        pop2 = ["AAAA", "AAAA"]
        
        result = compare_populations(pop1_sequences=pop1, pop2_sequences=pop2)
        assert result["fst"] == 0.0
        assert result["differentiation"] == "none"

    def test_genotype_matrix_comparison(self):
        """Test comparison using genotype matrices."""
        pop1_genotypes = [[0, 0], [0, 0]]
        pop2_genotypes = [[1, 1], [1, 1]]
        
        result = compare_populations(
            pop1_genotypes=pop1_genotypes,
            pop2_genotypes=pop2_genotypes
        )
        
        assert "pop1_stats" in result
        assert "pop2_stats" in result
        assert result["fst"] is None  # Fst requires sequences

    def test_missing_data_error(self):
        """Test that missing data raises error."""
        with pytest.raises(ValueError, match="Must provide"):
            compare_populations(pop1_sequences=["AAAA"])


class TestNeutralityTestSuite:
    """Test neutrality test suite."""

    def test_basic_neutrality_tests(self):
        """Test basic neutrality test suite."""
        seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        results = neutrality_test_suite(seqs)
        
        assert "tajimas_d" in results
        assert "nucleotide_diversity" in results
        assert "wattersons_theta" in results
        assert "segregating_sites" in results
        assert "pi_theta_ratio" in results
        assert "interpretation" in results
        assert "sample_size" in results

    def test_interpretation_negative_d(self):
        """Test negative Tajima's D interpretation."""
        # Sequences that favor low-frequency variants (negative D)
        seqs = ["AAAA", "AAAT", "AATT", "ATTT", "TTTT"]
        results = neutrality_test_suite(seqs)
        
        assert results["interpretation"] in [
            "negative_d", "strong_negative_d", "neutral"
        ]

    def test_interpretation_positive_d(self):
        """Test positive Tajima's D interpretation."""
        # Sequences with intermediate frequency variants (positive D)
        seqs = ["AAAA", "AAAA", "AAAT", "AAAT", "AATT"]
        results = neutrality_test_suite(seqs)
        
        assert results["interpretation"] in [
            "positive_d", "strong_positive_d", "neutral"
        ]

    def test_pi_theta_ratio(self):
        """Test pi/theta ratio calculation."""
        seqs = ["AAAA", "AAAT", "AATT"]
        results = neutrality_test_suite(seqs)
        
        assert "pi_theta_ratio" in results
        # Ratio should be close to 1 under neutrality
        assert results["pi_theta_ratio"] >= 0.0

