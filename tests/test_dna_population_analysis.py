"""Tests for population genetics analysis functions."""

from __future__ import annotations

import pytest

from metainformant.dna.population.analysis import (
    calculate_fst,
    calculate_nucleotide_diversity,
    calculate_tajima_d,
    calculate_wattersons_theta,
)


class TestCalculateNucleotideDiversity:
    """Test nucleotide diversity calculation."""

    def test_basic_diversity(self):
        """Test basic nucleotide diversity calculation."""
        seqs = ["AAAA", "AAAT", "AATT"]
        pi = calculate_nucleotide_diversity(seqs)
        assert pi > 0

    def test_identical_sequences(self):
        """Test with identical sequences."""
        seqs = ["AAAA", "AAAA"]
        pi = calculate_nucleotide_diversity(seqs)
        assert pi == 0.0


class TestCalculateFst:
    """Test Fst calculation."""

    def test_differentiated_populations(self):
        """Test with differentiated populations."""
        pop1 = ["AAAA", "AAAA"]
        pop2 = ["TTTT", "TTTT"]
        fst = calculate_fst(pop1, pop2)
        assert 0 <= fst <= 1

    def test_identical_populations(self):
        """Test with identical populations."""
        pop1 = ["AAAA", "AAAT"]
        pop2 = ["AAAA", "AAAT"]
        fst = calculate_fst(pop1, pop2)
        assert fst >= 0.0


class TestCalculateTajimaD:
    """Test Tajima's D calculation."""

    def test_basic_calculation(self):
        """Test basic Tajima's D calculation."""
        seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        result = calculate_tajima_d(seqs)
        # Returns tuple (d, p_value) or float depending on implementation
        if isinstance(result, tuple):
            d, p_value = result
            assert isinstance(d, (int, float))
        else:
            assert isinstance(result, float)

    def test_no_variation(self):
        """Test with no variation."""
        seqs = ["AAAA", "AAAA", "AAAA", "AAAA"]
        result = calculate_tajima_d(seqs)
        # Returns tuple (d, p_value) or float
        if isinstance(result, tuple):
            d, p_value = result
            assert d == 0.0
        else:
            assert result == 0.0


class TestCalculateWattersonsTheta:
    """Test Watterson's theta calculation."""

    def test_basic_calculation(self):
        """Test basic Watterson's theta calculation."""
        seqs = ["AAAA", "AAAT", "AATT"]
        theta = calculate_wattersons_theta(seqs)
        assert theta >= 0

    def test_no_variation(self):
        """Test with no variation."""
        seqs = ["AAAA", "AAAA"]
        theta = calculate_wattersons_theta(seqs)
        assert theta == 0.0
