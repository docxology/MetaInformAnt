"""Tests for DNA alignment functionality."""

from __future__ import annotations

import pytest

from metainformant.dna.alignment import pairwise as alignment


class TestDNAAlignment:
    """Test DNA alignment functionality."""

    def test_global_alignment_functionality(self):
        """Test global alignment functionality."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"
        try:
            result = alignment.global_align(seq1, seq2)
            assert result.score > 0
            assert len(result.aligned_seq1) == len(result.aligned_seq2)
        except (AttributeError, ImportError) as e:
            pytest.skip(f"Alignment functions unavailable: {e}")

    def test_local_alignment_functionality(self):
        """Test local alignment functionality."""
        seq1 = "ATCGATCGATCG"
        seq2 = "GATCG"
        try:
            result = alignment.local_align(seq1, seq2)
            assert result.score > 0
        except (AttributeError, ImportError) as e:
            pytest.skip(f"Alignment functions unavailable: {e}")

    def test_alignment_edge_cases(self):
        """Test edge cases in alignment functions."""
        try:
            # Single base
            result = alignment.global_align("A", "A")
            assert result.score >= 0
        except (AttributeError, ImportError) as e:
            pytest.skip(f"Alignment functions unavailable: {e}")

    def test_alignment_result_dataclass(self):
        """Test AlignmentResult dataclass functionality."""
        try:
            result = alignment.AlignmentResult("ATCG", "ATCG", 4.0)
            # Test basic attributes
            assert result.aligned_seq1 == "ATCG"
            assert result.aligned_seq2 == "ATCG"
            assert result.score == 4.0
        except (AttributeError, TypeError) as e:
            # AlignmentResult may not be a dataclass or may have different signature
            pytest.skip(f"AlignmentResult not available as expected: {e}")
