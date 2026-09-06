"""Tests for DNA alignment functionality."""

from __future__ import annotations

import pytest

from metainformant.dna.alignment import pairwise as alignment


class TestDNAAlignment:
    """Test DNA alignment functionality."""

    def test_global_alignment_functionality(self):
        """Global alignment of identical sequences reproduces the input."""
        result = alignment.global_align("ATCGATCG", "ATCGATCG")
        assert result.score == 8.0
        assert result.seq1_aligned == "ATCGATCG"
        assert result.seq2_aligned == "ATCGATCG"

    def test_global_alignment_introduces_gaps(self):
        """Gap penalties apply when lengths differ."""
        result = alignment.global_align("AAAA", "AA")
        assert result.score == 2 * 1 + 2 * -2
        assert len(result.seq1_aligned) == len(result.seq2_aligned)
        assert result.seq1_aligned == "AAAA"
        assert result.seq2_aligned == "--AA"

    def test_local_alignment_functionality(self):
        """Local alignment extracts the shared substring with no phantom bases."""
        result = alignment.local_align("TTTACGTTTA", "ACGT")
        assert result.score == 4.0
        assert result.seq1_aligned == "ACGT"
        assert result.seq2_aligned == "ACGT"

    def test_local_alignment_single_base(self):
        """A border-touching traceback must not append phantom characters."""
        result = alignment.local_align("A", "A")
        assert result.score == 1.0
        assert result.seq1_aligned == "A"
        assert result.seq2_aligned == "A"
        assert result.start_positions == (0, 0)

    def test_local_alignment_no_positive_score(self):
        """Local alignment with no matches returns an empty alignment."""
        result = alignment.local_align("AAAA", "TTTT", match=1, mismatch=-1, gap=-2)
        assert result.score == 0.0
        assert result.seq1_aligned == ""
        assert result.seq2_aligned == ""

    def test_alignment_edge_cases(self):
        """Empty sequences are rejected."""
        with pytest.raises(ValueError):
            alignment.global_align("", "ACGT")
        with pytest.raises(ValueError):
            alignment.local_align("ACGT", "")

    def test_alignment_result_dataclass(self):
        """AlignmentResult exposes aligned/seq1_aligned style accessors."""
        result = alignment.AlignmentResult("ATCG", "ATCG", 4.0)
        assert result.seq1_aligned == "ATCG"
        assert result.seq2_aligned == "ATCG"
        assert result.score == 4.0
        assert result.aligned_seq1 == "ATCG"
        assert result.aligned_seq2 == "ATCG"
