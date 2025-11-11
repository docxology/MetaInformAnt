"""Enhanced tests for DNA alignment functionality."""

import pytest
from metainformant.dna import alignment


class TestDNAAlignmentEnhanced:
    """Test enhanced DNA alignment functionality."""

    def test_alignment_identity_calculation(self):
        """Test alignment identity calculation."""
        # Perfect match
        result = alignment.AlignmentResult("ATCG", "ATCG", 4.0)
        identity = alignment.calculate_alignment_identity(result)
        assert identity == 100.0

        # Partial match with gaps
        result = alignment.AlignmentResult("ATCG", "AT-G", 3.0)
        identity = alignment.calculate_alignment_identity(result)
        # In "ATCG" vs "AT-G": 3 aligned positions (A-A, T-T, G-G), all match = 100%
        # But if counting original positions: 3 matches out of 4 = 75%
        # Function counts aligned positions (excluding gaps), so 3/3 = 100%
        assert identity == 100.0

        # No match (actually has 1 match: G-G in reverse complement)
        result = alignment.AlignmentResult("ATCG", "TGCA", 0.0)
        identity = alignment.calculate_alignment_identity(result)
        # "ATCG" vs "TGCA": A-T, T-G, C-C, G-A = 1 match out of 4 = 25%
        assert identity == 25.0

    def test_conserved_regions_finding(self):
        """Test finding conserved regions in alignments."""
        # Perfect conservation
        result = alignment.AlignmentResult("ATCGATCG", "ATCGATCG", 8.0)
        conserved = alignment.find_conserved_regions(result, min_length=3)
        assert len(conserved) >= 1
        assert conserved[0][0] == "ATCGATCG"  # Full sequence

        # Mixed conservation with gaps
        result = alignment.AlignmentResult("ATCGATCG", "AT-GATCG", 6.0)
        conserved = alignment.find_conserved_regions(result, min_length=2)
        # For "ATCGATCG" vs "AT-GATCG": finds "AT" (0-2) and "GATCG" (3-8) = 2 regions
        assert len(conserved) >= 2

        # No conserved regions
        result = alignment.AlignmentResult("ATCG", "TGCA", 0.0)
        conserved = alignment.find_conserved_regions(result, min_length=3)
        assert len(conserved) == 0

    def test_alignment_statistics(self):
        """Test comprehensive alignment statistics."""
        result = alignment.AlignmentResult("ATCGATCG", "AT-GATCG", 6.0)
        stats = alignment.alignment_statistics(result)

        assert stats["length1"] == 8
        assert stats["length2"] == 8
        # For "ATCGATCG" vs "AT-GATCG": 7 aligned positions (A-A, T-T, G-A, A-T, T-C, C-G, G-G)
        # All 7 positions match (A-A, T-T, G-G, A-A, T-T, C-C, G-G) = 7 matches
        assert stats["matches"] == 7
        assert stats["mismatches"] == 0
        assert stats["gaps1"] == 0
        assert stats["gaps2"] == 1
        assert stats["identity"] == 100.0  # 7/7 aligned positions match
        assert stats["score"] == 6.0

    def test_global_alignment_functionality(self):
        """Test global alignment functionality."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"
        result = alignment.global_align(seq1, seq2)

        assert result.score > 0
        assert len(result.aligned_seq1) == len(result.aligned_seq2)
        assert result.aligned_seq1.replace("-", "") == seq1
        assert result.aligned_seq2.replace("-", "") == seq2

    def test_local_alignment_functionality(self):
        """Test local alignment functionality."""
        seq1 = "ATCGATCGATCG"
        seq2 = "GATCG"
        result = alignment.local_align(seq1, seq2)

        assert result.score > 0
        assert "GATCG" in result.aligned_seq1.replace("-", "")
        assert "GATCG" in result.aligned_seq2.replace("-", "")

    def test_alignment_edge_cases(self):
        """Test edge cases in alignment functions."""
        # Empty sequences
        result = alignment.global_align("", "")
        assert result.score == 0.0
        assert result.aligned_seq1 == ""
        assert result.aligned_seq2 == ""

        # Single base
        result = alignment.global_align("A", "A")
        assert result.score > 0

        # Very different sequences
        result = alignment.global_align("AAAA", "TTTT")
        identity = alignment.calculate_alignment_identity(result)
        assert identity < 50.0  # Should be low similarity

    def test_alignment_result_dataclass(self):
        """Test AlignmentResult dataclass functionality."""
        result = alignment.AlignmentResult("ATCG", "ATCG", 4.0)

        # Test that it's a proper dataclass
        assert result.aligned_seq1 == "ATCG"
        assert result.aligned_seq2 == "ATCG"
        assert result.score == 4.0

        # Test that it's hashable (for use in sets/dicts if needed)
        result_set = {result}
        assert len(result_set) == 1

        # Test that it's comparable
        result2 = alignment.AlignmentResult("ATCG", "ATCG", 4.0)
        assert result == result2

        result3 = alignment.AlignmentResult("ATCG", "ATCG", 5.0)
        assert result != result3
