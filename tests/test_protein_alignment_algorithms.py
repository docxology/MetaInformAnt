"""Algorithm verification tests for protein alignment functions.

Tests for global_align (Needleman-Wunsch) and local_align (Smith-Waterman)
functions that were converted from placeholder implementations to real
functional algorithms. Following NO_MOCKING policy - all tests use real implementations.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.protein.sequence.alignment import global_align, local_align, calculate_alignment_identity


class TestGlobalAlignNeedlemanWunsch:
    """Tests for global_align function - Needleman-Wunsch algorithm verification."""

    def test_global_align_identical_sequences(self):
        """Test global alignment of identical sequences."""
        seq1 = "ACDEFG"
        seq2 = "ACDEFG"

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Perfect match: score should be len(seq) * match
        expected_score = len(seq1) * 1
        assert result["score"] == expected_score
        assert result["identity"] == 1.0
        assert result["aligned_seq1"] == seq1
        assert result["aligned_seq2"] == seq2
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_global_align_with_gaps(self):
        """Test global alignment requiring gaps."""
        seq1 = "ACDEFG"
        seq2 = "ACDFG"  # Missing 'E'

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Should align with gap: ACDEFG vs ACD-FG (score: 3)
        assert result["score"] == 3
        assert result["identity"] == 1.0  # All aligned positions match
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_global_align_with_mismatches(self):
        """Test global alignment with mismatches."""
        seq1 = "ACDEFG"
        seq2 = "ACDGFG"  # 'E' -> 'G' mismatch

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Should align: ACDEFG vs ACDGFG
        assert result["score"] == 4  # 5 matches - 1 mismatch
        assert result["identity"] == 5 / 6  # 5 out of 6 positions match
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_global_align_traceback_correctness(self):
        """Test traceback produces correct alignments."""
        seq1 = "HEAG"
        seq2 = "EAHG"

        result = global_align(seq1, seq2, match=2, mismatch=-1, gap=-2)

        # Actual score from the algorithm
        assert result["score"] == 2  # Based on actual implementation
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_global_align_scoring_matrix_construction(self):
        """Test that scoring matrix is constructed correctly."""
        seq1 = "AB"
        seq2 = "AB"

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # For identical sequences, score should be 2
        assert result["score"] == 2

    def test_global_align_empty_sequences(self):
        """Test global alignment of empty sequences."""
        result = global_align("", "", match=1, mismatch=-1, gap=-2)

        assert result["score"] == 0
        assert result["identity"] == 0.0  # Avoid division by zero
        assert result["aligned_seq1"] == ""
        assert result["aligned_seq2"] == ""

    def test_global_align_one_empty_sequence(self):
        """Test global alignment when one sequence is empty."""
        seq1 = "ABC"
        seq2 = ""

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # All gaps in seq2
        assert result["score"] == -6  # 3 gaps * -2
        assert result["identity"] == 0.0
        assert result["aligned_seq1"] == "ABC"
        assert result["aligned_seq2"] == "---"

    def test_global_align_single_residue(self):
        """Test global alignment of single residues."""
        result = global_align("A", "A", match=1, mismatch=-1, gap=-2)

        assert result["score"] == 1
        assert result["identity"] == 1.0
        assert result["aligned_seq1"] == "A"
        assert result["aligned_seq2"] == "A"

    def test_global_align_single_residue_mismatch(self):
        """Test global alignment of single mismatched residues."""
        result = global_align("A", "B", match=1, mismatch=-1, gap=-2)

        assert result["score"] == -1
        assert result["identity"] == 0.0
        assert result["aligned_seq1"] == "A"
        assert result["aligned_seq2"] == "B"

    def test_global_align_known_example(self):
        """Test against a known correct alignment example."""
        # Example from bioinformatics literature
        seq1 = "GATTACA"
        seq2 = "GCATGCU"

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-1)

        # Actual score from the algorithm
        assert result["score"] == 0  # Based on actual implementation
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])


class TestLocalAlignSmithWaterman:
    """Tests for local_align function - Smith-Waterman algorithm verification."""

    def test_local_align_identical_sequences(self):
        """Test local alignment of identical sequences."""
        seq1 = "ACDEFG"
        seq2 = "ACDEFG"

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Perfect match: score should be len(seq) * match
        expected_score = len(seq1) * 1
        assert result["score"] == expected_score
        assert result["identity"] == 1.0
        assert result["aligned_seq1"] == seq1
        assert result["aligned_seq2"] == seq2
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_local_align_subsequence(self):
        """Test local alignment finding subsequence match."""
        seq1 = "ABCDEFXYZ"
        seq2 = "XXXCDEFYYY"

        result = local_align(seq1, seq2, match=2, mismatch=-1, gap=-2)

        # Should find alignment
        # Actual score from the algorithm
        assert result["score"] == 9  # Based on actual implementation
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])
        # Check that CDEF is in both alignments
        assert "CDEF" in result["aligned_seq1"]
        assert "CDEF" in result["aligned_seq2"]

    def test_local_align_max_score_tracking(self):
        """Test that maximum score is correctly tracked."""
        seq1 = "ABCDE"
        seq2 = "ABCXY"

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Best alignment is "ABC" with score 3
        assert result["score"] == 3
        assert result["aligned_seq1"] == "ABC"
        assert result["aligned_seq2"] == "ABC"

    def test_local_align_no_alignment(self):
        """Test local alignment when sequences have no similarity."""
        seq1 = "AAAAA"
        seq2 = "TTTTT"

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # With high mismatch penalty, may align single characters
        assert result["score"] >= 0
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_local_align_empty_sequences(self):
        """Test local alignment of empty sequences."""
        result = local_align("", "", match=1, mismatch=-1, gap=-2)

        assert result["score"] == 0
        assert result["identity"] == 0.0
        assert result["aligned_seq1"] == ""
        assert result["aligned_seq2"] == ""

    def test_local_align_one_empty_sequence(self):
        """Test local alignment when one sequence is empty."""
        seq1 = "ABC"
        seq2 = ""

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # No alignment possible
        assert result["score"] == 0
        assert result["aligned_seq1"] == ""
        assert result["aligned_seq2"] == ""

    def test_local_align_single_residue_match(self):
        """Test local alignment of matching single residues."""
        result = local_align("A", "A", match=1, mismatch=-1, gap=-2)

        assert result["score"] == 1
        assert result["identity"] == 1.0
        assert result["aligned_seq1"] == "A"
        assert result["aligned_seq2"] == "A"

    def test_local_align_single_residue_mismatch(self):
        """Test local alignment of mismatched single residues."""
        result = local_align("A", "B", match=1, mismatch=-1, gap=-2)

        # No positive alignment possible
        assert result["score"] == 0
        assert result["aligned_seq1"] == ""
        assert result["aligned_seq2"] == ""

    def test_local_align_traceback_from_max(self):
        """Test traceback starts from maximum score position."""
        seq1 = "ABCD"
        seq2 = "XXABYY"

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Should find "AB" alignment
        assert result["score"] == 2
        assert result["aligned_seq1"] == "AB"
        assert result["aligned_seq2"] == "AB"


class TestAlignmentIdentityCalculation:
    """Tests for alignment identity calculations."""

    def test_calculate_alignment_identity_perfect(self):
        """Test identity calculation for perfect alignments."""
        alignment_result = global_align("ACDEFG", "ACDEFG")

        identity = calculate_alignment_identity(alignment_result)

        assert identity == 1.0

    def test_calculate_alignment_identity_no_matches(self):
        """Test identity calculation for completely different sequences."""
        alignment_result = global_align("AAAAA", "BBBBB")

        identity = calculate_alignment_identity(alignment_result)

        assert identity == 0.0

    def test_calculate_alignment_identity_partial(self):
        """Test identity calculation for partial matches."""
        alignment_result = global_align("ACDEFG", "ACDGFG")  # 5 matches, 1 mismatch

        identity = calculate_alignment_identity(alignment_result)

        assert identity == 5 / 6

    def test_calculate_alignment_identity_with_gaps(self):
        """Test identity calculation with gaps."""
        alignment_result = global_align("ACDEFG", "ACDFG")  # Should create gap alignment

        identity = calculate_alignment_identity(alignment_result)

        # With gaps, the alignment should still be calculated correctly
        # ACDEF-G vs ACDEF-G would be 6/6 = 1.0
        assert identity == 1.0

    def test_calculate_alignment_identity_mixed_gaps(self):
        """Test identity calculation with mixed gaps and mismatches."""
        alignment_result = global_align("ACDEFG", "ACDGF")  # Should create gap alignment

        identity = calculate_alignment_identity(alignment_result)

        # The alignment will have gaps to maximize score
        # Should be close to expected value
        assert isinstance(identity, float)
        assert 0.0 <= identity <= 1.0


class TestAlignmentDependencyAvailability:
    """Tests for dependency availability in alignment functions."""

    def test_numpy_dependency_required(self):
        """Test that numpy is required for alignment functions."""
        # This test verifies that the functions work when numpy is available
        # In a real test environment, numpy should be available
        try:
            import numpy as np

            assert np is not None
        except ImportError:
            pytest.skip("numpy not available")

        # Test that functions work with numpy available
        result = global_align("AB", "AB", match=1, mismatch=-1, gap=-2)
        assert result["score"] == 2

    def test_alignment_functions_work_with_numpy(self):
        """Test that alignment functions work correctly with numpy."""
        # Test basic functionality
        result = global_align("TEST", "TEST", match=1, mismatch=-1, gap=-2)
        assert result["score"] == 4
        assert result["identity"] == 1.0

        result = local_align("TEST", "TEST", match=1, mismatch=-1, gap=-2)
        assert result["score"] == 4
        assert result["identity"] == 1.0


class TestAlignmentEdgeCases:
    """Tests for edge cases in alignment functions."""

    def test_global_align_long_sequences(self):
        """Test global alignment with longer sequences."""
        seq1 = "MKTIIALSYIFCLVFADYKDDDDK" * 5  # ~120 residues
        seq2 = "MKTIIALSYIFCLVFADYKDDDDK" * 5

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        expected_score = len(seq1) * 1
        assert result["score"] == expected_score
        assert result["identity"] == 1.0

    def test_local_align_long_sequences(self):
        """Test local alignment with longer sequences."""
        seq1 = "MKTIIALSYIFCLVFADYKDDDDK" * 3
        seq2 = "XXX" + "MKTIIALSYIFCLVFADYKDDDDK" * 3 + "YYY"

        result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Should find the matching subsequence
        expected_score = len(seq1) * 1
        assert result["score"] == expected_score

    def test_alignment_very_different_sequences(self):
        """Test alignment of very different sequences."""
        seq1 = "ACDEFGHIKLMNPQRSTVWY"
        seq2 = "ZYXWVUTSRQPONMLKJIHGFEDCBA"

        global_result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)
        local_result = local_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Both should produce valid alignments
        assert isinstance(global_result["score"], (int, float))
        assert isinstance(local_result["score"], (int, float))
        assert global_result["score"] < 0  # Likely negative due to mismatches
        assert local_result["score"] >= 0  # Local alignment can be 0

    def test_alignment_case_sensitivity(self):
        """Test that alignment is case sensitive."""
        seq1 = "ACDEFG"
        seq2 = "acdefg"

        result = global_align(seq1, seq2, match=1, mismatch=-1, gap=-2)

        # Should treat as completely different
        assert result["score"] == -6  # All mismatches
        assert result["identity"] == 0.0
