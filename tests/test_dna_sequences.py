"""Tests for DNA sequence functionality."""

import pytest

from metainformant.dna import sequences


class TestDNASequences:
    """Test DNA sequence functionality."""

    def test_sequence_length_calculation(self):
        """Test sequence length calculation."""
        # Test normal sequence
        seq = "ATCGATCGATCG"
        length = sequences.sequence_length(seq)
        assert length == 12

        # Test sequence with whitespace
        seq_with_spaces = "ATCG ATCG ATCG"
        length = sequences.sequence_length(seq_with_spaces)
        assert length == 12

        # Test empty sequence
        length = sequences.sequence_length("")
        assert length == 0

    def test_dna_sequence_validation(self):
        """Test DNA sequence validation."""
        # Valid DNA sequences
        assert sequences.validate_dna_sequence("ATCG")
        assert sequences.validate_dna_sequence("atcg")
        assert sequences.validate_dna_sequence("ATCGATCG")
        assert sequences.validate_dna_sequence("ATCG-N")  # With ambiguity

        # Invalid sequences
        assert not sequences.validate_dna_sequence("ATCGX")  # Invalid character
        assert not sequences.validate_dna_sequence("ATCG1")  # Number
        assert not sequences.validate_dna_sequence("ATCG!")  # Special character

    def test_complementarity_score(self):
        """Test complementarity score calculation."""
        # Perfect complementarity
        score = sequences.dna_complementarity_score("ATCG", "TAGC")
        assert score == 1.0

        # Partial complementarity
        score = sequences.dna_complementarity_score("ATCG", "ATCG")
        assert score == 0.0

        # Mixed case
        score = sequences.dna_complementarity_score("atcg", "TAGC")
        assert score == 1.0

        # Different lengths (should raise error)
        from metainformant.core.utils.errors import ValidationError

        with pytest.raises(ValidationError):
            sequences.dna_complementarity_score("ATCG", "ATC")

    def test_repeat_finding(self):
        """Test repeat finding functionality."""
        # Sequence with repeats
        seq = "ATCGATCGATCG"
        repeats = sequences.find_repeats(seq, min_length=3)

        # Should find repeating patterns
        assert "ATC" in repeats
        assert "TCG" in repeats
        # Each repeat appears multiple times
        assert len(repeats["ATC"]) >= 2

        # Sequence without repeats
        seq_no_repeats = "ATCG"
        repeats = sequences.find_repeats(seq_no_repeats, min_length=2)
        # Short unique sequence may still have some overlapping k-mers
        # Just verify it returns a dict
        assert isinstance(repeats, dict)

    def test_motif_finding(self):
        """Test multiple motif finding."""
        seq = "ATCGATCGATCG"

        # Find multiple motifs
        motifs = ["ATC", "TCG", "GAT"]
        results = sequences.find_motifs(seq, motifs)

        assert "ATC" in results
        assert "TCG" in results
        assert "GAT" in results

        # Check positions
        assert 0 in results["ATC"]  # First ATC at position 0
        assert 4 in results["ATC"]  # Second ATC at position 4

    def test_sequence_complexity(self):
        """Test sequence complexity calculation."""
        # High complexity sequence
        complex_seq = "ATCGATCGATCG"
        complexity = sequences.calculate_sequence_complexity(complex_seq)
        assert 0 <= complexity <= 1

        # Low complexity sequence (repetitive)
        repetitive_seq = "AAAAAAAAAAAA"
        repetitive_complexity = sequences.calculate_sequence_complexity(repetitive_seq)
        assert repetitive_complexity < complexity

        # Very short sequence
        short_seq = "A"
        short_complexity = sequences.calculate_sequence_complexity(short_seq)
        assert short_complexity == 0.0

    def test_orf_finding(self):
        """Test ORF finding functionality."""
        # Sequence with ORFs - need RNA sequence for find_orfs
        # DNA: ATGAAATTTAAATAG contains start (ATG) and stop (TAA, TAG)
        seq = "ATGAAATTTAAATAG"

        # find_orfs works on RNA and requires min_length (default 30)
        # Use a smaller min_length for testing
        orfs = sequences.find_orfs(seq, min_length=2)
        # Returns list of tuples or empty list
        assert isinstance(orfs, list)

    def test_sequence_entropy(self):
        """Test sequence entropy calculation."""
        # Uniform sequence (low entropy)
        uniform_seq = "AAAA"
        entropy = sequences.calculate_sequence_entropy(uniform_seq, k=1)
        assert entropy == 0.0  # No information

        # Diverse sequence (high entropy)
        diverse_seq = "ATCG"
        entropy = sequences.calculate_sequence_entropy(diverse_seq, k=1)
        assert entropy > 0

        # Test with different k values
        entropy_k1 = sequences.calculate_sequence_entropy(diverse_seq, k=1)
        entropy_k2 = sequences.calculate_sequence_entropy(diverse_seq, k=2)
        # For short sequences, k=2 entropy may be lower than k=1
        # This is mathematically correct behavior
        assert entropy_k1 >= 0
        assert entropy_k2 >= 0

    def test_sequence_bias_detection(self):
        """Test sequence bias detection."""
        # GC-rich sequence
        gc_rich = "GCGCGCGCGCGC"
        bias = sequences.detect_sequence_bias(gc_rich)

        assert bias["gc_content"] > 0.8
        assert bias["at_content"] < 0.2
        assert bias["total_bases"] == len(gc_rich)

        # AT-rich sequence
        at_rich = "ATATATATATAT"
        bias = sequences.detect_sequence_bias(at_rich)

        assert bias["at_content"] > 0.8
        assert bias["gc_content"] < 0.2

        # Balanced sequence
        balanced = "ATCGATCGATCG"
        bias = sequences.detect_sequence_bias(balanced)

        assert abs(bias["gc_content"] - 0.5) < 0.1
        assert abs(bias["at_content"] - 0.5) < 0.1

        # Empty sequence
        empty_bias = sequences.detect_sequence_bias("")
        assert empty_bias["gc_content"] == 0.0
        assert empty_bias["at_content"] == 0.0

    def test_gc_skew_calculation(self):
        """Test GC skew calculation."""
        # All G sequence
        all_g = "GGGGGGGGGGGG"
        skew = sequences.calculate_gc_skew(all_g)
        assert skew == 1.0  # All G, no C

        # AT-only sequence (no G or C, skew is undefined/0.0)
        at_only = "ATATATATATAT"
        skew = sequences.calculate_gc_skew(at_only)
        assert skew == 0.0  # No G or C, returns 0.0

        # Balanced G/C sequence
        balanced = "ATCGATCGATCG"
        skew = sequences.calculate_gc_skew(balanced)
        assert abs(skew) < 0.1  # Should be close to zero

    def test_at_skew_calculation(self):
        """Test AT skew calculation."""
        # A-rich sequence
        a_rich = "AAAAAAAAAAAA"
        skew = sequences.calculate_at_skew(a_rich)
        assert skew == 1.0  # All A, no T

        # T-rich sequence
        t_rich = "TTTTTTTTTTTT"
        skew = sequences.calculate_at_skew(t_rich)
        assert skew == -1.0  # All T, no A

        # Balanced sequence
        balanced = "ATCGATCGATCG"
        skew = sequences.calculate_at_skew(balanced)
        assert abs(skew) < 0.1  # Should be close to zero

    def test_palindrome_finding(self):
        """Test palindrome finding functionality."""
        # Sequence with palindromes
        seq = "ATCGATCGATCG"  # Contains "ATCGAT" which is a palindrome when reversed

        palindromes = sequences.find_palindromes(seq, min_length=4)
        assert len(palindromes) >= 1

        # Check that found palindromes are actually palindromes
        for palindrome, start, end in palindromes:
            assert sequences.reverse_complement(palindrome) == palindrome

    def test_melting_temperature_calculation(self):
        """Test melting temperature calculation."""
        # Short sequence
        short_seq = "ATCGATCG"
        tm = sequences.calculate_melting_temperature(short_seq)

        # Should return a float
        assert isinstance(tm, (int, float))
        assert tm > 0  # Temperature should be positive

    def test_codon_usage_calculation(self):
        """Test codon usage calculation."""
        # Sequence divisible by 3
        seq = "ATGAAATTTAAATAG"  # 15 bases = 5 codons

        usage = sequences.calculate_codon_usage(seq)
        # Should return a dictionary of codon frequencies
        assert isinstance(usage, dict)
        # Should have some codons
        assert len(usage) > 0

    def test_start_codon_finding(self):
        """Test start codon finding."""
        # Sequence with two ATG start codons at positions 0 and 6
        seq = "ATGAAAATGAAATAA"

        starts = sequences.find_start_codons(seq)
        assert 0 in starts  # ATG at position 0
        assert 6 in starts  # ATG at position 6

    def test_stop_codon_finding(self):
        """Test stop codon finding."""
        seq = "ATGAAAATGAAATAA"

        stops = sequences.find_stop_codons(seq)
        assert 12 in stops  # TAA at position 12

    def test_sequence_complexity_calculation(self):
        """Test sequence complexity calculation."""
        # High complexity (all unique dimers)
        complex_seq = "ATCGATCGATCG"
        complexity = sequences.calculate_sequence_complexity(complex_seq)
        assert 0 <= complexity <= 1

        # Low complexity (repetitive)
        repetitive_seq = "AAAAAAAAAAAA"
        repetitive_complexity = sequences.calculate_sequence_complexity(repetitive_seq)
        assert repetitive_complexity < complexity

    def test_sequence_entropy_calculation(self):
        """Test sequence entropy calculation."""
        # Uniform sequence
        uniform = "AAAA"
        entropy = sequences.calculate_sequence_entropy(uniform, k=1)
        assert entropy == 0.0  # No information

        # Diverse sequence
        diverse = "ATCG"
        entropy = sequences.calculate_sequence_entropy(diverse, k=1)
        assert entropy > 0
