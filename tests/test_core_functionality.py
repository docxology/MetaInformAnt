"""Test core functionality that doesn't require external dependencies.

Comprehensive tests for basic functions that should work without scipy, sklearn, etc.
Following the no-mocking policy - all tests use real implementations.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.core.hash import sha256_bytes, sha256_file

# Test core modules
from metainformant.core.io import dump_json, load_json, open_text_auto
from metainformant.core.paths import expand_and_resolve, is_within
from metainformant.core.text import normalize_whitespace, slugify
from metainformant.dna.composition import gc_skew, melting_temperature
from metainformant.dna.distances import p_distance

# Test DNA analysis functions
from metainformant.dna.fastq import gc_content

# Test ecology functions
from metainformant.ecology.community import shannon_diversity, simpson_diversity

# Test math functions that don't need external dependencies
from metainformant.math.popgen import hardy_weinberg_genotype_freqs

# Test simulation functions (those that work with standard random)
from metainformant.simulation.sequences import generate_random_dna


class TestCoreIO:
    """Test core I/O functionality."""

    def test_json_operations(self):
        """Test JSON read/write operations."""
        test_data = {
            "species": ["E. coli", "S. cerevisiae"],
            "counts": [100, 200],
            "metadata": {"experiment": "test", "version": 1.0},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = Path(f.name)

        try:
            # Write and read JSON
            dump_json(test_data, temp_path)
            loaded_data = load_json(temp_path)

            assert loaded_data == test_data
            assert loaded_data["species"] == ["E. coli", "S. cerevisiae"]
            assert loaded_data["counts"] == [100, 200]
            assert loaded_data["metadata"]["experiment"] == "test"
        finally:
            temp_path.unlink(missing_ok=True)

    def test_text_processing(self):
        """Test text processing utilities."""
        messy_text = "  \n\t Multiple   spaces   and\nlines  \t "
        cleaned = normalize_whitespace(messy_text)

        assert cleaned == "Multiple spaces and lines"
        assert isinstance(cleaned, str)
        assert len(cleaned) < len(messy_text)

        # Test slugify
        slug = slugify("Test String With Spaces!")
        assert slug == "test-string-with-spaces"
        assert isinstance(slug, str)

    def test_hash_operations(self):
        """Test byte hashing for reproducibility."""
        data1 = b"Hello, World!"
        data2 = b"Hello, World!"
        data3 = b"Hello, Universe!"

        hash1 = sha256_bytes(data1)
        hash2 = sha256_bytes(data2)
        hash3 = sha256_bytes(data3)

        # Same data should have same hash
        assert hash1 == hash2

        # Different data should have different hash
        assert hash1 != hash3
        assert hash2 != hash3

        # Hashes should be reproducible strings
        assert isinstance(hash1, str)
        assert len(hash1) == 64  # SHA-256 hex digest length

    def test_path_operations(self):
        """Test path utility functions."""
        # Test expand and resolve
        home_path = expand_and_resolve("~")
        assert home_path.is_absolute()

        # Test is_within
        parent_dir = Path.cwd()
        child_dir = parent_dir / "subdirectory"

        assert is_within(child_dir, parent_dir) is True
        assert is_within(parent_dir, child_dir) is False


class TestDNAAnalysis:
    """Test DNA sequence analysis functions."""

    def test_gc_content_calculation(self):
        """Test GC content calculation."""
        # Test normal sequences
        assert gc_content("ATGC") == 0.5
        assert gc_content("AAAA") == 0.0
        assert gc_content("GGCC") == 1.0

        # Test edge cases
        assert gc_content("") == 0.0
        assert gc_content("NNNNN") == 0.0  # Should handle ambiguous bases
        assert gc_content("ATGCNatgc") == 0.5  # Should handle mixed case and N's

    def test_gc_skew_calculation(self):
        """Test GC skew calculation."""
        # Test balanced G/C
        assert gc_skew("GC") == 0.0

        # Test G-rich sequence
        assert gc_skew("GGGC") == 0.5

        # Test C-rich sequence
        assert gc_skew("GCCC") == -0.5

        # Test edge cases
        assert gc_skew("") == 0.0
        assert gc_skew("AT") == 0.0  # No G or C

    def test_melting_temperature(self):
        """Test melting temperature calculation."""
        # Test simple sequence
        seq = "ATGC"
        tm = melting_temperature(seq)
        assert tm > 0  # Should be positive temperature
        assert isinstance(tm, float)

        # GC-rich sequences should have higher melting temperature
        gc_rich = "GGGGCCCC"
        at_rich = "AAAATTTT"

        tm_gc = melting_temperature(gc_rich)
        tm_at = melting_temperature(at_rich)

        assert tm_gc > tm_at  # GC bonds are stronger

    def test_sequence_distances(self):
        """Test sequence distance calculations."""
        seq1 = "ATGCAT"
        seq2 = "ATGCAT"
        seq3 = "ATGCTG"  # One difference

        # Same sequences should have distance 0
        assert p_distance(seq1, seq2) == 0.0

        # Different sequences should have distance > 0
        dist = p_distance(seq1, seq3)
        assert 0 < dist <= 1
        # Two differences: position 4 (C->T) and position 5 (A->G) in "ATGCAT" vs "ATGCTG"
        assert abs(dist - 2 / 6) < 1e-10  # Should be 2 differences out of 6


class TestEcologyAnalysis:
    """Test community ecology analysis functions."""

    def test_diversity_metrics(self):
        """Test diversity calculations."""
        # Simple community
        abundances = [10, 5, 3, 2]  # 4 species

        shannon = shannon_diversity(abundances)
        simpson = simpson_diversity(abundances)

        assert shannon > 0
        assert 0 <= simpson <= 1

        # Equal abundances should give maximum diversity
        equal_abundances = [5, 5, 5, 5]
        shannon_equal = shannon_diversity(equal_abundances)
        assert shannon_equal > shannon  # Should be higher

        # Single species should give minimum diversity
        single_species = [20]
        assert shannon_diversity(single_species) == 0.0
        assert simpson_diversity(single_species) == 0.0  # No diversity with single species


class TestMathPopGen:
    """Test population genetics mathematical functions."""

    def test_hardy_weinberg_genotype_freqs(self):
        """Test Hardy-Weinberg genotype frequency calculations."""
        p = 0.6  # Frequency of A allele
        q = 0.4  # Frequency of a allele (should be 1-p)

        aa_freq, ab_freq, bb_freq = hardy_weinberg_genotype_freqs(p)

        # Should return frequencies for AA, Aa, aa
        assert abs(aa_freq - p**2) < 1e-10  # AA frequency
        assert abs(ab_freq - 2 * p * q) < 1e-10  # Aa frequency
        assert abs(bb_freq - q**2) < 1e-10  # aa frequency

        # Should sum to 1
        assert abs(aa_freq + ab_freq + bb_freq - 1.0) < 1e-10


class TestSimulation:
    """Test simulation functions that don't require external dependencies."""

    def test_sequence_generation(self):
        """Test random sequence generation."""
        # Test DNA generation
        dna_seq = generate_random_dna(length=100, gc_content=0.5)

        assert len(dna_seq) == 100
        assert all(base in "ATGC" for base in dna_seq)

        # Test GC content approximation (should be close to 0.5)
        gc_count = dna_seq.count("G") + dna_seq.count("C")
        gc_observed = gc_count / len(dna_seq)
        assert 0.3 < gc_observed < 0.7  # Allow some random variation


class TestIntegration:
    """Integration tests combining multiple modules."""

    def test_file_io_with_analysis(self):
        """Test file I/O combined with analysis functions."""
        # Create test data
        test_data = {
            "sequences": [
                {"id": "seq1", "sequence": "AAAAAATTTTTTT", "species": "E_coli"},  # Low GC (0.0)
                {"id": "seq2", "sequence": "ATGCATGCATGC", "species": "B_subtilis"},  # Medium GC (0.5)
                {"id": "seq3", "sequence": "GGGCCCGGGCCC", "species": "S_aureus"},  # High GC (1.0)
            ],
            "metadata": {"experiment": "test_diversity"},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = Path(f.name)

        try:
            # 1. Write data
            dump_json(test_data, temp_path)

            # 2. Read data
            loaded_data = load_json(temp_path)

            # 3. Extract sequences and analyze
            sequences = [item["sequence"] for item in loaded_data["sequences"]]
            gc_contents = [gc_content(seq) for seq in sequences]

            # 4. Verify analysis results
            assert len(gc_contents) == 3
            assert all(0 <= gc <= 1 for gc in gc_contents)

            # Each sequence should have its own distinct GC content
            assert len(set(gc_contents)) >= 2  # At least 2 different values

        finally:
            temp_path.unlink(missing_ok=True)


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_empty_input_handling(self):
        """Test functions handle empty inputs gracefully."""
        # Empty sequences
        assert gc_content("") == 0.0
        assert gc_skew("") == 0.0

        # Empty abundance data - test graceful handling
        try:
            shannon_diversity([])  # Might raise error or return special value
        except (ValueError, ZeroDivisionError):
            pass  # Acceptable to raise error for empty data

        try:
            p_distance("", "")  # Empty sequences
        except (ValueError, ZeroDivisionError):
            pass  # Acceptable to raise error

    def test_invalid_input_handling(self):
        """Test handling of invalid inputs."""
        from metainformant.core.errors import ValidationError
        
        # Invalid GC content
        with pytest.raises(ValidationError):
            generate_random_dna(length=100, gc_content=1.5)  # > 1.0

        with pytest.raises(ValidationError):
            generate_random_dna(length=100, gc_content=-0.1)  # < 0.0

        # Invalid length - negative length should raise error
        with pytest.raises((ValueError, ValidationError)):
            generate_random_dna(length=-1)


class TestReproducibility:
    """Test that functions produce reproducible results."""

    def test_deterministic_functions(self):
        """Test that deterministic functions are reproducible."""
        # Same input should always give same output
        seq = "ATGCGATCGATCG"

        # Test multiple calls
        gc1 = gc_content(seq)
        gc2 = gc_content(seq)
        assert gc1 == gc2

        skew1 = gc_skew(seq)
        skew2 = gc_skew(seq)
        assert skew1 == skew2

    def test_random_functions_with_seed(self):
        """Test that random functions are reproducible with seed."""
        import random

        # Use explicit RNG instances for reproducibility
        rng1 = random.Random(42)
        seq1 = generate_random_dna(length=100, gc_content=0.5, rng=rng1)

        rng2 = random.Random(42)
        seq2 = generate_random_dna(length=100, gc_content=0.5, rng=rng2)

        assert seq1 == seq2  # Should be identical with same seed

        # Different seeds should give different results
        np.random.seed(123)
        seq3 = generate_random_dna(length=100, gc_content=0.5)

        assert seq1 != seq3  # Should be different with different seed
