"""Tests for DNA mutation functions."""
from metainformant.dna.variation import mutations
from metainformant.dna.alignment import distances


def test_generate_point_mutations() -> None:
    """Test generating point mutations."""
    s = "AAAAAAAA"
    # Function takes num_mutations and mutation_rate (not seed)
    result = mutations.generate_point_mutations(s, num_mutations=2, mutation_rate=0.5)
    assert isinstance(result, str)
    assert len(result) == len(s)


def test_generate_point_mutations_output() -> None:
    """Test that point mutations produces valid output."""
    s = "ATCGATCG"
    result = mutations.generate_point_mutations(s, num_mutations=1, mutation_rate=1.0)
    assert isinstance(result, str)
    assert len(result) == len(s)
    # Result should only contain valid bases
    assert all(c in "ACGT" for c in result)


def test_hamming_distance() -> None:
    """Test Hamming distance calculation."""
    s1 = "ACGTACGT"
    s2 = "ACTTACGA"  # 2 differences at positions 2 and 7
    dist = distances.hamming_distance(s1, s2)
    assert dist == 2


def test_hamming_distance_identical() -> None:
    """Test Hamming distance for identical sequences."""
    s = "ACGTACGT"
    dist = distances.hamming_distance(s, s)
    assert dist == 0
