from metainformant.dna import codon


def test_codon_usage_counts_and_freqs() -> None:
    seq = "ATGAAATTTGGGCCC"
    counts = codon.codon_counts(seq)
    assert counts["ATG"] == 1
    assert counts.get("TAA", 0) == 0
    freqs = codon.codon_frequencies(seq)
    assert abs(sum(freqs.values()) - 1.0) < 1e-9
    # Only complete codons should be counted
    seq2 = "ATGAAAAT"  # 2 full codons + partial
    counts2 = codon.codon_counts(seq2)
    assert sum(counts2.values()) == 2
