from __future__ import annotations

from metainformant.dna import sequences


def test_reverse_complement_and_gc_content() -> None:
    s = "ACGTACGTTT"
    rc = sequences.reverse_complement(s)
    assert rc == "AAACGTACGT"  # reverse complement of s
    gc = sequences.gc_content(s)
    assert abs(gc - 0.4) < 1e-9


def test_kmer_counts_and_freqs() -> None:
    s = "ATATAT"
    counts = sequences.kmer_counts(s, k=2)
    assert counts["AT"] == 3 and counts["TA"] == 2
    freqs = sequences.kmer_frequencies(s, k=2)
    assert abs(freqs["AT"] - 3 / 5) < 1e-9
