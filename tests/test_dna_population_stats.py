from __future__ import annotations

from metainformant.dna import population


def test_nucleotide_diversity_two_sequences() -> None:
    seqs = ["AAAA", "AAAT"]
    pi = population.nucleotide_diversity(seqs)
    assert abs(pi - 0.25) < 1e-9


def test_tajimas_d_no_segregating_sites_is_zero() -> None:
    seqs = ["AAAA", "AAAA", "AAAA"]
    d = population.tajimas_d(seqs)
    assert d == 0.0


def test_fst_fixed_differences_is_one() -> None:
    pop1 = ["AAAA", "AAAA"]
    pop2 = ["TTTT", "TTTT"]
    fst = population.hudson_fst(pop1, pop2)
    assert abs(fst - 1.0) < 1e-9
