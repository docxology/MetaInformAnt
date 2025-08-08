from __future__ import annotations

from metainformant.dna import population


def test_snp_allele_frequencies_basic() -> None:
    # Two SNP sites across four individuals
    # Columns are sites, rows are individuals; 0/1 represent alleles
    genotype_matrix = [
        [0, 1],
        [0, 1],
        [1, 1],
        [1, 0],
    ]

    freqs = population.allele_frequencies(genotype_matrix)
    # First site: 2 of 4 are 1s -> freq1 = 0.5
    # Second site: 3 of 4 are 1s -> freq1 = 0.75
    assert freqs == [0.5, 0.75]


def test_observed_heterozygosity() -> None:
    # diploid genotypes encoded as pairs of alleles (0/1)
    genotypes = [
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 0),
    ]
    h_obs = population.observed_heterozygosity(genotypes)
    # 2 heterozygotes of 4 -> 0.5
    assert h_obs == 0.5


