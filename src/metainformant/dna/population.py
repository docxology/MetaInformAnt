from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple


def allele_frequencies(genotype_matrix: Sequence[Sequence[int]]) -> List[float]:
    """Compute frequency of allele '1' per site.

    genotype_matrix: rows are individuals, columns are sites; values 0/1.
    Returns list of frequencies per site.
    """
    if not genotype_matrix:
        return []
    num_sites = len(genotype_matrix[0])
    site_totals = [0] * num_sites
    for row in genotype_matrix:
        for j, val in enumerate(row):
            site_totals[j] += int(val == 1)
    n_individuals = len(genotype_matrix)
    return [total / n_individuals for total in site_totals]


def observed_heterozygosity(genotypes: Iterable[Tuple[int, int]]) -> float:
    """Proportion of heterozygous individuals among diploid genotypes.

    genotypes: iterable of tuples (a1, a2) with alleles 0/1
    """
    genotypes_list = list(genotypes)
    if not genotypes_list:
        return 0.0
    hetero = sum(1 for a1, a2 in genotypes_list if a1 != a2)
    return hetero / len(genotypes_list)


