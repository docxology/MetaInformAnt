"""calculate_ld_decay on sequences with missing data, and the F_ST of populations without comparable sites.

Zero-mocks. The vectorized LD decay is checked against a plain loop over site pairs on random alignments, and
against hand-computed values, including the case the old code got wrong: sequences with missing bases at
different sites used to be paired by their position in a list of the calls that remained.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from metainformant.dna.population import analysis, core


def _reference_ld_decay(sequences: list[str], max_distance: int = 0) -> list[tuple[int, float]]:
    """Mean r-squared by distance over the pairs of biallelic sites, scored on the sequences valid at both."""
    length = len(sequences[0])
    sites = [
        pos for pos in range(length) if len({seq[pos].upper() for seq in sequences if seq[pos].upper() in "ACGT"}) == 2
    ]
    by_distance: dict[int, list[float]] = {}
    for a, left in enumerate(sites):
        for right in sites[a + 1 :]:
            if max_distance > 0 and right - left > max_distance:
                continue
            pairs = [
                (seq[left].upper(), seq[right].upper())
                for seq in sequences
                if seq[left].upper() in "ACGT" and seq[right].upper() in "ACGT"
            ]
            if len(pairs) < 2:
                continue
            first = [x for x, _ in pairs]
            second = [y for _, y in pairs]
            major_a = max(sorted(set(first)), key=first.count)
            major_b = max(sorted(set(second)), key=second.count)
            p_a = sum(x != major_a for x in first) / len(pairs)
            p_b = sum(y != major_b for y in second) / len(pairs)
            if p_a in (0.0, 1.0) or p_b in (0.0, 1.0):
                continue
            p_ab = sum(x != major_a and y != major_b for x, y in pairs) / len(pairs)
            d = p_ab - p_a * p_b
            by_distance.setdefault(right - left, []).append(d * d / (p_a * (1 - p_a) * p_b * (1 - p_b)))
    return [(distance, sum(values) / len(values)) for distance, values in sorted(by_distance.items())]


def test_sequences_with_missing_bases_at_different_sites_are_not_paired_by_list_position() -> None:
    """Sites 0 and 1 are perfectly correlated where both are called (sequences 0-3); sequences 4 and 5 miss one.

    The old code paired the five calls at site 0 with the five at site 1 by their list positions, which matched
    sequence 5's T at site 0 with sequence 4's A at site 1 and reported r-squared 0.64.
    """
    sequences = ["AA", "AA", "TT", "TT", "NA", "TN"]
    assert analysis.calculate_ld_decay(sequences) == [(1, pytest.approx(1.0))]


def test_a_pair_needs_two_sequences_called_at_both_sites_and_both_alleles_among_them() -> None:
    assert analysis.calculate_ld_decay(["AN", "NA", "TN", "NT"]) == []  # no sequence is called at both sites
    # Only sequences 0 and 1 are called at both: site 1 shows one allele among them, so the pair is skipped
    assert analysis.calculate_ld_decay(["AA", "TA", "AN", "TN", "NT", "NA"]) == []


@pytest.mark.parametrize("seed", range(4))
def test_the_vectorized_ld_decay_matches_a_plain_loop_on_random_alignments_with_missing_data(seed: int) -> None:
    generator = random.Random(seed)
    alleles = [generator.sample("ACGT", 2) for _ in range(40)]  # a biallelic alignment ...
    sequences = [
        "".join(
            generator.choice("NN-") if generator.random() < 0.12 else generator.choice(pair)  # ... with 12 % missing
            for pair in alleles
        )
        for _ in range(14)
    ]
    for max_distance in (0, 7):
        got = analysis.calculate_ld_decay(sequences, max_distance=max_distance)
        expected = _reference_ld_decay(sequences, max_distance=max_distance)
        assert len(expected) >= (6 if max_distance else 20)  # real work was compared
        assert [distance for distance, _ in got] == [distance for distance, _ in expected]
        assert [value for _, value in got] == pytest.approx([value for _, value in expected], rel=1e-9)


def test_the_ld_decay_of_recombining_sequences_falls_with_distance() -> None:
    """The coalescent with recombination gives r-squared that decays over the sites between SNPs."""
    from metainformant.simulation.models.popgen import simulate_recombinant_haplotypes

    length, rho_per_site = 600, 0.05
    sequences_by_seed = []
    for seed in range(6):
        sample = simulate_recombinant_haplotypes(
            24, float(length), segregating_sites=80, rho_per_site=rho_per_site, rng=seed
        )
        # One base per integer position: the site of a mutation is the integer part of its position
        sites = np.floor(sample.positions).astype(int)
        keep = np.concatenate(([True], np.diff(sites) > 0))
        matrix = sample.haplotypes[:, keep]
        rows = []
        for haplotype in matrix:
            row = ["A"] * length
            for site, derived in zip(sites[keep], haplotype):
                row[site] = "T" if derived else "A"
            rows.append("".join(row))
        # a polymorphic column needs both bases: ancestral A and derived T, which the columns have
        sequences_by_seed.append(rows)
    near, far = [], []
    for sequences in sequences_by_seed:
        decay = analysis.calculate_ld_decay(sequences)
        near += [r2 for distance, r2 in decay if distance <= 40]
        far += [r2 for distance, r2 in decay if distance >= 300]
    assert np.mean(near) > 1.5 * np.mean(far)  # measured 0.18 within 40 sites against 0.07 beyond 300


def test_unaligned_sequences_are_rejected_instead_of_read_past_their_end() -> None:
    with pytest.raises(ValueError, match="aligned"):
        analysis.calculate_ld_decay(["ACGT", "ACG"])
    assert analysis.calculate_ld_decay(["ACGT"]) == [] and analysis.calculate_ld_decay([]) == []


def test_hudson_fst_without_a_site_both_populations_call_is_undefined_not_zero() -> None:
    """All-ambiguous sequences used to be reported as identical populations (F_ST of 0.0)."""
    with pytest.raises(ValueError, match="undefined"):
        core.hudson_fst(["NNNN", "NNNN"], ["ACGT", "ACGA"])
    with pytest.raises(ValueError, match="undefined"):
        core.hudson_fst(["AANN", "AANN"], ["NNCC", "NNCC"])  # no column is called in both
    assert core.hudson_fst(["AAAA"], ["AAAA"]) == pytest.approx(0.0)  # identical fixed populations: still 0
    assert core.hudson_fst(["AAAA"], ["GGGG"]) == pytest.approx(1.0)  # fixed differences: still 1
