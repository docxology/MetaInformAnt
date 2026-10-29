"""Advanced population genetics analysis utilities.

This module provides comprehensive tools for population genetics analysis
including demographic inference, selection tests, and population structure analysis.
"""

from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging
from metainformant.dna.annotation.gene_finding import GENETIC_CODES

logger = logging.get_logger(__name__)

#: How the p-value of a neutrality test is obtained: the normal approximation or a simulated null.
_P_VALUE_METHODS = ("normal", "simulation")
#: Default number of coalescent replicates behind a simulated p-value.
DEFAULT_NULL_SIMULATIONS = 2000


def calculate_fst(population1: List[str], population2: List[str]) -> float:
    """Calculate F_ST between two populations from aligned sequences.

    Uses the heterozygosity-based per-site estimator: at each polymorphic
    site ``F_ST = (H_T - H_S) / H_T``, with ``H_T`` the heterozygosity of
    the pooled populations and ``H_S`` the mean within-population
    heterozygosity, averaged over polymorphic sites (Nei's G_ST style
    estimator). This is not Hudson's estimator; for Hudson's F_ST based on
    between- versus within-population pairwise diversity use
    :func:`metainformant.dna.population.core.hudson_fst`.

    Args:
        population1: List of DNA sequences from population 1
        population2: List of DNA sequences from population 2

    Returns:
        F_ST value (0.0 to 1.0)

    Raises:
        ValueError: If populations have incompatible sequences

    Example:
        >>> pop1 = ["ATCG", "ATCG"]
        >>> pop2 = ["ATCG", "GCTA"]
        >>> fst = calculate_fst(pop1, pop2)
        >>> 0.0 <= fst <= 1.0
        True
    """
    if not population1 or not population2:
        return 0.0

    # Check sequence lengths
    seq_len = len(population1[0])
    if not all(len(seq) == seq_len for seq in population1 + population2):
        raise ValueError("All sequences must have the same length")

    total_sites = 0
    numerator_sum = 0.0

    for pos in range(seq_len):
        # Keep only unambiguous bases (ATCG, case-insensitive) so gaps and
        # ambiguous characters are not treated as alleles
        alleles_pop1 = [seq[pos].upper() for seq in population1 if seq[pos].upper() in "ATCG"]
        alleles_pop2 = [seq[pos].upper() for seq in population2 if seq[pos].upper() in "ATCG"]

        if not alleles_pop1 or not alleles_pop2:
            continue  # No usable data at this site

        # Count alleles
        allele_counts: Dict[str, int] = {}
        for allele in alleles_pop1 + alleles_pop2:
            allele_counts[allele] = allele_counts.get(allele, 0) + 1

        if len(allele_counts) < 2:
            continue  # Monomorphic site

        # Calculate allele frequencies
        n1, n2 = len(alleles_pop1), len(alleles_pop2)
        total_n = n1 + n2

        freq_pop1 = {}
        freq_pop2 = {}
        total_freq = {}

        for allele in allele_counts:
            count_pop1 = alleles_pop1.count(allele)
            count_pop2 = alleles_pop2.count(allele)

            freq_pop1[allele] = count_pop1 / n1 if n1 > 0 else 0
            freq_pop2[allele] = count_pop2 / n2 if n2 > 0 else 0
            total_freq[allele] = (count_pop1 + count_pop2) / total_n

        # Calculate F_ST contribution for this site
        ht = 1 - sum(f**2 for f in total_freq.values())  # Total heterozygosity
        hs = (sum(f**2 for f in freq_pop1.values()) + sum(f**2 for f in freq_pop2.values())) / 2
        hs = 1 - hs  # Average within-population heterozygosity

        if ht > 0:
            fst_site = (ht - hs) / ht
            numerator_sum += fst_site
            total_sites += 1

    if total_sites == 0:
        return 0.0

    # Mean F_ST across polymorphic sites with usable data
    return numerator_sum / total_sites


_SELECTION_METHODS = ("tajima_d", "fu_li_d", "mk_test")


def detect_selection(
    sequences: List[str],
    method: str = "tajima_d",
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Dict[str, float]:
    """Detect signatures of natural selection in population data.

    ``"tajima_d"`` and ``"fu_li_d"`` (Fu and Li's D*) return the statistic and its two-sided p-value, from
    :func:`calculate_tajima_d` and :func:`calculate_fu_li_d`; ``p_value_method`` selects the null behind the
    p-value as there: ``"normal"`` (default) is the normal approximation, ``"simulation"`` a Monte Carlo
    p-value against ``n_simulations`` neutral coalescent samples conditioned on the observed number of
    mutations, seeded by ``random_state`` (see :func:`calculate_tajima_d` for the details). Fewer than 4
    sequences, empty input included, give the statistic 0.0 and the p-value 1.0, as there. ``"mk_test"``
    returns the McDonald-Kreitman ``alpha`` and ``omega`` (:func:`mcdonald_kreitman_test`), whose
    significance is Fisher's exact test (:func:`mcdonald_kreitman_contingency`), not a normal or simulated
    null: a ``p_value_method`` other than ``"normal"`` is refused for it rather than ignored.

    Args:
        sequences: List of DNA sequences from a population
        method: Selection detection method ("tajima_d", "fu_li_d", "mk_test")
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (an integer, a ``random.Random`` or a
            ``numpy.random.Generator``); ``None`` draws fresh entropy

    Returns:
        Dictionary with selection statistics: ``{method: statistic, "p_value": p}`` for the first two
        methods, ``{"alpha": ..., "omega": ...}`` for ``"mk_test"``

    Raises:
        ValueError: If ``method`` or ``p_value_method`` is unknown, ``n_simulations`` is not a positive
            integer, the sequences are not aligned, or ``p_value_method="simulation"`` is combined with
            ``"mk_test"``

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> selection = detect_selection(seqs, method="tajima_d")
        >>> "tajima_d" in selection
        True
    """
    if method not in _SELECTION_METHODS:
        raise ValueError(f"Unknown selection detection method: {method}")
    _check_p_value_options(p_value_method, n_simulations)

    results = {}

    if method == "tajima_d":
        tajima_d, p_value = calculate_tajima_d(
            sequences, p_value_method=p_value_method, n_simulations=n_simulations, random_state=random_state
        )
        results["tajima_d"] = tajima_d
        results["p_value"] = p_value

    elif method == "fu_li_d":
        fu_li_d, p_value = calculate_fu_li_d(
            sequences, p_value_method=p_value_method, n_simulations=n_simulations, random_state=random_state
        )
        results["fu_li_d"] = fu_li_d
        results["p_value"] = p_value

    else:
        if p_value_method != "normal":
            raise ValueError(
                f"p_value_method={p_value_method!r} does not apply to the McDonald-Kreitman test, whose "
                "significance is Fisher's exact test (mcdonald_kreitman_contingency)"
            )
        alpha, omega = mcdonald_kreitman_test(sequences)
        results["alpha"] = alpha
        results["omega"] = omega

    return results


def _two_sided_normal_p(statistic: float) -> float:
    """Two-sided p-value of a standardized statistic under the normal approximation."""
    return math.erfc(abs(statistic) / math.sqrt(2.0))


def _check_p_value_options(p_value_method: str, n_simulations: int) -> None:
    """Reject an unknown ``p_value_method`` or a non-positive number of simulations."""
    if p_value_method not in _P_VALUE_METHODS:
        raise ValueError(f"p_value_method must be one of {_P_VALUE_METHODS}, got {p_value_method!r}")
    if isinstance(n_simulations, bool) or not isinstance(n_simulations, (int, np.integer)) or n_simulations < 1:
        raise ValueError(f"n_simulations must be a positive integer, got {n_simulations!r}")


def _neutrality_p_value(
    sequences: List[str],
    test: str,
    statistic: float,
    p_value_method: str,
    n_simulations: int,
    random_state: int | random.Random | np.random.Generator | None,
) -> float:
    """Two-sided p-value of one SFS test statistic by the chosen method."""
    if p_value_method == "simulation":
        # The coalescent kernel lives in the simulation domain; the integration adapter is the only
        # dna module that reaches it, and it is imported here so that dna works without simulation loaded.
        from ..integration.neutrality_null import simulated_neutrality_tests

        simulated = simulated_neutrality_tests(
            sequences, (test,), n_simulations=n_simulations, random_state=random_state
        )
        return simulated[test][1]
    return _two_sided_normal_p(statistic)


def calculate_tajima_d(
    sequences: List[str],
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Tuple[float, float]:
    """Calculate Tajima's D and its two-sided p-value.

    See :func:`metainformant.dna.population.core.tajimas_d`. ``p_value_method`` selects the null
    distribution behind the p-value:

    * ``"normal"`` (default): the normal approximation, which is rough for small samples.
    * ``"simulation"``: a Monte Carlo p-value against ``n_simulations`` samples drawn from the standard
      neutral coalescent (constant size, infinite sites, no recombination) with the same number of
      sequences. Every replicate is conditioned on the number of mutations observed in the sample
      rather than on an estimate such as ``theta = S / a_n``, with the genealogy drawn unweighted
      (Hudson's ``ms -s``; equivalently ``theta`` is integrated out with the scale-invariant prior). So
      no mutation rate is assumed or estimated, no replicate is monomorphic, and the replicate
      statistics are standardized with the same variance as the observed one. The p-value is twice the
      smaller one-sided tail, ``(1 + #{null at least as extreme}) / (n_simulations + 1)``, capped at 1,
      so it is never below ``2 / (n_simulations + 1)``. For neutral samples with about 60 mutations it
      is uniform to within sampling error; with very few mutations the null is discrete and the
      p-value conservative. Replicates are vectorized; ``random_state`` seeds them.

    Args:
        sequences: List of aligned DNA sequences
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed (an integer, a ``random.Random`` or a ``numpy.random.Generator``) of the
            simulated replicates; ``None`` draws fresh entropy, so pass a seed for reproducible p-values

    Returns:
        Tuple of (tajima_d, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Raises:
        ValueError: If ``p_value_method`` is unknown, ``n_simulations`` is not a positive integer, or
            the sequences are not aligned

    Example:
        >>> seqs = ["AAAA", "AAAA", "AACT", "AAAT"]
        >>> tajima_d, p_value = calculate_tajima_d(seqs)
        >>> round(tajima_d, 4)
        0.5916
    """
    _check_p_value_options(p_value_method, n_simulations)
    if len(sequences) < 4:
        return 0.0, 1.0  # Not enough sequences for meaningful test

    from . import core as population

    tajima_d = population.tajimas_d(sequences)
    return tajima_d, _neutrality_p_value(sequences, "tajima_d", tajima_d, p_value_method, n_simulations, random_state)


def calculate_fu_li_d(
    sequences: List[str],
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Tuple[float, float]:
    """Calculate Fu and Li's D* (no outgroup) and its two-sided p-value.

    See :func:`metainformant.dna.population.core.fu_and_li_d_star_from_sequences`. The p-value is
    obtained as in :func:`calculate_tajima_d`: by the normal approximation (default) or, with
    ``p_value_method="simulation"``, against ``n_simulations`` neutral coalescent samples conditioned on
    the observed number of mutations.

    Args:
        sequences: List of aligned DNA sequences
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (integer, ``random.Random`` or
            ``numpy.random.Generator``); ``None`` draws fresh entropy

    Returns:
        Tuple of (fu_li_d_star, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Raises:
        ValueError: If ``p_value_method`` is unknown, ``n_simulations`` is not a positive integer, or
            the sequences are not aligned

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTT"]
        >>> fu_li_d, p_value = calculate_fu_li_d(seqs)
        >>> isinstance(fu_li_d, float)
        True
    """
    _check_p_value_options(p_value_method, n_simulations)
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    d_star = population.fu_and_li_d_star_from_sequences(sequences)
    return d_star, _neutrality_p_value(sequences, "fu_li_d_star", d_star, p_value_method, n_simulations, random_state)


def _fisher_exact_two_sided_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for the 2x2 table ``[[a, b], [c, d]]``.

    Uses :func:`scipy.stats.fisher_exact` when SciPy is installed; otherwise
    computes the exact hypergeometric two-sided p-value using SciPy's
    definition: the sum of the probabilities of all tables with the same
    margins whose probability does not exceed the observed table's.
    """
    if min(a, b, c, d) < 0:
        raise ValueError("Counts must be non-negative")
    if a + b == 0 or c + d == 0 or a + c == 0 or b + d == 0:
        return 1.0  # A degenerate margin leaves nothing to test

    scipy_stats: Any  # None when SciPy is unavailable
    try:
        from scipy import stats as scipy_stats
    except ImportError:
        scipy_stats = None

    if scipy_stats is not None:
        return float(scipy_stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])

    n = a + b + c + d
    row1 = a + b
    col1 = a + c

    def _prob(k: int) -> float:
        return (math.comb(row1, k) * math.comb(n - row1, col1 - k)) / math.comb(n, col1)

    lo = max(0, col1 - (n - row1))
    hi = min(row1, col1)
    observed = _prob(a)
    tolerance = observed * 1e-9
    return min(
        1.0,
        sum(_prob(k) for k in range(lo, hi + 1) if _prob(k) <= observed + tolerance),
    )


def mcdonald_kreitman_contingency(sequences: List[str]) -> Dict[str, Any]:
    """Build the McDonald-Kreitman 2x2 contingency table from coding sequences.

    The final sequence is the outgroup (divergence reference); all preceding
    sequences form the ingroup population sample (polymorphism). Sequences
    must be aligned coding sequences of equal length; bases beyond the last
    complete codon are ignored, and codon columns containing ambiguous bases
    are skipped. A site polymorphic within the ingroup is counted only as
    polymorphism (standard site-based MK counting), never additionally as
    divergence.

    Returns:
        Dict with the synonymous/nonsynonymous polymorphism (``Ps``, ``Pn``)
        and divergence (``Ds``, ``Dn``) counts, the neutrality index
        ``NI = (Pn / Ps) / (Dn / Ds)`` (``None`` when the ratio is
        undefined), the two-sided Fisher exact p-value for the table
        ``[[Pn, Ps], [Dn, Ds]]``, and the ``alpha`` / ``omega`` summary
        returned by :func:`mcdonald_kreitman_test`.
    """
    if len(sequences) < 3:
        # Fewer than two ingroup sequences plus one outgroup carries no
        # joint polymorphism/divergence information: report an empty table.
        return {
            "Ps": 0,
            "Pn": 0,
            "Ds": 0,
            "Dn": 0,
            "neutral_ratio": None,
            "alpha": 0.0,
            "omega": 0.0,
            "fisher_p": 1.0,
        }

    seq_len = len(sequences[0])
    if any(len(seq) != seq_len for seq in sequences):
        raise ValueError("Sequences must be aligned (equal length)")

    ingroup = [seq.upper() for seq in sequences[:-1]]
    outgroup = sequences[-1].upper()
    code = GENETIC_CODES.get(1, GENETIC_CODES[1])
    bases = frozenset("ATCG")

    ps = pn = ds = dn = 0
    for i in range(0, seq_len - 2, 3):
        in_codons = [seq[i : i + 3] for seq in ingroup]
        out_codon = outgroup[i : i + 3]
        if any(set(codon) - bases for codon in in_codons) or set(out_codon) - bases:
            continue  # Ambiguous or gapped codon column is unusable

        distinct = set(in_codons)
        if len(distinct) > 1:
            # Polymorphic within the ingroup
            amino_acids = {code[codon] for codon in distinct}
            if len(amino_acids) > 1:
                pn += 1
            else:
                ps += 1
        elif in_codons[0] != out_codon:
            # Monomorphic in the ingroup but divergent from the outgroup
            if code[in_codons[0]] != code[out_codon]:
                dn += 1
            else:
                ds += 1

    fisher_p = _fisher_exact_two_sided_2x2(pn, ps, dn, ds)

    if ps > 0 and dn > 0:
        neutral_ratio = (pn / ps) * (ds / dn)
        alpha = 1.0 - neutral_ratio
    else:
        neutral_ratio = None
        alpha = 0.0

    omega = dn / ds if ds > 0 else 0.0

    return {
        "Ps": ps,
        "Pn": pn,
        "Ds": ds,
        "Dn": dn,
        "neutral_ratio": neutral_ratio,
        "alpha": alpha,
        "omega": omega,
        "fisher_p": fisher_p,
    }


def mcdonald_kreitman_test(sequences: List[str]) -> Tuple[float, float]:
    """Perform the McDonald-Kreitman test for selection on coding sequences.

    The last sequence is treated as the outgroup and the preceding aligned
    coding sequences as the ingroup population sample (see
    :func:`mcdonald_kreitman_contingency` for the exact counting rules).

    Returns:
        Tuple of ``(alpha, omega)``:

        - ``alpha = 1 - NI`` estimates the proportion of adaptive
          substitutions, where ``NI = (Pn/Ps) / (Dn/Ds)`` is the neutrality
          index; 0.0 when the ratio is undefined (no adaptive signal
          estimable).
        - ``omega = Dn / Ds`` is the nonsynonymous-to-synonymous divergence
          ratio; 0.0 when Ds == 0.

    Example:
        >>> seqs = ["ATGGCC", "ATGGCC", "ATGGCC"]
        >>> alpha, omega = mcdonald_kreitman_test(seqs)
        >>> isinstance(alpha, float) and isinstance(omega, float)
        True
    """
    results = mcdonald_kreitman_contingency(sequences)
    return results["alpha"], results["omega"]


def estimate_population_size(sequences: List[str], mutation_rate: float = 1e-8) -> Dict[str, float]:
    """Estimate effective population size using various methods.

    Args:
        sequences: List of DNA sequences
        mutation_rate: Per-site mutation rate

    Returns:
        Dictionary with population size estimates

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> sizes = estimate_population_size(seqs)
        >>> "theta_watterson" in sizes
        True
    """
    if not sequences:
        return {"theta_watterson": 0.0, "theta_pi": 0.0, "ne_estimate": 0.0}

    # Calculate theta estimates
    theta_w = calculate_wattersons_theta(sequences)
    theta_pi = calculate_nucleotide_diversity(sequences)

    # Estimate Ne using Watterson's theta
    # Ne = theta / (4 * mu * L) where L is sequence length
    seq_length = len(sequences[0])
    if mutation_rate > 0 and seq_length > 0:
        ne_estimate = theta_w / (4 * mutation_rate * seq_length)
    else:
        ne_estimate = 0.0

    return {
        "theta_watterson": theta_w,
        "theta_pi": theta_pi,
        "ne_estimate": ne_estimate,
        "sequence_length": seq_length,
        "sample_size": len(sequences),
    }


def detect_population_structure(sequences: List[str], k_max: int = 5) -> Dict[str, Any]:
    """Detect population structure using simple clustering.

    Args:
        sequences: List of DNA sequences
        k_max: Maximum number of populations to test

    Returns:
        Dictionary with structure analysis results

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> structure = detect_population_structure(seqs, k_max=2)
        >>> "clusters" in structure
        True
    """
    if len(sequences) < 4:
        return {"clusters": [], "k_optimal": 1}

    # Simple distance-based clustering
    # In practice, would use more sophisticated methods

    def _p_distance_simple(seq_a: str, seq_b: str) -> float:
        """Calculate p-distance between two sequences."""
        if len(seq_a) != len(seq_b):
            min_len = min(len(seq_a), len(seq_b))
            seq_a = seq_a[:min_len]
            seq_b = seq_b[:min_len]
        if len(seq_a) == 0:
            return 0.0
        mismatches = sum(1 for a, b in zip(seq_a, seq_b) if a != b)
        return mismatches / len(seq_a)

    # Calculate pairwise distances
    pairwise_dists: List[Tuple[int, int, float]] = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            dist = _p_distance_simple(sequences[i], sequences[j])
            pairwise_dists.append((i, j, dist))

    # Simple clustering based on distance threshold
    threshold = 0.1  # 10% divergence
    cluster_map: Dict[int, List[int]] = {}
    cluster_assignments: Dict[int, int] = {}

    cluster_id = 0
    for i in range(len(sequences)):
        assigned = False
        for cid, members in cluster_map.items():
            # Check if this sequence is close to cluster members
            avg_dist = sum(_p_distance_simple(sequences[i], sequences[j]) for j in members) / len(members)
            if avg_dist < threshold:
                members.append(i)
                cluster_assignments[i] = cid
                assigned = True
                break

        if not assigned:
            cluster_map[cluster_id] = [i]
            cluster_assignments[i] = cluster_id
            cluster_id += 1

    clusters = list(cluster_map.values())

    return {
        "clusters": clusters,
        "cluster_assignments": cluster_assignments,
        "k_optimal": len(clusters),
        "threshold_used": threshold,
    }


def calculate_ld_decay(sequences: List[str], max_distance: int = 0) -> List[Tuple[int, float]]:
    """Calculate linkage disequilibrium (r²) decay with physical distance.

    Computes pairwise r² between biallelic polymorphic sites and averages by distance. A site is
    biallelic when exactly two of A, C, G and T occur at it; the other characters (gaps, N, ...) are
    missing data. A pair of sites is scored on the sequences that carry an unambiguous base at both
    of them (and only when they are at least two, and each site still shows both of its alleles among
    them), so missing data never pair a sequence's allele at one site with another sequence's at the
    other.

    Args:
        sequences: List of aligned DNA sequences
        max_distance: Maximum distance to consider (0 = all)

    Returns:
        List of tuples (distance, average_r_squared), sorted by distance

    Raises:
        ValueError: If the sequences are not aligned

    Example:
        >>> seqs = ["ATCGATCG", "GCTAGCTA", "ATCGATCG", "GCTAGCTA"]
        >>> ld_decay = calculate_ld_decay(seqs)
        >>> isinstance(ld_decay, list)
        True
    """
    if len(sequences) < 2:
        return []
    seq_length = len(sequences[0])
    if any(len(seq) != seq_length for seq in sequences):
        raise ValueError("Sequences must be aligned (same length)")

    codes = np.array([np.frombuffer(seq.upper().encode("ascii", "replace"), dtype=np.uint8) for seq in sequences])
    counts = np.array([(codes == ord(base)).sum(axis=0) for base in "ACGT"])  # allele counts per site
    sites = np.flatnonzero((counts > 0).sum(axis=0) == 2)  # biallelic polymorphic sites
    if sites.size < 2:
        return []

    valid = np.isin(codes[:, sites], [ord(base) for base in "ACGT"])
    major = np.array([ord("ACGT"[index]) for index in counts[:, sites].argmax(axis=0)], dtype=np.uint8)
    present = valid.astype(float)  # sequence has an unambiguous base at the site
    minor = (valid & (codes[:, sites] != major)).astype(float)  # ... and it is not the major allele

    sums = np.zeros(seq_length)
    pairs = np.zeros(seq_length)
    for begin in range(0, sites.size, 256):  # blocks bound the memory of the pairwise matrices
        rows = slice(begin, begin + 256)
        n_both = present[:, rows].T @ present
        with np.errstate(divide="ignore", invalid="ignore"):  # pairs without data or without both alleles are dropped
            freq_a = (minor[:, rows].T @ present) / n_both
            freq_b = (present[:, rows].T @ minor) / n_both
            d = (minor[:, rows].T @ minor) / n_both - freq_a * freq_b
            r_squared = d * d / (freq_a * (1 - freq_a) * freq_b * (1 - freq_b))
        distance = sites[None, :] - sites[rows, None]
        keep = (distance > 0) & (n_both >= 2) & (freq_a > 0) & (freq_a < 1) & (freq_b > 0) & (freq_b < 1)
        if max_distance > 0:
            keep &= distance <= max_distance
        sums += np.bincount(distance[keep], weights=r_squared[keep], minlength=seq_length)
        pairs += np.bincount(distance[keep], minlength=seq_length)

    return [(int(dist), float(sums[dist] / pairs[dist])) for dist in np.flatnonzero(pairs)]


# Helper functions that delegate to existing implementations
def calculate_nucleotide_diversity(sequences: List[str]) -> float:
    """Calculate nucleotide diversity (π)."""
    from . import core as population

    return population.nucleotide_diversity(sequences)


def calculate_wattersons_theta(sequences: List[str]) -> float:
    """Calculate Watterson's θ."""
    from . import core as population

    return population.wattersons_theta(sequences)


def calculate_segregating_sites(sequences: List[str]) -> int:
    """Count segregating sites."""
    from . import core as population

    return population.segregating_sites(sequences)


def calculate_summary_statistics(
    sequences: List[str] | None = None,
    genotype_matrix: List[List[int]] | None = None,
    populations: List[int] | None = None,
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Dict[str, Any]:
    """Calculate comprehensive population genetics summary statistics.

    The neutrality tests of the sequences (Tajima's D and Fu and Li's D*) take their p-values from
    ``p_value_method`` as in :func:`calculate_tajima_d`: the normal approximation (default) or a Monte Carlo
    p-value against ``n_simulations`` neutral coalescent samples conditioned on the observed number of
    mutations, seeded by ``random_state``. Fewer than 4 sequences give the statistics 0.0 and the p-values
    1.0 there, and so here.

    Args:
        sequences: List of DNA sequences (for sequence-based statistics)
        genotype_matrix: Individuals x loci matrix of diploid 0/1/2 alternate-allele
            counts (for genotype-based statistics); other values count as missing
        populations: Population assignments for each individual. With sequences, ``fst`` is the
            heterozygosity-based F_ST (:func:`calculate_fst`) of the two populations, or the mean over all
            pairs of populations when there are more
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (an integer, a ``random.Random`` or a
            ``numpy.random.Generator``); ``None`` draws fresh entropy, so pass a seed for reproducible p-values

    Returns:
        Dictionary containing various population genetics statistics; ``fu_li_d``
        holds Fu and Li's D*, which needs no outgroup

    Raises:
        ValueError: If ``p_value_method`` is unknown, ``n_simulations`` is not a positive integer, the
            sequences are not aligned, or ``populations`` does not assign every sequence to a population

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> stats = calculate_summary_statistics(sequences=seqs)
        >>> isinstance(stats, dict)
        True
    """
    _check_p_value_options(p_value_method, n_simulations)
    p_options: Dict[str, Any] = {
        "p_value_method": p_value_method,
        "n_simulations": n_simulations,
        "random_state": random_state,
    }
    results: Dict[str, Any] = {}

    if sequences:
        logger.info(f"Calculating summary statistics for {len(sequences)} sequences")

        # Basic diversity measures
        results["nucleotide_diversity"] = calculate_nucleotide_diversity(sequences)
        results["watterson_theta"] = calculate_wattersons_theta(sequences)
        results["segregating_sites"] = calculate_segregating_sites(sequences)

        # Neutrality tests
        tajima_d, tajima_p = calculate_tajima_d(sequences, **p_options)
        results["tajima_d"] = tajima_d
        results["tajima_d_p_value"] = tajima_p

        fu_li_d, fu_li_p = calculate_fu_li_d(sequences, **p_options)
        results["fu_li_d"] = fu_li_d
        results["fu_li_d_p_value"] = fu_li_p

        # Sequence properties
        results["sequence_length"] = len(sequences[0]) if sequences else 0
        results["sample_size"] = len(sequences)

    if genotype_matrix:
        logger.info(f"Calculating summary statistics for genotype matrix with {len(genotype_matrix)} individuals")

        # Convert to numpy for easier computation
        genotypes = np.array(genotype_matrix)

        # Basic statistics
        results["num_individuals"] = len(genotypes)
        results["num_loci"] = genotypes.shape[1] if genotypes.size > 0 else 0

        # Observed heterozygosity: entries are diploid 0/1/2 alternate-allele
        # counts, so only 1 is heterozygous; other values count as missing.
        if genotypes.size > 0:
            called = np.isin(genotypes, (0, 1, 2)).sum(axis=0)
            heterozygous = (genotypes == 1).sum(axis=0)
            het_values = np.divide(heterozygous, called, out=np.zeros(genotypes.shape[1]), where=called > 0)

            results["mean_heterozygosity"] = float(np.mean(het_values))
            results["heterozygosity_variance"] = float(np.var(het_values))

    if populations and sequences:
        # Population differentiation: the mean over all pairs of populations (one pair gives that pair's value)
        if len(populations) != len(sequences):
            raise ValueError(f"populations has {len(populations)} entries for {len(sequences)} sequences")
        labels = sorted(set(populations), key=str)
        if len(labels) >= 2:
            members = {label: [sequences[j] for j, p in enumerate(populations) if p == label] for label in labels}
            fst_values = [
                calculate_fst(members[first], members[second])
                for index, first in enumerate(labels)
                for second in labels[index + 1 :]
            ]
            results["fst"] = sum(fst_values) / len(fst_values)

    return results


def compare_populations(pop1_data: Dict[str, Any], pop2_data: Dict[str, Any]) -> Dict[str, Any]:
    """Compare two populations based on their summary statistics.

    Args:
        pop1_data: Summary statistics for population 1
        pop2_data: Summary statistics for population 2

    Returns:
        Dictionary containing comparison metrics

    Example:
        >>> pop1 = {"nucleotide_diversity": 0.01, "tajima_d": 0.5}
        >>> pop2 = {"nucleotide_diversity": 0.015, "tajima_d": -0.3}
        >>> comparison = compare_populations(pop1, pop2)
        >>> isinstance(comparison, dict)
        True
    """
    comparison = {}

    # Compare nucleotide diversity
    if "nucleotide_diversity" in pop1_data and "nucleotide_diversity" in pop2_data:
        pi1 = pop1_data["nucleotide_diversity"]
        pi2 = pop2_data["nucleotide_diversity"]
        comparison["nucleotide_diversity_ratio"] = pi2 / pi1 if pi1 > 0 else float("inf")
        comparison["nucleotide_diversity_difference"] = pi2 - pi1

    # Compare Tajima's D
    if "tajima_d" in pop1_data and "tajima_d" in pop2_data:
        d1 = pop1_data["tajima_d"]
        d2 = pop2_data["tajima_d"]
        comparison["tajima_d_difference"] = d2 - d1

        # Classify selection patterns
        def classify_selection(d: float | None) -> str:
            if d is None:
                return "unknown"
            elif d > 2:
                return "balancing_selection"
            elif d < -2:
                return "positive_selection"
            else:
                return "neutral"

        comparison["pop1_selection_pattern"] = classify_selection(d1)
        comparison["pop2_selection_pattern"] = classify_selection(d2)

    # Compare F_ST if available
    if "fst" in pop1_data and "fst" in pop2_data:
        fst1 = pop1_data["fst"]
        fst2 = pop2_data["fst"]
        if fst1 is not None and fst2 is not None:
            comparison["fst_difference"] = fst2 - fst1
            comparison["fst_ratio"] = fst2 / fst1 if fst1 > 0 else float("inf")

    return comparison


def calculate_fay_wu_h(
    sequences: List[str],
    outgroup: Optional[str] = None,
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Tuple[float, Optional[float]]:
    """Calculate Fay and Wu's H and, with an outgroup, its two-sided p-value.

    ``h`` is :func:`metainformant.dna.population.core.fay_wu_h_from_sequences`
    (per-site theta_pi - theta_H). Its null distribution needs derived-allele
    frequencies, so the p-value is that of Zeng et al.'s (2006) normalized H
    (:func:`metainformant.dna.population.core.normalized_fay_wu_h`) and is None when no
    outgroup polarizes the sites of a polymorphic sample. A monomorphic sample has H = 0
    and p-value 1.0. With ``p_value_method="normal"`` (default) the p-value uses the normal
    approximation; with ``"simulation"`` it is the Monte Carlo p-value of
    :func:`calculate_tajima_d` for the outgroup-polarized (unfolded) spectrum, conditioned on the
    number of polarized mutations. Like the normal approximation, the simulated null assumes the
    outgroup polarizes every site correctly.

    Args:
        sequences: List of aligned DNA sequences
        outgroup: Optional outgroup sequence aligned to ``sequences``
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (integer, ``random.Random`` or
            ``numpy.random.Generator``); ``None`` draws fresh entropy

    Returns:
        Tuple of (fay_wu_h, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Raises:
        ValueError: If ``p_value_method`` is unknown, ``n_simulations`` is not a positive integer, or
            the sequences or outgroup are not aligned

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> h, p = calculate_fay_wu_h(seqs)
        >>> isinstance(h, float), p is None
        (True, True)
    """
    _check_p_value_options(p_value_method, n_simulations)
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    h = population.fay_wu_h_from_sequences(sequences, outgroup)
    if outgroup is not None:
        normalized_h = population.normalized_fay_wu_h(sequences, outgroup)  # also validates the outgroup
        if p_value_method == "simulation":
            from ..integration.neutrality_null import simulated_neutrality_tests

            simulated = simulated_neutrality_tests(
                sequences, ("fay_wu_h",), n_simulations=n_simulations, random_state=random_state, outgroup=outgroup
            )
            return h, simulated["fay_wu_h"][1]
        return h, _two_sided_normal_p(normalized_h)
    if population.segregating_sites(sequences) == 0:
        return h, 1.0
    return h, None


def calculate_fu_li_f(
    sequences: List[str],
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Tuple[float, float]:
    """Calculate Fu and Li's F* (no outgroup) and its two-sided p-value.

    See :func:`metainformant.dna.population.core.fu_and_li_f_star_from_sequences`. The p-value is
    obtained as in :func:`calculate_tajima_d`: by the normal approximation (default) or, with
    ``p_value_method="simulation"``, against ``n_simulations`` neutral coalescent samples conditioned on
    the observed number of mutations.

    Args:
        sequences: List of aligned DNA sequences
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (integer, ``random.Random`` or
            ``numpy.random.Generator``); ``None`` draws fresh entropy

    Returns:
        Tuple of (fu_li_f_star, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Raises:
        ValueError: If ``p_value_method`` is unknown, ``n_simulations`` is not a positive integer, or
            the sequences are not aligned

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> f_star, p = calculate_fu_li_f(seqs)
        >>> isinstance(f_star, float)
        True
    """
    _check_p_value_options(p_value_method, n_simulations)
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    f_star = population.fu_and_li_f_star_from_sequences(sequences)
    return f_star, _neutrality_p_value(sequences, "fu_li_f_star", f_star, p_value_method, n_simulations, random_state)


def neutrality_test_suite(
    sequences: List[str],
    outgroup: Optional[str] = None,
    *,
    p_value_method: str = "normal",
    n_simulations: int = DEFAULT_NULL_SIMULATIONS,
    random_state: int | random.Random | np.random.Generator | None = None,
) -> Dict[str, Any]:
    """Run a comprehensive suite of neutrality tests.

    ``p_value_method`` is passed to every test of the suite, which obtains its p-value by the normal
    approximation (default) or by a simulated neutral coalescent null, as described in
    :func:`calculate_tajima_d`. With an integer ``random_state`` every test sees the same simulated
    replicates (they depend only on the sample size and the number of mutations), so the p-values
    are reproducible and comparable across tests.

    Args:
        sequences: List of DNA sequences from a population
        outgroup: Optional outgroup sequence; without it ``fay_wu_h_p_value`` is None
        p_value_method: ``"normal"`` or ``"simulation"``
        n_simulations: Number of coalescent replicates per test for ``"simulation"`` (default 2000)
        random_state: Seed of the simulated replicates (integer, ``random.Random`` or
            ``numpy.random.Generator``); ``None`` draws fresh entropy

    Returns:
        Dictionary containing results from multiple neutrality tests. With ``p_value_method`` other
        than ``"normal"`` it also records ``p_value_method`` and ``n_simulations``.

    Raises:
        ValueError: If ``p_value_method`` is unknown or ``n_simulations`` is not a positive integer

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> results = neutrality_test_suite(seqs)
        >>> isinstance(results, dict)
        True
    """
    _check_p_value_options(p_value_method, n_simulations)
    results: Dict[str, Any] = {}

    if not sequences:
        return results

    # Basic statistics
    results["num_sequences"] = len(sequences)
    results["sequence_length"] = len(sequences[0]) if sequences else 0
    if p_value_method != "normal":
        results["p_value_method"] = p_value_method
        results["n_simulations"] = n_simulations

    p_options: Dict[str, Any] = {
        "p_value_method": p_value_method,
        "n_simulations": n_simulations,
        "random_state": random_state,
    }

    # Tajima's D
    try:
        tajima_d, tajima_p = calculate_tajima_d(sequences, **p_options)
        results["tajima_d"] = tajima_d
        results["tajima_d_p_value"] = tajima_p
    except Exception as e:
        logger.warning(f"Failed to calculate Tajima's D: {e}")
        results["tajima_d"] = None
        results["tajima_d_p_value"] = None

    # Fu and Li's D*
    try:
        fu_li_d_star, fu_li_d_star_p = calculate_fu_li_d(sequences, **p_options)
        results["fu_li_d_star"] = fu_li_d_star
        results["fu_li_d_star_p_value"] = fu_li_d_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's D*: {e}")
        results["fu_li_d_star"] = None
        results["fu_li_d_star_p_value"] = None

    # Fu and Li's F*
    try:
        fu_li_f_star, fu_li_f_star_p = calculate_fu_li_f(sequences, **p_options)
        results["fu_li_f_star"] = fu_li_f_star
        results["fu_li_f_star_p_value"] = fu_li_f_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's F*: {e}")
        results["fu_li_f_star"] = None
        results["fu_li_f_star_p_value"] = None

    # Fay and Wu's H
    try:
        fay_wu_h, fay_wu_h_p = calculate_fay_wu_h(sequences, outgroup, **p_options)
        results["fay_wu_h"] = fay_wu_h
        results["fay_wu_h_p_value"] = fay_wu_h_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fay and Wu's H: {e}")
        results["fay_wu_h"] = None
        results["fay_wu_h_p_value"] = None

    # Summary interpretation
    results["neutrality_summary"] = interpret_neutrality_results(results)

    return results


_NEUTRALITY_LABELS = {
    "tajima_d": ("balancing_selection_or_population_bottleneck", "positive_selection_or_population_expansion"),
    "fu_li_d_star": ("balancing_selection", "positive_selection_or_population_expansion"),
}


def interpret_neutrality_results(results: Dict[str, Any], alpha: float = 0.05) -> Dict[str, str]:
    """Interpret neutrality test results.

    A statistic gets a directional label only when it departs from neutrality:
    its ``<name>_p_value`` is below ``alpha``, or the results carry no p-value
    and the statistic is nonzero. Otherwise it is ``"neutral_evolution"``,
    meaning consistent with neutrality.

    Args:
        results: Results from neutrality_test_suite
        alpha: Significance level for the p-values

    Returns:
        Dictionary with interpretation of results
    """
    interpretation = {}
    for name, (positive, negative) in _NEUTRALITY_LABELS.items():
        value = results.get(name)
        p_value = results.get(f"{name}_p_value")
        if value is None:
            interpretation[name] = "calculation_failed"
        elif (p_value is not None and p_value >= alpha) or value == 0:
            interpretation[name] = "neutral_evolution"
        else:
            interpretation[name] = positive if value > 0 else negative
    return interpretation
