"""DNA mutation detection and analysis utilities.

This module provides tools for detecting, classifying, and analyzing mutations
in DNA sequences, including point mutations, insertions, deletions, and
evolutionary rate calculations.
"""

from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def calculate_mutation_rate(ancestral: str, derived: str) -> float:
    """Calculate mutation rate between ancestral and derived sequences.

    Args:
        ancestral: Ancestral DNA sequence
        derived: Derived DNA sequence

    Returns:
        Mutation rate (mutations per site)

    Raises:
        ValueError: If sequences have different lengths

    Example:
        >>> ancestral = "ATCG"
        >>> derived = "ATCG"
        >>> rate = calculate_mutation_rate(ancestral, derived)
        >>> rate == 0.0
        True
    """
    if len(ancestral) != len(derived):
        raise ValueError("Sequences must have equal length")

    if not ancestral:
        return 0.0

    mutations = 0
    for a, d in zip(ancestral.upper(), derived.upper()):
        if a != d:
            mutations += 1

    return mutations / len(ancestral)


def classify_mutations(ancestral: str, derived: str) -> Dict[str, int]:
    """Classify mutations by type between ancestral and derived sequences.

    Args:
        ancestral: Ancestral DNA sequence
        derived: Derived DNA sequence

    Returns:
        Dictionary with the counts ``transitions``, ``transversions``, ``synonymous``,
        ``nonsynonymous``, ``ambiguous`` and ``total``. ``total`` counts every position at which the
        sequences differ. A difference that involves a base other than A, C, G or T (``N``, a gap,
        an IUPAC code) has no transition or transversion type and is counted under ``ambiguous``
        only, so ``transitions + transversions + ambiguous == total``. ``synonymous`` and
        ``nonsynonymous`` read the sequence as coding from position 0 and judge each difference
        alone, as a single-base change of the ancestral codon; positions in an incomplete final
        codon, or in a codon that holds a non-ACGT base, are not classified, so the two can sum to
        less than ``transitions + transversions``.

    Raises:
        ValueError: If sequences have different lengths

    Example:
        >>> ancestral = "ATCG"
        >>> derived = "AGCT"
        >>> mutations = classify_mutations(ancestral, derived)
        >>> "transitions" in mutations
        True
    """
    if len(ancestral) != len(derived):
        raise ValueError("Sequences must have equal length")

    mutation_types = {
        "transitions": 0,
        "transversions": 0,
        "synonymous": 0,
        "nonsynonymous": 0,
        "ambiguous": 0,
        "total": 0,
    }

    # Transition pairs: A<->G, C<->T
    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    # Codon translation table (simplified)
    genetic_code = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    for i, (a, d) in enumerate(zip(ancestral.upper(), derived.upper())):
        if a != d:
            mutation_types["total"] += 1

            # A base outside A/C/G/T has no transition or transversion type
            if a not in "ACGT" or d not in "ACGT":
                mutation_types["ambiguous"] += 1
                continue

            # Classify transition vs transversion
            if (a, d) in transition_pairs:
                mutation_types["transitions"] += 1
            else:
                mutation_types["transversions"] += 1

            # Check if this affects protein coding (simplified)
            # Look at codon context
            codon_start = (i // 3) * 3
            if codon_start + 3 <= len(ancestral):
                ancestral_codon = ancestral[codon_start : codon_start + 3].upper()
                derived_codon = ancestral_codon[: i % 3] + d + ancestral_codon[i % 3 + 1 :]

                ancestral_aa = genetic_code.get(ancestral_codon)
                derived_aa = genetic_code.get(derived_codon)
                if ancestral_aa is None or derived_aa is None:
                    continue  # a non-ACGT base elsewhere in the codon: the amino acid is unknown

                if ancestral_aa == derived_aa:
                    mutation_types["synonymous"] += 1
                else:
                    mutation_types["nonsynonymous"] += 1

    return mutation_types


def generate_point_mutations(
    sequence: str,
    num_mutations: int | None = None,
    mutation_rate: float = 0.001,
    *,
    rng: random.Random | None = None,
) -> str:
    """Generate point mutations in a DNA sequence.

    ``num_mutations`` selects between two modes. With ``num_mutations`` given, exactly that many
    distinct positions are sampled uniformly at random and ``mutation_rate`` is not used. Without
    it, each position mutates independently with probability ``mutation_rate``, so the number of
    mutations is binomial. Every mutated position changes to a different nucleotide.

    Args:
        sequence: Original DNA sequence
        num_mutations: Exact number of distinct positions to mutate, or None to draw the
            mutations from ``mutation_rate``
        mutation_rate: Probability that each position mutates, within [0, 1]; applies only when
            ``num_mutations`` is None
        rng: Random generator for reproducibility. The default draws from Python's global
            ``random`` state, which ``random.seed`` controls

    Returns:
        Mutated DNA sequence (upper case)

    Raises:
        ValueError: If num_mutations is negative or exceeds the sequence length, or, when
            num_mutations is None, mutation_rate is outside [0, 1]

    Example:
        >>> seq = "ATCGATCG"
        >>> mutated = generate_point_mutations(seq, num_mutations=2)
        >>> len(mutated) == len(seq)
        True
    """
    if not sequence:
        return sequence

    mutated = list(sequence.upper())
    source: Any = random if rng is None else rng

    if num_mutations is None:
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError("mutation_rate must be within [0, 1]")
        positions = [pos for pos in range(len(mutated)) if source.random() < mutation_rate]
    else:
        if num_mutations < 0:
            raise ValueError("num_mutations must be non-negative")
        if num_mutations > len(mutated):
            raise ValueError("num_mutations cannot exceed the sequence length")
        positions = source.sample(range(len(mutated)), num_mutations)

    nucleotides = ["A", "C", "G", "T"]
    for pos in positions:
        current = mutated[pos]
        # Choose different nucleotide
        possible = [n for n in nucleotides if n != current]
        mutated[pos] = source.choice(possible)

    return "".join(mutated)


def _binomial_count(trials: int, probability: float, source: Any) -> int:
    """Draw a Binomial(``trials``, ``probability``) count from ``source.random()``.

    The count is the number of successes among ``trials`` independent Bernoulli trials, found by
    jumping from one success to the next: the number of failures before a success is geometric,
    ``floor(ln(1 - U) / ln(1 - p))`` for ``U`` uniform on [0, 1). The draw is exact and costs about
    ``trials * probability + 1`` random numbers rather than ``trials``, and it uses only ``random()``,
    so a seeded generator gives the same count on every Python version (``Random.binomialvariate``
    exists only from Python 3.12).

    Args:
        trials: Number of Bernoulli trials, non-negative
        probability: Success probability of each trial, within [0, 1]
        source: Object with a ``random()`` method, a ``random.Random`` or the ``random`` module

    Returns:
        Number of successes, between 0 and ``trials``
    """
    if trials <= 0 or probability <= 0.0:
        return 0
    if probability >= 1.0:
        return trials

    log_failure = math.log1p(-probability)  # ln(1 - p), negative
    count = 0
    last_success = -1  # index of the most recent success; trials are indexed 0 .. trials - 1
    while True:
        failures_before_next = math.log1p(-source.random()) / log_failure
        if failures_before_next >= trials - last_success - 1:
            return count  # the next success would fall beyond the last trial
        last_success += int(failures_before_next) + 1
        count += 1


def simulate_sequence_evolution(
    sequence: str,
    generations: int,
    mutation_rate: float = 0.001,
    *,
    rng: random.Random | None = None,
) -> List[str]:
    """Simulate sequence evolution under mutation pressure.

    Every generation the number of point mutations is drawn from Binomial(L, ``mutation_rate``),
    with L the sequence length, so each site mutates with probability ``mutation_rate`` per
    generation and a lineage accumulates ``mutation_rate * L * generations`` substitutions on
    average. That many distinct positions, chosen uniformly, then change to a nucleotide different
    from the current one (as in :func:`generate_point_mutations`). Mutations are not reverted or
    excluded: a position may mutate again in a later generation, and sequence length is fixed
    (no insertions or deletions).

    Args:
        sequence: Starting DNA sequence
        generations: Number of generations to simulate
        mutation_rate: Probability that a site mutates in one generation, within [0, 1]
        rng: Random generator that supplies the mutation counts, positions and replacement bases.
            The default draws from Python's global ``random`` state, which ``random.seed`` controls;
            a private ``random.Random`` makes the lineage reproducible without touching it

    Returns:
        List of sequences, one per generation. The first is ``sequence`` as given and the rest are
        upper case; consecutive sequences differ at exactly the number of positions drawn for that
        generation. An empty sequence or a negative ``generations`` gives ``[sequence]``.

    Raises:
        ValueError: If ``mutation_rate`` is outside [0, 1] (NaN included)

    Example:
        >>> seq = "ATCG"
        >>> evolution = simulate_sequence_evolution(seq, generations=3, rng=random.Random(1))
        >>> len(evolution) == 4  # Including original
        True
    """
    if not 0.0 <= mutation_rate <= 1.0:
        raise ValueError("mutation_rate must be within [0, 1]")
    if not sequence or generations < 0:
        return [sequence]

    source: Any = random if rng is None else rng
    evolutionary_lineage = [sequence]

    for _ in range(generations):
        current_seq = evolutionary_lineage[-1]
        n_mutations = _binomial_count(len(current_seq), mutation_rate, source)
        evolutionary_lineage.append(generate_point_mutations(current_seq, num_mutations=n_mutations, rng=rng))

    return evolutionary_lineage


def find_mutational_hotspots(sequence: str, window_size: int = 10) -> List[Tuple[int, float]]:
    """Score sliding windows by a simple GC-content and base-repeat heuristic.

    The score of a window is ``0.6 * GC fraction + 0.4 * (adjacent identical bases / window size)``.
    The weights are a fixed heuristic, not an empirical mutability model (no CpG, context-dependent
    or polymerase-slippage rates are fitted), so a high score marks a GC-rich or repetitive window,
    not a measured excess of mutations.

    Args:
        sequence: DNA sequence to analyze
        window_size: Size of sliding window for analysis

    Returns:
        List of tuples (window start position, heuristic score in [0, 1]), one per window; an empty
        list for an empty sequence, a non-positive window or a window longer than the sequence

    Example:
        >>> seq = "ATCGATCGATCG"
        >>> hotspots = find_mutational_hotspots(seq, window_size=4)
        >>> isinstance(hotspots, list)
        True
    """
    if not sequence or window_size <= 0:
        return []

    hotspots = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size].upper()

        # Calculate GC content (simple hotspot predictor)
        gc_count = window.count("G") + window.count("C")
        gc_content = gc_count / window_size

        # Calculate repeat content
        repeat_score = 0
        for j in range(window_size - 1):
            if window[j] == window[j + 1]:
                repeat_score += 1

        # Combined hotspot score
        hotspot_score = (gc_content * 0.6) + (repeat_score / window_size * 0.4)

        hotspots.append((i, hotspot_score))

    return hotspots


def calculate_substitution_matrix(seq1: str, seq2: str) -> Dict[Tuple[str, str], int]:
    """Calculate nucleotide substitution matrix between two sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Dictionary mapping (from_nuc, to_nuc) pairs to counts

    Raises:
        ValueError: If sequences have different lengths

    Example:
        >>> seq1 = "ATCG"
        >>> seq2 = "AGCT"
        >>> matrix = calculate_substitution_matrix(seq1, seq2)
        >>> ("T", "G") in matrix
        True
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")

    substitutions: Dict[Tuple[str, str], int] = {}

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a != b:
            pair = (a, b)
            substitutions[pair] = substitutions.get(pair, 0) + 1

    return substitutions


def detect_selection_signatures(sequence: str, reference_sequences: List[str]) -> Dict[str, float]:
    """Summarize how far a sequence has diverged from reference sequences.

    ``diversity`` is the mean proportion of differing sites between ``sequence`` and each reference
    of the same length (:func:`calculate_mutation_rate`; a reference of another length cannot be
    compared site by site and is ignored) and ``conservation`` is ``1 - diversity``. These are
    divergence measures, not tests of selection: nothing here separates synonymous from
    nonsynonymous change (see :func:`classify_mutations`), models drift or estimates dN/dS.

    Args:
        sequence: Query sequence
        reference_sequences: List of reference sequences

    Returns:
        Dictionary with ``diversity`` and ``conservation``

    Raises:
        ValueError: If no reference has the same length as ``sequence`` (an empty list included),
            because the divergence is then undefined, not zero

    Example:
        >>> query = "ATCGATCG"
        >>> refs = ["ATCGATCG", "ATCGATCG"]
        >>> signatures = detect_selection_signatures(query, refs)
        >>> "diversity" in signatures
        True
    """
    diversities = [
        calculate_mutation_rate(sequence, reference)
        for reference in reference_sequences
        if len(reference) == len(sequence)
    ]
    if not diversities:
        raise ValueError("reference_sequences must contain at least one sequence of the same length as the query")

    diversity = sum(diversities) / len(diversities)
    return {"diversity": diversity, "conservation": 1.0 - diversity}


def generate_mutant_library(sequence: str, positions: List[int], mutations: List[str]) -> List[str]:
    """Generate a library of mutant sequences.

    Args:
        sequence: Original DNA sequence
        positions: List of positions to mutate
        mutations: List of mutation strings (one per position)

    Returns:
        List of mutant sequences

    Example:
        >>> seq = "ATCG"
        >>> mutants = generate_mutant_library(seq, [1, 3], ["G", "T"])
        >>> len(mutants) == 2
        True
    """
    if len(positions) != len(mutations):
        raise ValueError("Positions and mutations lists must have equal length")

    mutants = []

    for pos, mutation in zip(positions, mutations):
        if 0 <= pos < len(sequence):
            mutant = list(sequence)
            mutant[pos] = mutation.upper()
            mutants.append("".join(mutant))

    return mutants


def analyze_mutation_spectrum(sequence: str, reference: str) -> Dict[str, Any]:
    """Analyze the spectrum of mutations between sequences.

    Args:
        sequence: Derived sequence
        reference: Reference sequence

    Returns:
        Dictionary with ``mutation_types`` (:func:`classify_mutations`), ``mutation_density`` (differing
        positions per site), ``ti_tv_ratio``, ``sequence_length`` and ``substitution_matrix``
        (:func:`calculate_substitution_matrix`). ``ti_tv_ratio`` is ``inf`` when there are transitions
        but no transversions and ``nan`` when there are neither, because the ratio is then undefined.

    Example:
        >>> ref = "ATCGATCG"
        >>> seq = "AGCGATCG"
        >>> spectrum = analyze_mutation_spectrum(seq, ref)
        >>> "mutation_types" in spectrum
        True
    """
    if len(sequence) != len(reference):
        raise ValueError("Sequences must have equal length")

    mutation_types = classify_mutations(reference, sequence)

    # Calculate mutation density
    mutation_density = mutation_types["total"] / len(sequence) if sequence else 0

    # Calculate transition/transversion ratio; 0 / 0 is undefined, not infinite
    if mutation_types["transversions"] > 0:
        ti_tv_ratio = mutation_types["transitions"] / mutation_types["transversions"]
    elif mutation_types["transitions"] > 0:
        ti_tv_ratio = float("inf")
    else:
        ti_tv_ratio = float("nan")

    spectrum = {
        "mutation_types": mutation_types,
        "mutation_density": mutation_density,
        "ti_tv_ratio": ti_tv_ratio,
        "sequence_length": len(sequence),
        "substitution_matrix": calculate_substitution_matrix(reference, sequence),
    }

    return spectrum
